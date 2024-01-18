function [traj_max,traj_mean,xl_max,xl_mean,P_max,P_mean,traj_sample_iwmax,xn_traj] = ...
    particleFilter(dynModel,measModel,odometry,y,...
    x0_nonLin,x0_lin,P0_lin,Q,R,N_P,dt,sparseFeatures,makePlots)
% PARTICLEFILTER - Run Rao-Blackwellized particle filter
%
% Syntax:
%   [traj_max,traj_mean,xl_max,xl_mean,P_max,P_mean,traj_sample,xn_traj] = 
%       particleFilter(dynModel,measModel,odometry,y,...
%           x0_nonLin,x0_lin,P0_lin,Q,R,N_P,dt,sparseFeatures,makePlots)
%
% In:
%   dynModel    - Dynamical model function handle, xn = @(xn,dx,dt,Q)
%   measModel   - Measurement model function handle, y = @(xn,xl) + r
%   odometry    - Dynamical model control observations [N_T x n_o]
%   y           - Observations [N_T x n_y]
%   x0_nonLin   - Initial non-linear state [nNonLin x 1]
%   x0_lin      - Initial linear states [nLin x N_P] or [nLin x 1]
%   P0_lin      - Initial linear state cov [nLin x nLin]
%   Q           - Process noise cov [nw x nw x N_T] or [nw x nw]
%   R           - Measurement noise cov [n_y x n_y]
%   N_P         - Number of particles
%   dt          - Time between two time steps [N_T x 1] or scalar
%   sparseFeatures - Use an EKF observation model
%   makePlots   - Visualization function (default: none)
%
% Out:
%   traj_max    - Highest-weight trajectory (nonlinear)
%   traj_mean   - Weighted-mean trajectory (nonlinear) 
%   xl_max      - Maximum-weight linear state (~map)
%   xl_mean     - Weighted-mean linear state (~map)
%   P_max       - Covariance of maximum-weight linear state
%   P_mean      - Covariance of weighted-mean linear state
%   traj_sample_iwmax - Trajectory of particle with highest weight at t = N_T
%   xn_traj     - Trajectories of nonlinear states over time
%
% Description:
%   Run particle filter for a conditionally linear or conditionally
%   linearized state space model. See [1] for details.
%
% References:
%
%   [1] Manon Kok, Arno Solin, and Thomas B. Schon. Rao-Blackwellized 
%       Particle Smoothing for Simultaneous Localization and Mapping.
%       pre-print: https://arxiv.org/abs/2306.03953
%
% See also:
%   particleSmoother
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

%% Initialise weights, states and covariance matrices

% Initial weights
w = 1/N_P * ones(1,N_P);
logw = log(w);
  
% Initial states (nonlinear) and initial means (linear states)
xn = repmat(x0_nonLin,1,N_P); % Nonlinear states
if size(x0_lin,2) > 1 % Conditionally linear Gaussian states
    xl = x0_lin;
else
    xl = repmat(x0_lin,1,N_P); 
end
  
% Initial covariance matrices for linear states
P = repmat(P0_lin,[1,1,N_P]);

%% Parameters and settings
% Extract some parameters
nNonLin = size(x0_nonLin,1);
N_T = size(y,1);
  
% Allow for both time-varying and constant Q
if size(Q,3) == 1 
    Q = repmat(Q,[1 1 N_T-1]);
end
  
% Allow for both time-varying and constant time step
if length(dt) == 1 
    dt = dt * ones(N_T-1,1);
end
  
% Defaults
if nargin < 11 || isempty(sparseFeatures), sparseFeatures = false; end
if nargin < 12 || isempty(makePlots), makePlots = []; end

% Jitter to use if Cholesky decomposition fails due to numerical instability
jitter = 1e-3; 
  
%% Preallocate trajectories
traj_max = nan(nNonLin,N_T); % Maximum-weight trajectory
traj_mean = nan(nNonLin,N_T); % Weighted-mean trajectory
yhattraj = nan(size(y,2), N_T); % Predicted measurement by the maximum-weight particle
xn_traj = zeros(nNonLin, N_P, N_T); % Collection of all trajectories
xn_traj(:,:,1) = xn;
ai = zeros(N_P,1); % Sampled ancestors
    
%% Filter recursion
for t=1:N_T
    % Particle filter prediction
    xn_ = xn; % Copy old nonlinear states
    if t ~= 1 % Don't do a prediction at the very first time instance
        for i = 1:N_P
            % Draw ancestor index...
            ai(i) = sample(w); 
            % ... and propagate that nonlinear state through dynamics
            xn(:,i) = dynModel(xn_(:,ai(i)),odometry(t-1,:),dt(t-1),Q(:,:,t-1)); 
        end
        % The linear states remain equal since no dynamics but shuffle /
        % duplicate based on sampled ancestors
        xl = xl(:,ai);
        P = P(:,:,ai);
        
        % Save trajectory with shuffled ancestor indices, e.g. to visualise
        % path degeneracy
        xn_traj(:,:,t) = xn;
        xn_traj(:,:,1:t-1) = xn_traj(:,ai,1:t-1);
    end
    
    % Compute the importance weights
    yt = y(t,:); % Measurements at time t
    if ~sparseFeatures
        dy = measModel(xn); % Compute gradient of cond linear measurement model
    end
    for i=1:N_P 
        if sparseFeatures
            % Linearize measurement model
            [yhat,dy] = measModel(xn(:,i), xl(:,i));
            % Compute innovations and their covariances
            e = yt' - yhat;
            SS = dy*P(:,:,i)*dy' + R;
            % Strip away those that are not observed
            ind = ~isnan(yt);
            e = e(ind);
            SS = SS(ind,ind);
        else
            % Compute innovations and their covariances
            dyt = squeeze(dy(i,:,:));
            e = yt' - dyt * xl(:,i);
            SS = dyt * P(:,:,i) * dyt' + R;
        end

        % Compute the log weights
        [cS,flag] = chol(SS,'lower');
        if flag>0
            cS = chol(SS+jitter*eye(size(SS,1)),'lower');
        end
        v = cS\e;
        logw(i) = -sum(log(diag(cS))) - .5*(v'*v) - .5*numel(e)*log(2*pi);
    end

    % Normalize by log-sum-exp trick
    c = max(logw);
    lse = c + log(sum(exp(logw - c)));
    w = exp(logw - lse);  
    
    % Store trajectories
    [~,iw_max] = max(w);
    traj_max(:,t) = xn(:,iw_max);   % Store maximum-weight particle
    traj_mean(:,t) = sum(xn.*w,2);  % Store weighted-mean particle

    % Kalman filter measurement update
    for i = 1:N_P
        if sparseFeatures
            % Linearize measurement
            [yhat,dy] = measModel(xn(:,i), xl(:,i));
            % Compute innovations and their covariances
            e = yt' - yhat;
            SS = dy*P(:,:,i)*dy' + R;
            % Strip away those that are not observed
            ind = ~isnan(yt);
            e = e(ind);
            SS = SS(ind,ind);
            % Compute Kalman gain
            [cS,flag] = chol(SS,'lower');
            if flag>0
                cS = chol(SS+jitter*eye(size(SS,1)),'lower');
            end
            K = P(:,:,i)*((dy(ind,:)'/cS')/cS);
        else
            % Use the dy computed previously to compute predicted
            % measurement
            dyt = squeeze(dy(i,:,:)); 
            yhat = dyt * xl(:,i);
            % Compute innovations and their covariances
            e = yt' - yhat;    
            SS = dyt * P(:,:,i) * dyt' + R;    
            % Compute Kalman gain
            [cS,flag] = chol(SS,'lower');
            if flag>0
                cS = chol(SS+jitter*eye(size(SS,1)),'lower');
            end
            K = P(:,:,i)*((dyt'/cS')/cS);
        end
        % Update linear states and covariances
        xl(:,i)  = xl(:,i) + K*e;  
        P(:,:,i) = P(:,:,i) - K*SS*K'; 
        
        % Predicted measurements of highest weight particle, store for plotting
        if i == iw_max
            yhattraj(:,t) = yhat;
        end
    end
    
    % Plot results using
    % 1) all particles for the nonlinear states
    % 2) the highest weight map and its covariance
    % 3) the sequence of maximum weight particles up to this time
    % 4) predicted measurements along the trajectory of the highest
    % weight particles
    % 5) all trajectories
    % 6) the sequence of weighted mean particles up to this time
    % 7) all linear states and their covariances
    if ~isempty(makePlots)
      makePlots(xn,xl(:,iw_max),P(:,:,iw_max),traj_max,yhattraj,xn_traj,traj_mean,xl,P)
    end
end

%% Extract final map and trajectory
% Map of highest weight particle and its covariance
xl_max = xl(:,iw_max);
P_max = P(:,:,iw_max);
  
% Weighted mean map and covariance
xl_mean = sum(xl.*w,2);
P_mean = zeros(length(xl_mean));
for i = 1:N_P
    P_mean = w(i) * ( P(:,:,i) + (xl_mean - xl(:,i)) * (xl_mean - xl(:,i))');
end
  
% Trajectory of particle with highest weight at last time instance
traj_sample_iwmax = squeeze(xn_traj(:,iw_max,:));

end