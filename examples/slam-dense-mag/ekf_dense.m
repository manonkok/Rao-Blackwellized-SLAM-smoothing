function [xf_traj,qnb_traj,Pf_traj] = ekf_dense(dynModel,measModel,odometry,y,...
            x0,q0,P0,Q,R,dt)
%% EKF_DENSE - Run EKF for comparison for the conditionally linear problem
%
% Syntax:
%   [xf_traj,qnb_traj,Pf_traj] = ekf_dense(dynModel,measModel,odometry,y,...
%            x0,P0,Q,R,dt)
%
% In:
%   dynModel    - Dynamical model function handle, xn = @(xn,dx,dt,Q)
%   measModel   - Measurement model function handle, y = @(xn,xl) + r
%   odometry    - Dynamical model control observations [N_T x n_o]
%   y           - Observations [N_T x n_y]
%   x0          - Initial state [nStates x 1]
%   P0          - Initial covariance [nStates x nStates]
%   Q           - Process noise cov [nw x nw x N_T] or [nw x nw]
%   R           - Measurement noise cov [n_y x n_y]
%   dt          - Time between two time steps [N_T x 1] or scalar

% Out:
%   xf_traj     - Filtered mean estimates
%   qnb_traj    - Filtered mean orientation estimates
%   PF_traj     - Filtered covariance matrices
%
% Description:
%   Run EKF for dense, conditionally linear state space model.
%
% References:
%   [1] Frida Viset, Rudy Helmons and Manon Kok, An Extended Kalman Filter 
%       for Magnetic Field SLAM Using Gaussian Process Regression, 
%       MDPI Sensors, 22(8), 2833, 2022.
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

%% Initialise states and covariance matrices
xp = x0;
Pp = P0;
q_nb = q0;

%% Parameters and settings
% Extract some parameters
nStates = size(x0,1);
N_T = size(y,1);
  
% Allow for both time-varying and constant Q
if size(Q,3) == 1 
    Q = repmat(Q,[1 1 N_T-1]);
end
  
% Allow for both time-varying and constant time step
if length(dt) == 1 
    dt = dt * ones(N_T-1,1);
end

% Jitter to use if Cholesky decomposition fails due to numerical instability
jitter = 1e-3; 

iOri = 4:6;
  
%% Preallocate trajectories
xf_traj = nan(nStates,N_T); % Filtered mean estimates
Pf_traj = nan(nStates,nStates,N_T); % Filtered covariance matrices
qnb_traj = nan(4,N_T); % Filtered mean orientation estimates   

%% Filter recursion
for t=1:N_T
    % EKF prediction
    if t ~= 1 % Don't do a prediction at the very first time instance
        % Propagate state through dynamics
        [xp,q_nb,F,G] = dynModel(xf,q_nb,odometry(t-1,:)); 
        Qt = dt(t-1) * Q(:,:,t-1); 
        Pp = F * Pf * F' + G * Qt * G';
    end
    
    % Compute innovations and their covariances
    yt = y(t,:); % Measurements at time t
    [yhat,dy] = measModel(xp,q_nb); % Compute gradient of cond linear measurement model
    e = yt' - yhat;
    SS = dy*Pp*dy' + R;
        
    % Compute Kalman gain
    [cS,flag] = chol(SS,'lower');
    if flag>0
        cS = chol(SS+jitter*eye(size(SS,1)),'lower');
    end
    K = Pp*((dy'/cS')/cS);
        
    % Update states and covariances
    xf  = xp + K*e;  
    Pf = Pp - K*SS*K';
    Pf = 0.5 * (Pf + Pf'); % To keep the covariance matrices symmetric

    % Relinearise orientation
    q_nb = qLeft(expq(xf(iOri)/2)) * q_nb;
    xf(iOri) = zeros(3,1);

    % Save state, covariance and linearisation point orientation
    xf_traj(:,t) = xf;
    Pf_traj(:,:,t) = Pf;
    qnb_traj(:,t) = q_nb;
end

end