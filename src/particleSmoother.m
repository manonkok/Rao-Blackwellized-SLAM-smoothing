function [XNK,XLK,PK] = particleSmoother(dynModel,measModel,dynResNorm,odometry,y,...
    x0_nonLin,x0_lin,P0_lin,Q,R,N_P,N_K,dt,sparseFeatures,makePlots)
% PARTICLESMOOTHER - Run Rao-Blackwellized particle smoother
%
% Syntax:
%   [XNK,XLK,PK] = particleSmoother(dynModel,measModel,dynRes,odometry,y,...
%    x0_nonLin,x0_lin,P0_lin,Q,R,N_P,N_K,dt,sparseFeatures,makePlots)
%
% In:
%   dynModel    - Dynamical model function handle, xn = @(xn,dx,dt,Q)
%   measModel   - Measurement model function handle, y = @(xn,xl) + r
%   dynResNorm  - Function handle to compute normalised residual of dynamic
%   model, eDyn = @(xnk,xn,odometry,dt,Q);
%   odometry    - Dynamical model control observations [N_T x n_o]
%   y           - Observations [N_T x n_y]
%   x0_nonLin   - Initial non-linear state [nNonLin x 1]
%   x0_lin      - Initial linear states [nLin x N_P] or [nLin x 1]
%   P0_lin      - Initial linear state cov [nLin x nLin]
%   Q           - Process noise cov [nw x nw x N_T] or [nw x nw]
%   R           - Measurement noise cov [n_y x n_y]
%   N_P         - Number of particles
%   N_K         - Number of smoother iterations
%   dt          - Time between two time steps [N_T x 1] or scalar
%   sparseFeatures - Use an EKF observation model
%   makePlots   - Visualization function (default: none)
%
% Out:
%   XNK         - k sampled trajectories of the nonlinear states
%   XLK         - k sampled means of the parameter estimates (maps)
%   PK          - Covariances of the k sampled parameter estimates
%
% Description:
%   Run particle smoother for a conditionally linear or conditionally
%   linearized state space model. See [1] for details.
%
% References:
%
%   [1] Manon Kok, Arno Solin, and Thomas B. Schon. Rao-Blackwellized 
%       Particle Smoothing for Simultaneous Localization and Mapping.
%       pre-print: https://arxiv.org/abs/2306.03953
%
% See also:
%   particleFilter
%
% Copyright:
%   2023-   Manon Kok and Arno Solin
  
%% Parameters and settings
% Extract some parameters
nNonLin = size(x0_nonLin,1);
nLin = size(x0_lin,1);
N_T = size(y,1);
ny = size(y,2);
  
% Allow for both time-varying and constant Q
if size(Q,3) == 1 
    Q = repmat(Q,[1 1 N_T-1]);
end
  
% Allow for both time-varying and constant time step
if length(dt) == 1 
    dt = dt * ones(N_T-1,1);
end
  
% Defaults
if nargin < 14 || isempty(sparseFeatures), sparseFeatures = false; end
if nargin < 15 || isempty(makePlots), makePlots = []; end

% Jitter to use if Cholesky decomposition fails due to numerical instability
jitter = 1e-2; 
 
%% Preallocate trajectories
xn_traj = zeros(nNonLin, N_P, N_T);
% xl_traj = zeros(nLin, N_P, N_T);
    
ai = zeros(N_P,1);
AI = nan(N_P,N_T);
  
% Allocate space for outputs
XNK = nan(nNonLin,N_T,N_K);
XLK = nan(nLin,N_K);
PK = nan(nLin,nLin,N_K);

%% Run N_K iterations of the conditional Rao-Blackwellised particle filter 
% with ancestor sampling for SLAM. Note that the first iteration is a simple 
% particle filter, in the later iterations we take a reference trajectory 
% of the previous iteration into account.
for k = 1:N_K 
    % Initialise weights, states and covariance matrix (Steps 1-4 of Alg. 2)
    % Initial states (nonlinear)
    xn = repmat(x0_nonLin,1,N_P); % Initialise nonlinear states with x0
    if k~=1
        % If it's not the first smoother iteration, overwrite the N_P'th
        % particle with the state from the sampled trajectory
        xn(:,N_P) = xnk(:,1);
    end
    
    % Initial means (linear states)
    if size(x0_lin,2) > 1 % Conditionally linear Gaussian states
        xl = x0_lin;
    else
        xl = repmat(x0_lin,1,N_P);
    end
    P  = repmat(P0_lin,[1,1,N_P]); % Initial covariance matrix

    % Initial weights
    w = 1/N_P * ones(1,N_P); % weights
    logw = log(w);
    
    % Add reference state to trajectory
    if k~=1
        xn_traj(:,N_P,:) = reshape(xnk,[nNonLin,1,N_T]);
    end

    % Store states at time instance 1 for all particles
    xn_traj(:,:,1) = xn;

    % Pre-compute Ck for the sampled trajectory (Step 5 of Alg. 2)
    if k~=1 && ~sparseFeatures
        dy_xnk = measModel(xnk);
    end
    
    % Filter recursion
    for t=1:N_T
        % Resampling, time update, compute ancestor weights reference
        % trajectory (Steps 7-11 of Alg. 2)
        if t ~= 1 % Don't do a prediction at the very first time instance
            % Particle filter prediction
            xn_pred = zeros(size(xn)); % Pre-allocate predicted nonlinear states
            xl_pred = zeros(size(xl)); % Pre-allocate predicted linear state mean
            P_pred = zeros(size(P)); % Pre-allocate predicted linear state covariance
            for i = 1:N_P-1 % For all particles except the N_P'th one ...
                % ... draw ancestor index...
                ai(i) = sample(w);
                % ... and propagate that nonlinear state through dynamics
                xn_pred(:,i) = dynModel(xn(:,ai(i)),odometry(t-1,:),dt(t-1),Q(:,:,t-1)); 
            end
            % The linear states remain equal since no dynamics but shuffle /
            % duplicate based on sampled ancestors
            xl_pred(:,1:end-1) = xl(:,ai(1:end-1));
            P_pred(:,:,1:end-1) = P(:,:,ai(1:end-1));
            
            % Draw an ancestor sample of the N_P'th particle (Step 10 of
            % Alg. 2)
            if k==1 % If it is the first iteration, then we don't have a 
                % reference trajectory, so just do the same for the N_P'th
                % particle as for the other ones. 
                % Draw ancestor index...
                ai(N_P) = sample(w); 
                % ... and propagate through dynamics
                xn_pred(:,N_P) = dynModel(xn(:,ai(N_P)),odometry(t-1,:),dt(t-1),Q(:,:,t-1)); 
                % The linear state remains equal since no dynamics but shuffle /
                % duplicate based on sampled ancestor
                xl_pred(:,N_P) = xl(:,ai(N_P));
                P_pred(:,:,N_P) = P(:,:,ai(N_P)); 
            else % If this is not the first smoother iteration ... 
                % Vector with log ancestor probabilities of reference
                % trajectory
                paNtLog = zeros(N_P,1);
                
                % Precompute measurement model for trajectory
                if ~sparseFeatures
                    dy = dy_xnk(t:N_T,:,:); % Extract relevant part of dy_xnk for this time instance
                    if size(dy,3) ~= 1 % If measurements are not scalar ... 
                        dy = reshape(permute(dy,[2 1 3]),[ny*(N_T-t)+ny nLin]);
                    end
                end

                % Compute ancestor probabilities of reference trajectory
                xnkt = xnk(:,t);
                for i = 1:N_P  
                    % Compute the probability of transitioning from
                    % x^i_{t-1} (so xn(:,i) since this is before
                    % prediction) to x'_t (so xnk(:,t))
                    if isempty(dynResNorm)
                        eDyn = (xnkt - xn(:,i) - odometry(t-1,:)')'/...
                            chol(dt(t-1)*Q(:,:,t-1),'lower');
                    else
                        eDyn = dynResNorm(xnkt,xn(:,i),odometry(t-1,:),dt(t-1),Q(:,:,t-1));
                    end
                    % Omit constant factor since it's just a constant
                    logwDyn = -1/2*(eDyn*eDyn'); 
                
                    % Compute the probability that future measurements at
                    % locations x'_{t:T} of the reference trajectory match 
                    % the map of the current particle 
                    
                    if ~sparseFeatures
                        % Compute residual and its covariance of joint
                        % prediction into the future
                        SS = dy*P(:,:,i)*dy'+kron(eye(N_T-t+1),R); 
                        yt = reshape(y(t:end,:)',ny*(N_T-t+1),1);
                        e = yt - dy*xl(:,i);
                    else
                        % Compute residual and its covariance of joint
                        % prediction into the future
                        e = [];
                        dy = [];
                        RS = [];
                        for ti=t:N_T
                            % Observation
                            yt = y(ti,:);
                            ind = ~isnan(yt);
                    
                            % Linearize measurement
                            [yhat,dyi] = measModel(xnk(:,ti), xl(:,i));
                  
                            % Stack innovation mean and Jacobians
                            ei = yt' - yhat;
                            e = [e; ei(ind)]; %#ok
                            dy = [dy; dyi(ind,:)]; %#ok
                            RS = blkdiag(RS,R(ind,ind));
                        end
                        
                        % Innovation covariance approx
                        SS = dy*P(:,:,i)*dy' + RS; 

                    end

                    % Solve using Cholesky
                    [cS,flag] = chol(SS,'lower');
                    if flag>0
                        cS = chol(SS+jitter*eye(size(SS,1)),'lower');
                    end
                    v = cS\e;

                    % Compute the measurement weight
                    %logwMeas = -.5*log(det(SS)) - .5*(e'*(SS\e)) =
                    logwMeas = -sum(log(diag(cS))) - .5*(v'*v) - size(e,1)/2*log(2*pi);
                    
                    % Combine weights
                    paNtLog(i) = log(w(i)) + logwDyn + logwMeas;
                end

                % Normalize weights by log-sum-exp trick
                c = max(paNtLog);
                lse = c + log(sum(exp(paNtLog - c)));
                paNt = exp(paNtLog - lse);

                AI(:,t) = paNt; % Save to double check weights, can be removed later
                ai(N_P) = sample(paNt); % ... and sample. 
                xn_pred(:,N_P) = xnkt;
                xl_pred(:,N_P) = xl(:,ai(N_P)); % Use the sampled map
                P_pred(:,:,N_P) = P(:,:,ai(N_P)); % and its covariance.
            end
        
            % Now overwrite current states xn and xl with predicted
            % versions (we need the previous xn's in the computation of the
            % dynamic weights and the previous xl's in the computation of
            % the predictions so we can't do it earlier
            xn = xn_pred;
            xl = xl_pred;
            P = P_pred;
        
            % Save and adapt history (part of Step 11 of Alg. 2)
            xn_traj(:,:,t) = xn;
            xn_traj(:,:,1:t-1) = xn_traj(:,ai,1:t-1);
            
        end

        % Compute importance weights (step 12) 
        yt = y(t,:); % Measurements at time t
        if ~sparseFeatures
            dy = measModel(xn); % Compute gradient of cond linear measurement model
        end
        for i=1:N_P 
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
            else
                % Compute innovations and their covariances
                dyi = squeeze(dy(i,:,:));
                e = yt' - dyi * xl(:,i);
                SS = dyi * P(:,:,i) * dyi' + R;
            end

            % Compute the weights (actually use log-weights)
            % w(i) = exp(-.5*log(det(SS)) - .5*(e'*(SS\e)) - .5*numel(e)*log(2*pi));

            % Compute the log weights
            [cS,flag] = chol(SS,'lower');
            if flag>0
                cS = chol(SS+jitter*eye(size(SS,1)),'lower');
            end            
            v = cS\e;
            logw(i) = -sum(log(diag(cS))) - .5*(v'*v) - .5*numel(e)*log(2*pi);
        end
        
        % Normalize weights (actually use log-weights)
        %w = w ./ sum(w); 

        % Normalize by log-sum-exp trick
        c = max(logw);
        lse = c + log(sum(exp(logw - c)));
        w = exp(logw - lse);  

        % Update parameter means and covariances (step 13)
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
                dyi = squeeze(dy(i,:,:));
                yhat = dyi * xl(:,i);
                % Compute innovations and their covariances
                e = yt' - yhat;
                SS = dyi * P(:,:,i) * dyi' + R;
                % Compute Kalman gain
                [cS,flag] = chol(SS,'lower');
                if flag>0
                    cS = chol(SS+jitter*eye(size(SS,1)),'lower');
                end
                K = P(:,:,i)*((dyi'/cS')/cS);
            end
            % Update linear states and covariances
            xl(:,i)  = xl(:,i) + K*e;
            P(:,:,i) = P(:,:,i) - K*SS*K';
        end
    end
    
    % After computing the nonlinear states and parameters for t = 1...N_T,
    % sample a trajectory and its corresponding map (both mean and
    % covariance) (step 15)
    ak = sample(w);
    xnk = permute(xn_traj(:,ak,:),[1 3 2]);
    xlk = xl(:,ak);
    Pk = P(:,:,ak);
    
    % Add the sampled trajectory and its corresponding map to the outputs
    XNK(:,:,k) = xnk;        
    XLK(:,k) = xlk;
    PK(:,:,k) = Pk;

    % Plot results using
    % 1) sampled nonlinear trajectory
    % 2) sampled linear trajectory
    % 3) smoother iteration
    if ~isempty(makePlots)
      makePlots(xnk,xlk,k,XNK,XLK,PK)
    end    
    
    % Report
    fprintf('Particle smoother iteration %i/%i done.\n',k,N_K)
end

