function [traj_max,traj_mean] = particleFilterLocalization(dynModel,measModel,odometry,y,...
    x0_nonLin,Q,R,N_P,dt,makePlots)
% PARTICLEFILTERLOCALIZATION - Particle filter for localization
%
% Syntax:
%   [traj_max,traj_mean,xl_max,traj_sample,xn_traj] = 
%     particleFilter(dynModel,measModel,odometry,y,...
%                 x0_nonLin,x0_lin,P0_lin,Q,R,N_P,sparseFeatures,makePlots)
%
% In:
%   dynModel    - Dynamical model function handle, xn = @(xn,dx,Q)
%   measModel   - Measurement model function handle, w = @(y,xn)
%   odometry    - Dynamical model control observations [T x c]
%   y           - Observations [T x m]
%   x0_nonLin   - Initial non-linear state [dn x 1]
%   x0_lin      - Initial linear states [dl x N_P]
%   P0_lin      - Initial linear state cov [dl x dl]
%   Q           - Process noise cov (can be different for each step)
%   R           - Measurement noise cov (not used)
%   N_P         - Number of particles
%   dt          - Time between two time steps (can be different for each)
%   makePlots   - Visualization function (default: none)
%
% Out:
%   traj_max    - Highest-weight trajectory (nonlinear)
%   traj_mean   - Weighted-mean trajectory (nonlinear) 
%
% Description:
%   Run particle filter for the magnetic terrain matching problem. This
%   code combines apporoaches from [1]-[2].
%
% References:
%   [1] Arno Solin, Simo Särkkä, Juho Kannala, and Esa Rahtu (2016). 
%       Terrain navigation in the magnetic landscape: Particle filtering 
%       for indoor positioning. In Proceedings of the European Navigation 
%       Conference (ENC). Pages 1–9. IEEE.
%   [2] Manon Kok and Arno Solin (2018). Scalable magnetic field SLAM in 
%       3D using Gaussian process maps. In Proceedings of the International 
%       Conference on Information Fusion (FUSION). Pages 1353–1360. 
%       Cambridge, UK.
%
% See also:
%   particleFilter
%
% Copyright:
%   2022-   Manon Kok and Arno Solin
%

%%

  % Initial weights
  w = 1/N_P * ones(1,N_P);
  
  % Initial states
  if size(x0_nonLin,2) > 1  
    xn = x0_nonLin;
  else
    xn = repmat(x0_nonLin,1,N_P); % Nonlinear states
  end
  
  % Extract some parameters
  nNonLin = size(x0_nonLin,1);
  N_T = size(y,1);
  
  % Allow for both time-varying and constant Q
  if size(Q,3) == 1 
      Q = repmat(Q,[1 1 N_T]);
  end
  
  % Allow for both time-varying and constant time step
  if length(dt) == 1 
      dt = dt * ones(N_T,1);
  end
    
  % Preallocate trajectories
  traj_max = nan(nNonLin,N_T);
  traj_mean = nan(nNonLin,N_T);
  yhattraj = nan(size(y,2), N_T);
  xn_traj = zeros(nNonLin, N_P, N_T);
  xn_traj(:,:,1) = xn;
  ai = zeros(N_P,1);
    
  % For each measurement / time step
  for t=1:N_T
      
    % Copy old ones
    xn_ = xn;
      
    % Particle filter prediction
    if t ~= 1
        for i = 1:N_P
            % Draw ancestor index...
            ai(i) = sample(w); 
            % ... and propagate through dynamics
            xn(:,i) = dynModel(xn_(:,ai(i)),odometry(t-1,:),dt(t-1),Q(:,:,t-1)); % Draw sample
        end
    end
    
    % Save trajectory with shuffled ancestor indices, e.g. to visualise
    % path degeneracy
    if t ~= 1
        xn_traj(:,:,t) = xn;
        xn_traj(:,:,1:t-1) = xn_traj(:,ai,1:t-1);
    end
    
    % Compute the importance weights
    yt = y(t,:);
    
    % Compute the weights
    w = measModel(yt,xn);
    
    % Check if the filter diverged
    if sum(w)<=1e-12
        disp(['Weights filter close to zero at t=',num2str(t),' !!!']);
    end
    
    % Normalize weights
    w = w ./ sum(w); 
    
    % Save trajectories
    [~,iw_max] = max(w);
    traj_max(:,t) = xn(:,iw_max);  
    traj_mean(:,t) = sum(xn.*w,2);

    % Save trajectory
    traj_max(:,t) = xn(:,iw_max);
    
    % Plot results
    if ~isempty(makePlots)
      makePlots(xn,traj_max,yhattraj,xn_traj,traj_mean)
    end
  end
  
end