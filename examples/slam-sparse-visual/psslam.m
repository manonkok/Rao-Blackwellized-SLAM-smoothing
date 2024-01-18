function [traj_mean,m,P,rmse_path,rmse_th,rmse_map] = ...
    psslam(Y,u,map,groundTruth,initPos,initTheta,initMapVar,noiseVar,guessMapVar,f,fp,fw,visualization,savevideo,N_K,N_P)
%% PSSLAM - Perform sparse particle smoother based simultaneous localization and mapping (SLAM)
%
% Syntax:
%   [traj_mean,m,P,rmse_path,rmse_th,rmse_map] = 
%       psslam(Y,u,map,groundTruth,initPos,initTheta,initMapVar,noiseVar,guessMapVar,f,fp,fw,visualization,savevideo,N_K,N_P)
%
% In:
%   Y            - Observations
%   u            - Control inputs
%   map          - Initial map
%   groundTruth  - Ground truth data for path and map
%   initPos      - Initial position of the robot
%   initTheta    - Initial orientation of the robot
%   initMapVar   - Variance for initial map estimation
%   noiseVar     - Noise variance for measurements
%   guessMapVar  - Variance for map point guesses
%   f            - Focal length of the camera
%   fp           - Principal point of the camera
%   fw           - Field of view width of the camera
%   visualization- Flag for enabling visualization
%   savevideo    - Flag for saving the visualization as a video
%   N_K          - Number of smoother samples to generate (default: 10)
%   N_P          - Number of particles to use on each run (default: 100)
%
% Out:
%   traj_mean    - Mean trajectory estimate
%   m            - Map estimate (currently undefined)
%   P            - Covariance matrix of the estimate (currently undefined)
%   rmse_path    - Root mean square error of the path
%   rmse_th      - Root mean square error of the orientation
%   rmse_map     - Root mean square error of the map
%
% Description:
%   This function performs sparse (visual) SLAM using the particle
%   smoothing approach in the paper. It takes sensor observations, control
%   inputs, and initial conditions to estimate the trajectory of the robot
%   and the map of the environment. The function can also visualize the
%   SLAM process and save it as a video. The output includes the mean
%   trajectory estimate, map estimate, and the associated root mean square
%   errors for the path, orientation, and map.
%
% See also:
%   particleSmoother, measurement, calc_rmses
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

%% Prepare visualization

  % Save current random seed
  scurr = rng;

  vidObj = [];

  % Colormap for map points
  rng(0,'twister')
  cmap = parula(size(map,2));
  cmap = cmap(randperm(size(cmap,1)),:);

  if visualization
    figure(1); clf
    set(gcf,'color','w')
    
    % Write video
    if savevideo
        filename = 'loop-ps';
        vidObj = VideoWriter([filename '.mp4'],'MPEG-4');
        vidObj.FrameRate = 10;
        open(vidObj);
    end
  end

  
%% Call external particle filter

  rng(scurr)

  % Number of smoother iterations
  if nargin<15
    N_K = 10;
  end

  % Number of particles
  if nargin<16
    N_P = 100;
  end

  % Dynamic and measurement model
  dynModel = @(xn,dx,dt,Q) xn + dx' + sqrt(dt*Q)*randn(size(xn,1),1);
  measModel = @(xn,xl) measurement([xn(1:3); xl],f,fp,fw,true);
  
  % Odometry (control) and measurements
  odometry = u;
  y = Y';
  dt = 1;
  
  % Initial states
  x0_nonLin = [initPos; initTheta];
  x0_lin = map(:) + sqrt(guessMapVar)*randn(numel(map),N_P);
  P0_lin = initMapVar*eye(size(x0_lin,1));
  Q = blkdiag(0.1^2*eye(2),0.001^2);
  R = noiseVar*eye(size(y,2));
  
  % Use the sparse features flag
  sparseFeatures = true;
  
  % Plotting function 
  if visualization
    makePlots = @(xnk,xlk,k,XNK,XLK,PK) ...
      plot_visual_slam_progress(XNK,k,[],[],XLK,map,groundTruth,f,fp,fw,cmap,vidObj);
  else
    makePlots = [];  
  end
 
  % Run smoother
  [XNK,XLK] = particleSmoother(dynModel,measModel,[],odometry,y,...
    x0_nonLin,x0_lin,P0_lin,Q,R,N_P,N_K,dt,sparseFeatures,makePlots);
  
  % Mean path and map over all samples
  xnk = mean(XNK(:,:,2:end),3);
  xlk = mean(XLK(:,2:end),2);

  % RMSE based on map alignment
  [rmse_path,rmse_th,rmse_map,Z_traj,Z_map] = calc_rmses(map',reshape(xlk,2,[])', ...
    map',reshape(xlk,2,[])',groundTruth',xnk');  

  % Outputs undefined
  m = []; P = []; traj_mean = [];
      
  % Close file if video saved
  if savevideo
    close(vidObj);
  end
    
