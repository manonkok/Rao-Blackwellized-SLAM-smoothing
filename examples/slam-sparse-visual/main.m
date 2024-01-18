%% VISUAL SLAM EXAMPLE IN 2D
%
% This example demonstrates the use of the filtering vs. smoothing
% approaches presented in the paper (see Sec. 5.3). The model corresponds
% to a pinhole camera model observing landmark points while moving in a 
% 2D space.
% 
% The codes allow to visualize the progression of the particle filter in
% real-time, while the smoother output is given as samples of the final
% smoothing distribution.
%
% Copyright:
%   2023-   Manon Kok and Arno Solin
  
%% Dependencies and options

  clear, close all, clc

  % Add paths
  addpath ../../tools/
  addpath ../../src/

  % Choose seed for reproducability
  seed = 42;

  % Set noise scales
  initMapVar = 4^2;
  noiseVar = .1^2;
  guessMapVar = 1^2;

  % Load data and ground-truth to compare to
  rng(seed,'twister')
  [Y,u,map,groundTruth,initPos,initTheta,f,fp,fw] = load_data();

  % Show animation
  showAnimation = true;
  saveVideo = false;


%% Run filter

  rng(seed,'twister')
  
  % Run the PF SLAM
  [x_traj,m,P,rmse_path,rmse_th,rmse_map] = ...
        pfslam(Y,u,map,groundTruth,initPos,initTheta,initMapVar, ...
               noiseVar,guessMapVar,f,fp,fw,showAnimation,saveVideo);

  % Report RMSE (aligned with Procrustes first)
  fprintf('PF SLAM:\n  RMSE (path): %.2f\n  RMSE (map): %.2f\n', ...
        rmse_path,rmse_map)


%% Run smoother

  rng(seed,'twister')    
  
  % Run the PS SLAM (with less samples)
  [x_traj,m,P,rmse_path,rmse_th,rmse_map] = ...
        psslam(Y,u,map,groundTruth,initPos,initTheta,initMapVar, ...
               noiseVar,guessMapVar,f,fp,fw,showAnimation,saveVideo,10,10);

  % Report RMSE (aligned with Procrustes first)
  fprintf('PS SLAM:\n  RMSE (path): %.2f\n  RMSE (map): %.2f\n', ...
        rmse_path,rmse_map)

