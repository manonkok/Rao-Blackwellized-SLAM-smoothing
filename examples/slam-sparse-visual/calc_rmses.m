function [rmse_path,rmse_th,rmse_map,Z_traj,Z_map] = calc_rmses(truth,estimate,map,map_est,traj,traj_est)
%% CALC_RMSES - RMSE error metrics for the 2d visual SLAM case
%
% Syntax:
%   [rmse_path,rmse_th,rmse_map,Z_traj,Z_map] = calc_rmses(truth,estimate,map,map_est,traj,traj_est)
%
% In:
%   truth    - True points to align to (n x 2)
%   estimate - Estimated points to align to (n x 2)
%   map      - True map feature point locations (m x 2)
%   map_est  - Estimated map feature points (m x 2)
%   traj     - True trajectory (translation + angle) (k x 3)
%   traj_est - Estimated trajectory (translation + angle) (k x 3)
%   
% Out:
%
%   rmse_path - RMSE in 2d path (RMSE of distance offsets)
%   rmse_th   - RMSE in angle
%   rmse_map  - RMSE in 2d map point locations (RMSE of distance offset)
%
% Description:
%   Calculate RMSE metrics for the visual SLAM in 2D. The error metrics are
%   based on the assumption that absolute scale is unobservable and thus
%   the estimates are first aligned based on mathcing 'truth' and 'estimate'
%   with PROCRUSTES analysis (match scale, rotation, translation).
%
%   One might typically want to do the matching by supplying the final map
%   points and grund-truth map points as the 'truth' and 'estimate',
%   alternatively the path estimate can be used.
%

%%

  % Match scale, rotation, and translation with procrustes
  [D, Z, TRANSFORM] = procrustes(truth, estimate);
  
  % Transform trajectory
  Z = TRANSFORM.b * traj_est(:,1:2) * TRANSFORM.T + TRANSFORM.c(1,:);
  
  % Transform angle
  Z_th = traj(:,3) * nan; % Not used/implemented
  
  % Use the same transformation to match the map points
  Z_map =  TRANSFORM.b * map_est * TRANSFORM.T + TRANSFORM.c(1,:);
  
  % RMSE (path)
  d = sqrt(sum((traj(:,1:2)-Z).^2,2));
  rmse_path = sqrt(mean(d.^2));
  
  % RMSE (angle)
  rmse_th = sqrt(mean((unwrap(traj(:,3))-Z_th).^2));
      
  % RMSE (map)
  d = sqrt(sum((map-Z_map).^2,2));
  rmse_map = sqrt(mean(d.^2));
  
  % Return also aligned
  Z_traj = [Z Z_th];
  
  
  