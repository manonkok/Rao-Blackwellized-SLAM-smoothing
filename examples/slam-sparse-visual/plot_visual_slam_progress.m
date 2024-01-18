function plot_visual_slam_progress(xn,xl_max,traj_mean,traj_max,xl,map,groundTruth,f,fp,fw,cmap,vidObj)
%% PLOT_VISUAL_SLAM_PROGRESS - Visualize the progress of a visual SLAM algorithm
%
% Syntax:
%   plot_visual_slam_progress(xn,xl_max,traj_mean,traj_max,xl,map,groundTruth,f,fp,fw,cmap,vidObj)
%
% In:
%   xn          - Nonlinear state trajectories over time
%   xl_max      - Maximum-weight linear state (~map) or current iteration index for smoother
%   traj_mean   - Mean trajectory estimate
%   traj_max    - Highest-weight trajectory estimate
%   xl          - Linear state trajectories over time
%   map         - Initial map
%   groundTruth - Ground truth data for path and map
%   f           - Focal length of the camera
%   fp          - Principal point of the camera
%   fw          - Field of view width of the camera
%   cmap        - Color map for the points in the map
%   vidObj      - VideoWriter object for saving the visualization as a video (optional)
%
% Description:
%   This function visualizes the progress of a visual SLAM (Simultaneous
%   Localization and Mapping) algorithm. It shows the estimated
%   trajectories, map points, and the robot's current view. The function
%   can plot results from both particle filter and smoother
%   implementations. It can also save the visualizations as a video if a
%   VideoWriter object is provided.
%
% See also:
%   plot_map, calc_rmses
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

  % Load background
  [img,~,alpha] = imread('furniture.png');

  % Gridding
  xi = linspace(-5,5,size(img,1));
  yi = xi;  
  
  % Flip
  %img = img(end:-1:1,:,:);
  alpha = alpha(end:-1:1,:);
  
  figure(1); clf; hold on
  
  % Align
  %{
  [rmse_path,rmse_th,rmse_map,Z_traj,Z_map] = calc_rmses(map',reshape(xl_max,2,[])', ...
      map',reshape(xl_max,2,[])',groundTruth',traj_mean');
  [rmse_path,rmse_th,rmse_map,Z_traj_max,Z_map] = calc_rmses(map',reshape(xl_max,2,[])', ...
      map',reshape(xl_max,2,[])',groundTruth',traj_max');
  %}
  
  % Background
  image(xi,yi,(255-.3*alpha))
  colormap(gray), caxis([0 255])
  
  % Are we plotting filter or smoother results?
  if numel(xl_max)==1 % Plotting smoother

    % Smoother iteration
    k = xl_max;

    % Plot ground-truth
    plot(groundTruth(1,:),groundTruth(2,:),'--r')
    
    % Plot trajectory (all samples)
    for i=1:k

      % Use Procrustes for alignment (rotation and scale lost)
      [~,~,~,Z_traj,Z_map] = calc_rmses(map',reshape(xl(:,i),2,[])', ...
        map',reshape(xl(:,i),2,[])',groundTruth',xn(:,:,i)');

      % Plot aligned
      plot(Z_traj(:,1),Z_traj(:,2),'-k')

      % Plot map points as scatter
      for j=1:size(Z_map,1)
        plot(Z_map(j,1),Z_map(j,2),'.','color',cmap(j,:))
      end

      % Plot unaligned
      %plot(xn(1,:,i),xn(2,:,i),'-k')
      
    end

    % Plot final view
    plot_map([Z_traj(end,1:2)'; xn(end,end,k); reshape(Z_map',[],1)],f,fp,fw,cmap);
   
  else % Plotting filter

    % How far are we?
    t = find(~isnan(traj_mean(1,:)),1,'last');

    % Plot unaligned
    plot_map([traj_mean(:,t); xl_max],f,fp,fw,cmap);
  
    % Plot ground-truth
    plot(groundTruth(1,:),groundTruth(2,:),'--r')
  
    % Plot trajectory (mean)
    plot(traj_mean(1,1:t),traj_mean(2,1:t),'-k')
  
    % Plot trajectory (max)
    %plot(Z_traj_max(1:t,1),Z_traj_max(1:t,2),'-r')
  
    % Points
    plot(xn(1,:),xn(2,:),'.k')

    % Plot map points as scatter
    for i=1:2:(size(xl,1)/2)
      plot(xl(i,:),xl(i+1,:),'.','color',cmap(i,:))
    end
    
  end
    
  axis xy image off
  xlim([-5 5]), ylim([-5 5])
  
  drawnow
  
  % Write to video
  if isobject(vidObj)
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
