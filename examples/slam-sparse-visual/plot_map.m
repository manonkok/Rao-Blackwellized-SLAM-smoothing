function plot_map(x,f,fp,fw,cmap,for_paper)
%% PLOT_MAP - Visualize the map and robot's pose in a 2D space
%
% Syntax:
%   plot_map(x,f,fp,fw,cmap,for_paper)
%
% In:
%   x         - State vector containing robot pose and map points
%   f         - Focal length of the camera
%   fp        - Principal point of the camera
%   fw        - Field of view width of the camera
%   cmap      - Color map for the points in the map
%   for_paper - Flag to adjust plot for paper publication (optional)
%
% Description:
%   This function visualizes the map and the robot's pose in a 2D space. It
%   plots the field of view of the camera, the robot's current position and
%   orientation, and the map points. The function can also highlight points
%   that are not visible in the camera's field of view. An optional
%   parameter 'for_paper' can be set to true for a plot style suitable for
%   publications.
%
% Copyright:
%   2023-   Manon Kok and Arno Solin
 
  if nargin<6, for_paper = false; end

  % Euler angle to rot
  eul2rot = @(th) [cos(th) -sin(th); sin(th) cos(th)];

  % Extract pose from state
  p = x(1:2);
  th = x(3);
  R = eul2rot(th);

  % Extract map from state
  map = reshape(x(4:end),2,[]);

  % Field-of-view
  fov = 2*atan2(fw,f);
  
  % Draw frustrum
  scale = .5;
  frust = @(x,th) scale*eul2rot(th)*[0 -tan(fov/2) tan(fov/2) 0; 0 1 1 0] + x;

  % Projection (observation)    
  K = [f fp; 0 1];
  u = K*[R' -R'*p]*[map; ones(1,size(map,2))];
  uu = u(1,:)./u(2,:);
  
  % Plot scene
  %figure(1); clf; hold on
  hold on
  
    % Plot FOV
    xf = eul2rot(th)*[0 -20*tan(fov/2) 20*tan(fov/2) 0; 0 20 20 0] + p;
    h=fill(xf(1,:),xf(2,:),1);
    set(h,'FaceColor',[.5 .5 1],'LineStyle','--','FaceAlpha',.05)
  
    % Show map
    if for_paper
      for i=1:size(map,2)
        scatter(map(1,i),map(2,i),32,cmap(i,:),'filled') 
      end
    else
      scatter(map(1,:),map(2,:),32,cmap,'filled')
    end
    
    % Mark those not visible in camera
    ind = u(2,:)<0 | abs(uu)>fw;
    plot(map(1,ind),map(2,ind),'xk','MarkerSize',10)
    
    % Plot pose and frustrum
    x = frust(p,th);
    plot(p(1),p(2),'ok','MarkerSize',10,'MarkerFaceColor','k')
    plot(x(1,:),x(2,:),'-k','LineWidth',2)
    axis equal
    xlim([-5 5])
    ylim([-5 5])
    box on
    set(gca,'layer','top')

