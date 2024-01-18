function [y,dy,not_visible] = measurement(x,f,fp,fw,onlyLin)
%% MEASUREMENT - Compute the projection of points in a map given a robot pose
%
% Syntax:
%   [y,dy,not_visible] = measurement(x,f,fp,fw,onlyLin)
%
% In:
%   x         - State vector containing robot pose and map points
%   f         - Focal length of the camera
%   fp        - Principal point of the camera
%   fw        - Field of view width of the camera
%   onlyLin   - Flag to return only linear derivatives (optional)
%
% Out:
%   y         - Projected points on the image plane
%   dy        - Derivative of the projection w.r.t. the state
%   not_visible - Indicator for points not visible in the camera view
%
% Description:
%   This function computes the projection of 2D points (map) onto the image 
%   plane of a camera, given the pose of the robot. The state vector 'x' 
%   includes the robot pose (position and orientation) and the map points. 
%   The function returns the projected points 'y', their derivatives 'dy', 
%   and a boolean array 'not_visible' indicating which points are not visible 
%   in the camera's field of view. Optionally, if 'onlyLin' is true, 
%   'dy' will only contain derivatives with respect to the linear part of 
%   the state.
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

  % Euler angle to rot
  eul2rot = @(th) [cos(th) -sin(th); sin(th) cos(th)];

  % Extract pose from state
  p = x(1:2);
  th = x(3);
  R = eul2rot(th);

  % Extract map from state
  map = reshape(x(4:end),2,[]);

  % Form camera matrix  
  K = [f fp; 0 1];

  % Do projection
  u = K*[R' -R'*p]*[map; ones(1,size(map,2))];

  % Output measurements  
  y = (u(1,:)./u(2,:))';
  
  % Output derivatives
  dy = nan(size(y,1),size(x,1));
  
  % These cannot be oberved (outside fov or behind cam)
  not_visible = u(2,:)'<0 | abs(y)>fw;
  
  % Derivatives
  div = (map(2,:)*cos(th) - p(2)*cos(th) - map(1,:)*sin(th) + p(1)*sin(th)).^2;
  
  % diff(y,p1)
  dy(:,1) = -(f*(map(2,:) - p(2)))'./div(:);
  
  % diff(y,p2)
  dy(:,2) =  (f*(map(1,:) - p(1)))'./div(:);
  
  % diff(y,th)
  dy(:,3) = (f*(map(1,:).^2 - 2*map(1,:)*p(1) + map(2,:).^2 - 2*map(2,:)*p(2) + p(1)^2 + p(2)^2))'./div(:);
  
  % Then the map...
  % diff(y,m1)
  dym1 = (f*(map(2,:) - p(2)))./div;

  % diff(y,m1)
  dym2 = -(f*(map(1,:) - p(1)))./div;

  % Assign
  dy(:,4:2:end) = diag(dym1(:));
  dy(:,5:2:end) = diag(dym2(:));

  % Discard derivatives that are possibly not needed
  if nargin>4 && onlyLin
    dy = dy(:,4:end);
  end
  
  % Derivations:
  
  % diff(y,p1)
  % -(f*(m2 - p2))/(m2*cos(th) - p2*cos(th) - m1*sin(th) + p1*sin(th))^2

  % diff(y,p2)
  % (f*(m1 - p1))/(m2*cos(th) - p2*cos(th) - m1*sin(th) + p1*sin(th))^2

  % diff(y,th)
  % (f*(m1^2 - 2*m1*p1 + m2^2 - 2*m2*p2 + p1^2 + p2^2))/(m2*cos(th) - p2*cos(th) - m1*sin(th) + p1*sin(th))^2

  % diff(y,m1)
  % (f*(m2 - p2))/(m2*cos(th) - p2*cos(th) - m1*sin(th) + p1*sin(th))^2
 
  % diff(y,m2)
  % -(f*(m1 - p1))/(m2*cos(th) - p2*cos(th) - m1*sin(th) + p1*sin(th))^2
 