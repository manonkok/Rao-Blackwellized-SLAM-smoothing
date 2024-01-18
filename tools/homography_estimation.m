function [A,C] = homography_estimation(X,Y)
%% HOMOGRAPHY_ESTIMATION - Estimate the coefficients in the transform
%
% Syntax:
%   [A,C] = homography_estimation(X,Y)
%
% In:
%   X - Inputs (n x 2, points in camer image)
%   Y - Outputs (n x 2, )
%
% Out:
%   A - Transformation matrix
%   C - Transformation vector
%
% Description:
%   This function estimates the coefficient matrices A and C in the
%   following kind of homographic projections (also known as Moebius 
%   transforms in some special cases). This is the transform that
%   compensates for the prespective error in the web camera images, 
%   and maps the camera-observed coordinates to real-world coordinates.
%   
%   The homographic projection is
%       y = A*x / (c*x),
%   where y = [y1, y2, 1]' are the true coordinates and x = [x1, x2, 1]  
%   are the coordinates seen by the camera. Some details on this type
%   of projections are available on http://en.wikipedia.org/wiki/Homography
%   
%   Details on how the parameters are estimated can be found in [1] and
%   described on http://xenia.media.mit.edu/~cwren/interpolator/
%
% References:
%   [1] A. Criminisi, I. Reid, A. Zisserman (1997). A Plane Measuring 
%       Device. Available online.
%
% Author:
%   2015 Arno Solin
  
  B = [ Y ones(size(Y,1),1) zeros(size(Y,1),3) -Y(:,1:2).*X(:,[1 1]) ...
        zeros(size(Y,1),3) Y ones(size(Y,1),1) -Y(:,1:2).*X(:,[2 2])];
  B = reshape (B', 8 , [] )';
  D = reshape (X', [] , 1 );
  l = B \ D;
  A = reshape([l(1:6)' 0 0 1 ],3,3)';
  C = [l(7:8)' 1];
  
  
%%

  % Debugging
  %X = ijxy(:,3:4);
  %Y = ijxy(:,1:2);
  
  %X = [X; X];
  %Y = [Y; Y];

  % This is the projection we want to have
  %    y = A*x / (c*x);
  % where y are the true coordinates and x are the coordinates seen by the
  % camera.

% A =
% 
%    -0.0074   -0.0012    3.8683
%    -0.0026    0.0157   -2.3910
%          0         0    1.0000
% 
% C =
% 
%     0.0005    0.0043    1.0000
% 

%%
% 
%   Xp=xy(:,1);
%   Yp=xy(:,2);
%   X=ij(:,1);
%   Y=ij(:,2);
% 
%   B = [ X Y ones(size(X)) zeros(size(X,1),3)        -X.*Xp -Y.*Xp ...
%         zeros(size(X,1),3)        X Y ones(size(X)) -X.*Yp -Y.*Yp ];
%   B = reshape (B', 8 , 8 )';
%   D = [ Xp , Yp ];
%   D = reshape (D', 8 , 1 );
%   l = inv(B' * B) * B' * D;
%   A = reshape([l(1:6)' 0 0 1 ],3,3)';
%   C = [l(7:8)' 1];
% 
%   %t=A*[x;y;1]/(C*[x;y;1]);
%   
%   % Make grid
%   x = [linspace(min(xy(:,1)),max(xy(:,1)),5) nan];
%   y = [linspace(min(xy(:,2)),max(xy(:,2)),5) nan];
%   [X,Y] = meshgrid(x,y);
%   
%   % Transform
%   foo = [X(:) Y(:) ones(size(X(:)))]/A;
%   foo = bsxfun(@rdivide,foo,foo(:,3))';
%   I = reshape(foo(1,:),size(X));
%   J = reshape(foo(2,:),size(X));
