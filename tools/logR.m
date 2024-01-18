function phi = logR(R)
%% LOGR - Computes the vector representation of the matrix logarithm of a rotation matrix
%
% Syntax:
%   phi = logR(R)
%
% In:
%   R  - A 3x3 rotation matrix
%
% Out:
%   phi - The 3-dimensional vector representation of the matrix logarithm of R
%
% Description:
%   This function computes the matrix logarithm of a 3x3 rotation matrix R 
%   and then extracts the components of the corresponding vector phi. The 
%   matrix logarithm of a rotation matrix represents the axis-angle form 
%   of the rotation. The output phi is the axis of rotation scaled by the 
%   rotation angle, often used in robotics and 3D transformations.
%
% Examples:
%   phi = logR([1 0 0; 0 1 0; 0 0 1]);
%   phi = logR([0 -1 0; 1 0 0; 0 0 1]);
%
% See also:
%   expm, logm
%

phicross = logm(R);
phi = [phicross(3,2) ; phicross(1,3) ; phicross(2,1)];

end