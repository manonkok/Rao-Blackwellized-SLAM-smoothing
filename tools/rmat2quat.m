function q = rmat2quat(R)
% RMAT2QUAT - Transforms a rotation matrix R to a quaternion q
%
% Syntax:
%   q = rmat2quat(R)
%
% In:
%   R  - A single 3x3 rotation matrix or an array of 3x3xN rotation matrices
%
% Out:
%   q  - Quaternion or quaternions corresponding to the input rotation matrix 
%        or matrices. Represented as 4x1 vector for a single rotation matrix, 
%        or a 4xN matrix for multiple rotation matrices.
%
% Description:
%   This function transforms a 3x3 rotation matrix (or multiple matrices) 
%   into the corresponding quaternion representation. Due to numerical 
%   challenges, the transformation is performed via the rotation vector. 
%   The function assumes the input is a single rotation matrix or an array 
%   of rotation matrices, and outputs a quaternion for each input matrix.
%
% Examples:
%   q = rmat2quat([1 0 0; 0 1 0; 0 0 1]);
%   qs = rmat2quat(cat(3, [1 0 0; 0 1 0; 0 0 1], [0 -1 0; 1 0 0; 0 0 1]));
%
% See also:
%   expq, logR

if size(R,3)>1
  q = zeros(4,size(R,3));
  for i=1:size(R,3)
    q(:,i) = rmat2quat(R(:,:,i));
  end
else
  phi = logR(R);
  q = expq(phi/2);
end

end
