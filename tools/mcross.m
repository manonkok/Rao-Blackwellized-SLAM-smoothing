function M = mcross(v)
%% MCROSS - Computes the matrix cross product of vector v: M = [v x]
%
% Syntax:
%   M = mcross(v)
%
% In:
%   v  - Vectors, represented as Nx3 matrix. For a single vector, both 1x3 
%        and 3x1 formats are accepted.
%
% Out:
%   M  - The matrix cross product of each vector provided. For a single 
%        vector, M is a 3x3 matrix. For multiple vectors, M is an array of 
%        3x3 matrices, with size 3x3xN.
%
% Description:
%   This function computes the matrix cross product for a single 3D vector 
%   or a batch of 3D vectors. The matrix cross product M of a vector v is a 
%   matrix such that the cross product of v with another vector w can be 
%   computed as M*w. This is useful in various applications in physics and 
%   engineering, particularly in the representation of angular velocities.
%

% If v = 3x1, transpose
if size(v,2) == 1
    v = v';
end

% Determine N
N = size(v,1);

% Determine matrix cross product
if N == 1
    M = [0,-v(3),v(2);...
        v(3),0,-v(1);...
        -v(2),v(1),0];
else
    % Note that the matrix is filled via the unvec operation
    M = reshape([zeros(N,1),v(:,3),-v(:,2),...
        -v(:,3),zeros(N,1),v(:,1),...
        v(:,2),-v(:,1),zeros(N,1)]',3,3,N);
end

end

