function pL = qLeft(q)
%% QLEFT - Computes the left multiplication matrices of a vector of quaternions
%
% Syntax:
%   pL = qLeft(q)
%
% In:
%   q  - Quaternions, represented as Nx4 matrix. For a single quaternion, 
%        both 1x4 and 4x1 formats are accepted.
%
% Out:
%   pL - The left multiplication matrices for each quaternion. For a single 
%        quaternion, pL is a 4x4 matrix. For multiple quaternions, pL is an 
%        array of 4x4 matrices, with size 4x4xN.
%
% Description:
%   This function computes the left multiplication matrices for a single 
%   quaternion or a vector of quaternions. Each quaternion is represented 
%   as a row vector in Nx4 format, and the function outputs a 4x4 matrix 
%   for each quaternion which represents its left multiplication matrix.
%
% Examples:
%   pL = qLeft([0.7071, 0, 0.7071, 0]);
%   pLs = qLeft([1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0]);
%
% See also:
%   mcross, multiprod

% Determine left matrix multiplication matrix
if any(size(q) == 1)
    if size(q) == [1 4]
        q = q';
    end
    pL = [q(1) , -q(2:4)' ; q(2:4) , q(1)*eye(3)+mcross(q(2:4))];
else
    N = size(q,1); % Determine N
    pL = [reshape(q(:,1),1,1,N), reshape(-q(:,2:4)',[1,3,N]) ; ...
            reshape(q(:,2:4)',[3,1,N]), ...
            multiprod(reshape(q(:,1),[1 1 N]),repmat(eye(3),[1 1 N]))+mcross(q(:,2:4))];
end