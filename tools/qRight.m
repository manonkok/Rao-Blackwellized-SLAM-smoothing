function qR = qRight(q)
%% QRIGHT - Computes the right multiplication matrices of a vector of quaternions
%
% Syntax:
%   qR = qRight(q)
%
% In:
%   q  - Quaternions, represented as Nx4 matrix. For a single quaternion, 
%        both 1x4 and 4x1 formats are accepted.
%
% Out:
%   qR - The right multiplication matrices for each quaternion. For a single 
%        quaternion, qR is a 4x4 matrix. For multiple quaternions, qR is an 
%        array of 4x4 matrices, with size 4x4xN.
%
% Description:
%   This function computes the right multiplication matrices for a single 
%   quaternion or a vector of quaternions. Each quaternion is represented 
%   as a row vector in Nx4 format, and the function outputs a 4x4 matrix 
%   for each quaternion which represents its right multiplication matrix.
%
% Examples:
%   qR = qRight([0.7071, 0, 0.7071, 0]);
%   qRs = qRight([1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0]);
%
% See also:
%   mcross, multiprod

if any(size(q) == 1)
    if size(q,1) == 1
        q = q';
    end
    qR = [q(1) , -q(2:4)' ; q(2:4) , q(1)*eye(3)-mcross(q(2:4))];
else
    N = size(q,1); % Determine N
    qR = [reshape(q(:,1),1,1,N), reshape(-q(:,2:4)',[1,3,N]) ; ...
        reshape(q(:,2:4)',[3,1,N]), ...
        multiprod(reshape(q(:,1),[1 1 N]),repmat(eye(3),[1 1 N]))-mcross(q(:,2:4))];
end

