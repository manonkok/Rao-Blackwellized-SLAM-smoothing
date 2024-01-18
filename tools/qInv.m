function q = qInv(q)
%% QINV - Computes the quaternion inverse of q
%
% Syntax:
%   q = qInv(q)
%
% In:
%   q  - Quaternions, represented as Nx4 matrix. For a single quaternion, 
%        both 1x4 and 4x1 formats are accepted.
%
% Out:
%   q  - The inverse of each quaternion provided. The output format is the 
%        same as the input format.
%
% Description:
%   This function computes the inverse of one or more quaternions. The 
%   quaternion inverse is used in operations like quaternion division or 
%   rotating a vector by the inverse rotation. The function accepts a single 
%   quaternion or an array of quaternions, and it outputs the inverse in the 
%   same format as the input.
%
% Examples:
%   qInv = qInv([0.7071, 0, 0.7071, 0]);
%   qInvs = qInv([1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0]);
%

if any(size(q) == 1)
    q(2:4) = -q(2:4);
else
    q = [q(:,1) -q(:,2) -q(:,3) -q(:,4)];
end

end