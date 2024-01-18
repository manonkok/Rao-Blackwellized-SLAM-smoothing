function lq = logq(q)
%% LOGQ - Computes the quaternion logarithm of the quaternion q
%
% Syntax:
%   lq = logq(q)
%
% In:
%   q  - A single quaternion, represented as either a 1x4 or 4x1 vector.
%
% Out:
%   lq - The quaternion logarithm of 'q'.
%
% Description:
%   This function computes the quaternion logarithm of a given quaternion. 
%   The quaternion logarithm is used in various applications involving 
%   rotations, particularly in interpolating rotations.
%
% Examples:
%   lq = logq([0.7071, 0, 0.7071, 0]);
%   lq = logq([1, 0, 0, 0]);
%
% See also:
%   expq

if any(size(q) == 1)
    if q(1) < 0
        q = -q;
    end
    na = acos(q(1));
    lq = na*q(2:4) / (sin(na) + (na==0));
else
    q(q(:,1)<=0,:) = -q(q(:,1)<=0,:);
    na = acos(q(:,1));
    lq = repmat(na,[1 3]) .* q(:,2:4) ./ repmat((sin(na) + (na==0)),[1 3]);   
end