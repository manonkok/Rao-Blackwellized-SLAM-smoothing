function e = quat2euler(q)
%% QUAT2EULER - Converts quaternions into Euler angles
%
% Syntax:
%   e = quat2euler(q)
%
% In:
%   q  - Quaternions, represented as Nx4 matrix. For a single quaternion, 
%        both 1x4 and 4x1 formats are accepted.
%
% Out:
%   e  - Euler angles corresponding to the quaternions, in degrees. The 
%        angles are in the order of yaw, pitch, and roll.
%
% Description:
%   This function converts one or more quaternions into their corresponding 
%   Euler angles. Euler angles are a method for representing the orientation 
%   of a rigid body with a sequence of three rotations around different 
%   axes. The output is in degrees and represents yaw, pitch, and roll angles.
%
% Examples:
%   e = quat2euler([0.7071, 0, 0.7071, 0]);
%   es = quat2euler([1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0]);
%

if any(size(q) == 1)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
else
    q0 = q(:,1); q1 = q(:,2); q2 = q(:,3); q3 = q(:,4);
end

e = 180/pi*[atan2( (2*q2.*q3 - 2*q0.*q1) , (2*q0.^2 + 2*q3.^2 - 1 ) ) , ...
    - real(asin(2*q1.*q3 + 2*q0.*q2)) , ...
    atan2( (2*q1.*q2 - 2*q0.*q3) , (2*q0.^2 + 2*q1.^2 -1 ) ) ];

% e = 180/pi*[atan( (2*q2.*q3 - 2*q0.*q1) ./ (2*q0.^2 + 2*q3.^2 - 1 ) ) , ...
%     - real(asin(2*q1.*q3 + 2*q0.*q2)) , ...
%     atan( (2*q1.*q2 - 2*q0.*q3) ./ (2*q0.^2 + 2*q1.^2 -1 ) ) ];

% [yaw,pitch,roll] = quat2angle( qInv(q) );
% e2 = 180/pi*[roll, pitch, yaw];

end
