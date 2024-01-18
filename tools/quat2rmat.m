function R = quat2rmat(q)
%% QUAT2RMAT - Converts quaternions into rotation matrices
%
% Syntax:
%   R = quat2rmat(q)
%
% In:
%   q  - Quaternions, represented as Nx4 matrix. For a single quaternion, 
%        both 1x4 and 4x1 formats are accepted.
%
% Out:
%   R  - Rotation matrices corresponding to the quaternions. If multiple 
%        quaternions are provided, R is a 3x3xN array of rotation matrices.
%
% Description:
%   This function converts one or more quaternions into their corresponding 
%   rotation matrices. A quaternion is represented as [q0, q1, q2, q3], and 
%   the function supports a batch of quaternions in an Nx4 matrix format. 
%   The output is a 3x3 rotation matrix for each quaternion, which describes 
%   the rotation in 3D space.
%
% Examples:
%   R = quat2rmat([0.7071, 0, 0.7071, 0]);
%   Rs = quat2rmat([1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0]);
%

if any(size(q) == 1)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    R = [q0^2+q1^2-q2^2-q3^2, 2*q1*q2 - 2*q0*q3 ,2*q1*q3 + 2*q0*q2 ; ...
        2*q1*q2 + 2*q0*q3, q0^2 - q1^2 + q2^2 - q3^2, 2*q2*q3 - 2*q0*q1 ; ...
       2*q1*q3 - 2*q0*q2, 2*q2*q3 + 2*q0*q1, q0^2 - q1^2 - q2^2 + q3^2];
else
    nQuat = size(q,1);
    q0 = q(:,1); q1 = q(:,2); q2 = q(:,3); q3 = q(:,4);
    R = [q0.^2+q1.^2-q2.^2-q3.^2, 2*q1.*q2 - 2*q0.*q3 ,2*q1.*q3 + 2*q0.*q2 , ...
        2*q1.*q2 + 2*q0.*q3, q0.^2 - q1.^2 + q2.^2 - q3.^2, 2*q2.*q3 - 2*q0.*q1 , ...
       2*q1.*q3 - 2*q0.*q2, 2*q2.*q3 + 2*q0.*q1, q0.^2 - q1.^2 - q2.^2 + q3.^2];
    R = permute(reshape(R',[3 3 nQuat]),[2 1 3]);
end

end