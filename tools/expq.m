function eq = expq(phi)
%% EXPQ - Computes the quaternion exponential of the vector phi
%
% Syntax:
%   eq = expq(phi)
%
% In:
%   phi  - A single 3-dimensional vector or an Nx3 matrix of vectors
%
% Out:
%   eq   - Quaternion exponential of 'phi'. Represented as a 4x1 vector for 
%          a single input vector, or an Nx4 matrix for multiple vectors.
%
% Description:
%   This function computes the quaternion exponential of a 3-dimensional 
%   vector or a batch of such vectors. The quaternion exponential maps a 
%   vector in 3D space to the quaternion space. This is used in the context 
%   of rotations, where the input vector represents a rotation axis scaled 
%   by the rotation angle. The output quaternion represents the same rotation.
%

if any(size(phi) == 1)
    mag_phi = norm(phi);
    norm_phi = phi./(mag_phi + (mag_phi == 0));
    if size(norm_phi) == [1 3]
        norm_phi = norm_phi';
    end

    eq = [cos(mag_phi) ; norm_phi*sin(mag_phi)];
    if eq(1) < 0
        eq = - eq;
    end
else
    mag_phi = sqrt(phi(:,1).^2 + phi(:,2).^2 + phi(:,3).^2);
    norm_phi = phi./repmat((mag_phi + (mag_phi == 0)),[1 3]);
    eq = [cos(mag_phi) , norm_phi.*repmat(sin(mag_phi),[1 3])];
    eq(eq(:,1)<=0,:) = -eq(eq(:,1)<=0,:);
end