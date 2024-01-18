function J=JacobianPhi3D(x,N_m,xl,xu,yl,yu,zl,zu,Indices)
%% JACOBIANPHI3D - Computes the Jacobian matrix for basis functions in 3D space
%
% Syntax:
%   J = JacobianPhi3D(x, N_m, xl, xu, yl, yu, zl, zu, Indices)
%
% In:
%   x       - Input points for the Jacobian computation (3xN matrix)
%   N_m     - Number of basis functions
%   xl, xu  - Lower and upper bounds of the x-dimension
%   yl, yu  - Lower and upper bounds of the y-dimension
%   zl, zu  - Lower and upper bounds of the z-dimension
%   Indices - Indices for the basis functions
%
% Out:
%   J       - Jacobian matrix of the basis functions at the input points 
%             (3x3xN_mxN matrix)
%
% Description:
%   This function computes the Jacobian matrix of basis functions for a 
%   squared exponential kernel in 3D space. It is used in Gaussian Process 
%   regression and other applications where spatial basis functions are 
%   involved. The function takes input points, the number of basis functions, 
%   domain bounds, and basis function indices to compute the Jacobian.
%
% Copyright:
%   2022   Frida Viset

j=Indices;

N=size(x,2);
J=zeros(3,3,N_m,N);

a=[xl yl zl];
b=[xu yu zu];

f=zeros(N_m,3);
for d=1:3
    f(:,d)=(pi*j(:,d))./(b(d)-a(d));
end


for i=1:N
    
    s=zeros(N_m,3);
    c=zeros(N_m,3);
    for d=1:3
        core=pi*j(:,d)*(x(d,i)-a(d))./(b(d)-a(d));
        mult=1./sqrt(0.5*(b(d)-a(d)));
        s(:,d)=sin(core)*mult;
        c(:,d)=cos(core)*mult;
    end
    
    J(1,1,:,i)=-f(:,1).^2.*s(:,1).*s(:,2).*s(:,3);
    J(1,2,:,i)=f(:,1).*f(:,2).*c(:,1).*c(:,2).*s(:,3);
    J(1,3,:,i)=f(:,1).*f(:,3).*c(:,1).*s(:,2).*c(:,3);
    J(2,1,:,i)=f(:,2).*f(:,1).*c(:,1).*c(:,2).*s(:,3);
    J(2,2,:,i)=-f(:,2).^2.*s(:,1).*s(:,2).*s(:,3);
    J(2,3,:,i)=f(:,2).*f(:,3).*s(:,1).*c(:,2).*c(:,3);
    J(3,1,:,i)=f(:,3).*f(:,1).*c(:,1).*s(:,2).*c(:,3);
    J(3,2,:,i)=f(:,3).*f(:,2).*s(:,1).*c(:,2).*c(:,3);
    J(3,3,:,i)=-f(:,3).^2.*s(:,1).*s(:,2).*s(:,3);

end

end