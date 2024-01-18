function [eigenval,eigenfun,eigenfun_dx,NN] = domain_cartesian_dx(m,d,L)
%% domain_cartesian_dx - Laplace operator eigendecomposition in a hypercube
% 
% Syntax:
%  [eigenfun,eigenval,NN] = domain_cartesian(m,d,L)
%
% In:
%   m  - Number of eigenfunctions
%   d  - Dimensionality
%   L  - Domain boundary [-L1,L1]x[-L2,L2]x...x[-Ln,Ln]
%      
% Out:
%   eigenval    - Function handle: eigenval(n)
%   eigenfun    - Function handle: eigenfun(n,x)
%   eigenfun_dx - Function handle: eigenfun_dx(n,x)
%   NN          - Indices to evaluate the handles at
% 
% Description:
%   This code returns the eigendecomposition of the Laplacian in
%   Cartesian coordinates (x1,x2,...) in [-L1,L1]x[-L2,L2]x...
%   with respect to indices (n1,n2,...). The function is vectorized
%   with respect to both the location x and the indices n.
%
% Copyright (c) 2014-2023 Arno Solin

  % If domain is not centered, center it
  if size(L,1)>1
    L = (max(L,[],1)-min(L,[],1))/2;
  end

  % This is stupid, but at least we should get enough 
  % of basis function to choose from
  N = ceil(m^(1/d)*L/min(L));
  
  % Combined eigenfunction indices (checked the numbers)
  NN = ndgridm(N);

  % Define eigenvalues of the negative Laplacian in ND 
  % s.t. Dirichlet boundary. This forms an orthonormal basis.
  eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
    
  % Sort and take only m most important eigenfunctions
  [~,ind] = sort(eigenval(NN)); NN = NN(ind(1:m),:);  

  % Define eigenfunction of the negative Laplacian in ND 
  % s.t. Dirichlet boundary. This forms an orthonormal basis.
  eigenfun = @(n,x) laplace_eig_cart_dirichlet(n,x,L);

  % Define derivative of the eigenfunction of the negative 
  % Laplacian in ND s.t. Dirichlet boundary 
  eigenfun_dx = @(n,x,di) laplace_eig_cart_dirichlet_dx(n,x,di,L);
  
  
end


function [v]=laplace_eig_cart_dirichlet(n,x,L)
%% laplace_eig_cart_dirichlet - Laplace operator eigenfunctions in a hypercube
% 
% Syntax:
%  [v] = laplace_eig_cart_dirichlet(n,x,L)
%
% In:
%   n  - Eigenfunction indices
%   x  - Spatial locations [x1 x2]
%   L  - Domain boundary [-L1,L1]x[-L2,L2]x...x[-Ln,Ln]
%      
% Out:
%   v - The evaluated value
% 
% Description:
%   This code calculates the eigenvectors of the Laplacian in
%   Cartesian coordinates (x1,x2,...) in [-L1,L1]x[-L2,L2]x...
%   with respect to indices (n1,n2,...). The function is vectorized
%   with respect to both the location x and the indices n.
%
%   The corresponding eigenvalues can be calculated by
% 
%     eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
%
% Copyright (C) 2012 Arno Solin
%

  % Allocate space
  v = ones(size(x,1),size(n,1));

  % Evaluate eigenfunctions
  for j=1:size(n,2)
    for i=1:size(n,1)
      v(:,i) = v(:,i) .* 1./sqrt(L(j)) .* ...
        sin(pi*n(i,j).*(x(:,j)+L(j))/(2*L(j)));
    end
  end
  
  
%   % Evaluate eigenfunctions
%   if size(x,2)==1
%       for i=1:numel(n)
%           v(:,i) = sqrt(1./L)*sin(pi*n(i)*(x(:)+L)/2/L);
%       end
%   else
%       for i=1:size(n,1)
%           % Eigenfunctions for x in Omega and n = (n1,n2,...nn), ni = 1,2,...,Ni
%           v(:,i) = prod(bsxfun(@times,sqrt(1./L), ...
%               sin(pi*bsxfun(@times,n(i,:)./L,bsxfun(@plus,x,L))/2)),2);
%           if all(n(i,:)==0)
%               v(:,i) = ones(size(x,1),1);
%           end
%       end
%   end

end

function [v]=laplace_eig_cart_dirichlet_dx(n,x,di,L)
%% laplace_eig_cart_dirichlet_dx - Derivative of Laplace operator eigenfunctions in a hypercube
% 
% Syntax:
%  [v] = laplace_eig_cart_dirichlet_dx(n,x,di,L)
%
% In:
%   n  - Eigenfunction indices
%   x  - Spatial locations [x1 x2]
%   di - Differentiate w.r.t. this dimension d/dx_i
%   L  - Domain boundary [-L1,L1]x[-L2,L2]x...x[-Ln,Ln]
%      
% Out:
%   v - The evaluated value
% 
% Description:
%   This code calculates the eigenvectors of the Laplacian in
%   Cartesian coordinates (x1,x2,...) in [-L1,L1]x[-L2,L2]x...
%   with respect to indices (n1,n2,...). The function is vectorized
%   with respect to both the location x and the indices n.
%
%   The corresponding eigenvalues can be calculated by
% 
%     eigenval = @(n) sum((pi*bsxfun(@rdivide,n,2*L)).^2,2);
%
% Copyright (C) 2015 Arno Solin
%

  % Allocate space
  v = ones(size(x,1),size(n,1));

  % Evaluate eigenfunctions
  for j=1:size(n,2)
    if j==di % Differentiate w.r.t this dimension
      for i=1:size(n,1)
        %v(:,i) = v(:,i) .* 1./sqrt(L(j)) .* ...
        %  cos(pi*n(i,j).*(x(:,j)+L(j))/(2*L(j)));
        
        % This is the actual derivative  
        v(:,i) = v(:,i) .* pi*n(i,j)/(2*L(j)*sqrt(L(j))) .* ...
          cos(pi*n(i,j).*(x(:,j)+L(j))/(2*L(j)));
      
        %v(:,i) = v(:,i).*prod(bsxfun(@times,sqrt(1./L(j)), ...
        %   cos(pi*bsxfun(@times,n(i,j)./L(j),bsxfun(@plus,x(:,j),L(j)))/2)),2);
      end
    else % Do not differentiate w.r.t this dimension
      for i=1:size(n,1)
        v(:,i) = v(:,i) .* 1./sqrt(L(j)) .* ...
          sin(pi*n(i,j).*(x(:,j)+L(j))/(2*L(j)));
        %v(:,i) = v(:,i).*prod(bsxfun(@times,sqrt(1./L(j)), ...
        %   sin(pi*bsxfun(@times,n(i,j)./L(j),bsxfun(@plus,x(:,j),L(j)))/2)),2);       
        %if all(n(i,j)==0)
        %  %v(:,i) = ones(size(x,1),1);
        %end
      end
    end
  end
  
end

function NN = ndgridm(N)
%% ndgridm - Expand index hypercude
%
% Syntax:
%  [NN] = ndgridm(N)
%
% In:
%   N  - Vector of max indices
%      
% Out:
%   NN - Matrix of index combinations
%
% Description:
%   A more felxible variant of 'ndgrid'. This functions gives combinations
%   of the indices in N such that, for example, for a 2D case we get
%   (1,1),(1,2),...,(1,N2),(2,1),...,(2,N2),(...),(N1,1),...,(N1,N2).
%   This function works for any dimension.
%
% Copyright (C) 2014 Arno Solin
%

  % Allocate space for indices
  NN = zeros(prod(N),numel(N));

  % For each level/diemsion
  if numel(N)==1

     % The lowest level
     NN(:,1) = (1:N)';
     
  else

    % This level
    n = 1:N(1);

    % Recursive call
    nn = ndgridm(N(2:end));

    % Assign values
    NN(:,1)     = kron(n,ones(1,prod(N(2:end))))';
    NN(:,2:end) = repmat(nn,[N(1) 1]);
   
  end

end
