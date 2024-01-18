function [f,df,y] = gp_rnd_scalar_potential_fast(x,m,LL,theta)
%% GP_RND_SCALAR_POTENTIAL_FAST - Generates random draws from a GP for scalar potential field
%
% Syntax:
%   [f, df, y] = gp_rnd_scalar_potential_fast(x, m, LL, theta)
%
% In:
%   x     - Input points for the GP (Nx3 matrix for 3D points)
%   m     - Number of basis functions
%   LL    - Domain limits for the input points
%   theta - Hyperparameters of the GP ([linSigma2, lengthScale, magnSigma2, sigma2])
%
% Out:
%   f     - GP's function values at the input points (scalar potential)
%   df    - Derivatives of the GP's function values at the input points
%   y     - GP's derivative values with added measurement noise
%
% Description:
%   This function generates random draws from a GP for simulating scalar 
%   potential fields, defined over a specified domain with input points 'x', 
%   number of basis functions 'm', and domain limits 'LL'. The GP's behavior 
%   is governed by hyperparameters specified in 'theta', including linear 
%   sigma, length scale, signal variance, and noise variance. The function 
%   outputs the GP function values, their derivatives, and the derivative 
%   values with added measurement noise.
%
% See also:
%   domain_cartesian_dx
%
% References:
%   [1] Arno Solin and Simo Särkkä (2020). Hilbert space methods for 
%       reduced-rank Gaussian process regression. Statistics and Computing, 
%       30(2):419–446. 
%   [2] Arno Solin, Manon Kok, Niklas Wahlström, Thomas B. Schön, and 
%       Simo Särkkä (2018). Modeling and interpolation of the ambient 
%       magnetic field by Gaussian processes. IEEE Transactions on Robotics 
%       (T-RO), 34(4):1112–1127. 
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

%% Scale inputs

  % Translate x from x \in [Lmin,Lmax] to x \in [-L,L]
  x  = bsxfun(@minus,x,mean(LL,1));
  LL = range(LL)/2;


%% Set up basis functions

  % The eigenbasis
  [eigenval,eigenfun,eigenfun_dx,NN] = domain_cartesian_dx(m,3,LL);
  
  % The eigenvalues
  lambda = eigenval(NN);
  
  % Evaluate Phi for the observations
  Phi   = eigenfun(NN,x);
  dPhix = eigenfun_dx(NN,x,1);
  dPhiy = eigenfun_dx(NN,x,2);
  dPhiz = eigenfun_dx(NN,x,3);
  
  % Account for the linear covariance function
  Phi   = [x Phi]; 
  dPhix = [ones(size(x,1),1) zeros(size(x,1),2) dPhix];
  dPhiy = [zeros(size(x,1),1) ones(size(x,1),1) zeros(size(x,1),1) dPhiy];
  dPhiz = [zeros(size(x,1),2) ones(size(x,1),1) dPhiz];
  dPhi  = [dPhix; dPhiy; dPhiz];
  
  
  
    
%% Set up spectral density of the covariance function 
  
  % The dimensionality of the inputs
  d = size(x,2);

  % The spectral density of the squared exponential covariance function
  Sse = @(w,lengthScale,magnSigma2) ...
            magnSigma2*sqrt(2*pi)^d*lengthScale^d*exp(-w.^2*lengthScale^2/2);   
    
  % Spectral density
  S = @(w,linSigma2,lengthScale,magnSigma2) ...
      [linSigma2; linSigma2; linSigma2; Sse(w(1:end),lengthScale,magnSigma2)];
  
%% Simulate random draw from the GP   
  
  % Extract hyper-parameters
  linSigma2   = theta(1);
  lengthScale = theta(2);
  magnSigma2  = theta(3);
  sigma2      = theta(4);
  
  % Solve GP with optimized hyperparameters and 
  % return predictive mean and variance
  k = S(sqrt(lambda),linSigma2,lengthScale,magnSigma2);
  foo = diag(sqrt(k))*randn(numel(k),1);
  f = Phi*foo;
  df = [dPhix*foo dPhiy*foo dPhiz*foo];
  
  % Measurement noise
  y = df + sqrt(sigma2)*randn(size(df));


