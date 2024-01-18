function [f,y] = gp_rnd_SE1D_fast(x,m,LL,theta)
%% GP_RND_SE1D_FAST - Random draws from a GP with a SE covariance function
%
% Syntax:
%   [f, y] = gp_rnd_SE1D_fast(x, m, LL, theta)
%
% In:
%   x     - Input points for the GP (Nx1 vector)
%   m     - Number of basis functions
%   LL    - Domain limits for the input points
%   theta - Hyperparameters of the GP ([lengthScale, magnSigma2, sigma2])
%
% Out:
%   f     - GP's function values at the input points
%   y     - GP's function values with added measurement noise
%
% Description:
%   This function generates random draws from a GP with a squared exponential 
%   covariance function. The GP is defined over a specified domain with input 
%   points 'x', number of basis functions 'm', and domain limits 'LL'. The 
%   GP's behavior is governed by hyperparameters specified in 'theta', 
%   including length scale, signal variance, and noise variance. The function 
%   outputs both the noise-free GP function values and the values with added 
%   measurement noise.
%
%   Note: The '1D' in the function name referes to single-output and 
%   not input dimensionality (can be misleading).
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

  % The dimensionality of the inputs
  d = size(x,2);

  % The eigenbasis
  [eigenval,eigenfun,~,NN] = domain_cartesian_dx(m,d,LL);
  
  % The eigenvalues
  lambda = eigenval(NN);
  
  % Evaluate Phi for the observations
  Phi   = eigenfun(NN,x);  
    
%% Set up spectral density of the covariance function 

  % The spectral density of the squared exponential covariance function
  S = @(w,lengthScale,magnSigma2) ...
            magnSigma2*sqrt(2*pi)^d*lengthScale^d*exp(-w.^2*lengthScale^2/2);   
  
%% Simulate random draw from the GP   
  
  % Extract hyper-parameters
  lengthScale = theta(1);
  magnSigma2  = theta(2);
  sigma2      = theta(3);
  
  % Solve GP with optimized hyperparameters and 
  % return predictive mean and variance
  k = S(sqrt(lambda),lengthScale,magnSigma2);
  foo = diag(sqrt(k))*randn(numel(k),1);
  f = Phi*foo;
  
  % Measurement noise
  y = f + sqrt(sigma2)*randn(size(f));


