function [Eft,dEft,Varft,theta,lik,dVarft] = gp_scalar_potential_fast(x,y,xt,m,LL,theta,opt,predcf)
%% gp_scalar_potential_fast - Scalar potential GP (reduced-rank)
%
% Syntax:
%   [Eft,dEft,Varft,theta,lik,dVarft] = gp_scalar_potential_fast(x,y,xt,m,LL,theta,opt,predcf)
%
% In:
%   x      - Training inputs (n x 3)
%   y      - Training target gradients (n x 3)
%   xt     - Test inputs (n x 3)
%   m      - Number of basis functions
%   LL     - Domain boundaries [x1_min x2_min x3_min; x1_max x2_max x3_max]
%   theta  - Hyperparameters (linSigma2,lengthScale,magnSigma2,sigma2)
%   opt    - Optimize hyperparamters (default: true)
%   predcf - The covariance functions used in prediction (default: [1 2])
%
% Out:
%   Eft    - Test target posterior mean
%   dEft   - Test target gradient posterior mean
%   Varft  - Test target posterior marginal variance
%   theta  - ML estimate of optimized hyperparameters
%   lik    - Negative log marginal likelihood
%   dVarft - Test target gradient posterior marginal variance
%   
% Description:
%   Minimum working example of the proposed method.
%   The model is as follows:
%
%     f(x) ~ GP(0,k_lin.(x,x') + k_SE(x,x'))
%     y_i = \nabla f(x_i) + \epsilon_i, 
%
%   This code follows the presentation in the notes, and is not optimized
%   for speed. For more exotic domains or other covariance functions, only
%   the spectral density and functions eigenval/eigenfun need to be
%   replaced.
%
% Copyright: 
%   (c) Arno Solin, 2015
% 

%% Check inputs and set up model

  % Domain boundaries
  if nargin < 5 || isempty(LL)
    rangex = range(x);
    pm = 0.1*min(rangex(rangex>0));
    LL = [min(x)-pm; max(x)+pm];
  end
  
  % Initial hyper-parameters
  if nargin < 6 || isempty(theta)
    theta = [800 0.1*mean(range(x)) 1 1];
  end

  % Optimize hyperparameters
  if nargin < 7 || isempty(opt)
    opt = [true true true true];
  end

  % Show only contribution of the SE field
  if nargin < 8 || isempty(predcf)
    predcf = [1 2];
  end

  % Optimize all parameters by default
  if numel(opt)==1
    opt = [opt opt opt opt];
  end
  
%% Scale inputs

  % Translate x from x \in [Lmin,Lmax] to x \in [-L,L]
  x  = bsxfun(@minus,x,mean(LL,1));
  xt = bsxfun(@minus,xt,mean(LL,1));
  LL = range(LL)/2;
  
  
%% Set up basis functions

  % The eigenbasis
  [eigenval,eigenfun,eigenfun_dx,NN] = domain_cartesian_dx(m,3,LL);
  
  % The eigenvalues
  lambda = eigenval(NN);
  
  % Evaluate Phi for the observations
  dPhix = eigenfun_dx(NN,x,1);
  dPhiy = eigenfun_dx(NN,x,2);
  dPhiz = eigenfun_dx(NN,x,3);
  
  % Evaluate the Phi for test inputs
  Phit   = eigenfun(NN,xt);
  dPhixt = eigenfun_dx(NN,xt,1);
  dPhiyt = eigenfun_dx(NN,xt,2);
  dPhizt = eigenfun_dx(NN,xt,3);
  
  % Account for the linear covariance function
  dPhix  = [ones(size(x,1),1) zeros(size(x,1),2) dPhix];
  dPhiy  = [zeros(size(x,1),1) ones(size(x,1),1) zeros(size(x,1),1) dPhiy];
  dPhiz  = [zeros(size(x,1),2) ones(size(x,1),1) dPhiz];
  
  % Which components to predict
  Phit   = [any(predcf==1)*xt any(predcf==2)*Phit];
  dPhixt = [any(predcf==1)*[ones(size(xt,1),1) zeros(size(xt,1),2)]                     any(predcf==2)*dPhixt];  
  dPhiyt = [any(predcf==1)*[zeros(size(xt,1),1) ones(size(xt,1),1) zeros(size(xt,1),1)] any(predcf==2)*dPhiyt];  
  dPhizt = [any(predcf==1)*[zeros(size(xt,1),2) ones(size(xt,1),1)]                     any(predcf==2)*dPhizt];  
  
  
%% Set up spectral density of the covariance function 
  
  % The dimensionality of the inputs
  d = size(x,2);

  % The spectral density of the squared exponential covariance function
  Sse = @(w,lengthScale,magnSigma2) ...
            magnSigma2*sqrt(2*pi)^d*lengthScale^d*exp(-w.^2*lengthScale^2/2);   
  
  % Derivative w.r.t. lengthScale  
  dSse{1} = @(w,lengthScale,magnSigma2) ...
            (d/lengthScale - lengthScale*w.^2).* ...
            Sse(w,lengthScale,magnSigma2);
  
  % Derivative w.r.t. magnSigma2  
  dSse{2} = @(w,lengthScale,magnSigma2) ...
            Sse(w,lengthScale,magnSigma2)/magnSigma2;
  
  % Spectral density
  S = @(w,linSigma2,lengthScale,magnSigma2) ...
      [linSigma2; linSigma2; linSigma2; Sse(w(1:end),lengthScale,magnSigma2)];
  dS{1} = @(w,linSigma2,lengthScale,magnSigma2) ...
      [1; 1; 1; 0*w];
  dS{2} = @(w,linSigma2,lengthScale,magnSigma2) ...
      [0; 0; 0; dSse{1}(w,lengthScale,magnSigma2)];
  dS{3} = @(w,linSigma2,lengthScale,magnSigma2) ...
      [0; 0; 0; dSse{2}(w,lengthScale,magnSigma2)];
  
  % Pre-calculate basis functions
  Phi    = [dPhix; dPhiy; dPhiz];
  PhiPhi = Phi'*Phi;        % O(nm^2)
  Phiy   = Phi'*y(:);
 

%% Optimize hyper-parameters  
  
  if any(opt==true)

    % Optimize hyperparameters
    w0 = log(theta);
    fun = @(w) opthyperparams(w,y(:),lambda,Phiy,PhiPhi,S,dS,theta,opt);
    options = optimset('Display','iter','GradObj','on', ...
                       'TolX',1e-6,'TolFun',1e-3,'DerivativeCheck','off');
  
    % Warning off
    warning off
                   
    % Optimize
    [w,lik,~,output] = fminunc(fun,w0(opt),options);

    % Warning on
    warning on
    
    % The paramters
    theta(opt) = exp(w);
  
  else
    
    % Evaluate likelihood (neg.)
    lik = opthyperparams(log(theta),y(:),lambda,Phiy,PhiPhi,S,dS,theta,true(1,4));
    
  end
  
  
%% Solve the GP regression problem
  
  % Extract hyper-parameters
  linSigma2   = theta(1);
  lengthScale = theta(2);
  magnSigma2  = theta(3);
  sigma2      = theta(4);
  
  % Report hyper-parameters
  fprintf(['Hyper-parameters:\n' ...
      '  linSigma2   : %.4f\n' ...
      '  lengthScale : %.4f\n' ...
      '  magnSigma2  : %.4f\n' ...
      '  sigma2      : %.4f\n'],theta)
  
  % Solve GP with optimized hyperparameters and 
  % return predictive mean and variance
  k = S(sqrt(lambda),linSigma2,lengthScale,magnSigma2);
  L = chol(PhiPhi + diag(sigma2./k),'lower');
  foo = (L'\(L\Phiy));
  Eft = Phit*foo;
  Eftx = dPhixt*foo;
  Efty = dPhiyt*foo;
  Eftz = dPhizt*foo;
  Varft = sigma2*sum((Phit/L').^2,2);

  % Also return variance for components
  if nargout>5
    dVarft = sigma2*[sum((dPhixt/L').^2,2) ...
                     sum((dPhiyt/L').^2,2) ...
                     sum((dPhizt/L').^2,2)];  
  end
  
  % Concatenate
  dEft = [Eftx Efty Eftz];

  % Report model  
  fprintf('The linear model is as follows: (%.3f, %.3f, %.3f). \n', ...
      foo(1:3))
  
end
  
function [e,eg] = opthyperparams(w,y,lambda,Phiy,PhiPhi,S,dS,theta,opt)

  % Initialize
  e  = nan;
  eg = nan(1,numel(theta));
  
  % Extract parameters
  theta(opt)  = exp(w);
  linSigma2   = theta(1); % Linear model scale parameter
  lengthScale = theta(2); % Length scale parameter 
  magnSigma2  = theta(3); % Magnitude parameter
  sigma2      = theta(4); % Measurement noise variance
  
  % Evaluate the spectral density
  k = S(sqrt(lambda),linSigma2,lengthScale,magnSigma2);
  
  % Number of n=observations and m=basis functions
  n = numel(y);
  m = size(Phiy,1);
  
  % Calculate the Cholesky factor
  [L,p] = chol(PhiPhi + diag(sigma2./k),'lower');  % O(m^3)
  
  % Check if pos. def
  if p>0, return; end
  
  % Evaluate all parts
  v = L\Phiy; % Phiy = (Phi'*y);
  yiQy = (y'*y - v'*v)/sigma2;
  logdetQ = (n-m)*log(sigma2) + sum(log(k)) + 2*sum(log(diag(L)));
  
  % Return approx. negative log marginal likelihood
  e = .5*yiQy + .5*logdetQ + .5*n*log(2*pi);
  
  % Precalculate
  vv = L'\v;
  LLk = L'\(L\diag(1./k)); % O(m^3)
  
  % Return if no derivatives requested
  if nargout==1 || isnan(e), return; end
  
  % For the covariance function hyperparameters
  for j=1:numel(theta)-1
    
    % Should we skip this?
    if opt(j)==false, continue; end
    
    % Evaluate the partial derivative
    dk = dS{j}(sqrt(lambda),linSigma2,lengthScale,magnSigma2);
        
    % Evaluate parts
    dlogdetQ = sum(dk./k) - sigma2*sum(diag(LLk).*dk./k);
    dyiQy = -vv'*diag(dk./k.^2)*vv;
        
    % The partial derivative
    eg(j) = .5*dlogdetQ + .5*dyiQy;
        
    % Account for the log-transformed values
    %eg(j) = exp(w(j))*eg(j);
    
  end
  
  % Gradient of noise magnitude, sigma2
  if opt(end)==true
  
    % Evaluate parts
    dlogdetQ = (n-m)/sigma2 + sum(diag(LLk));
    dyiQy    = vv'*diag(1./k)*vv/sigma2 - yiQy/sigma2;
  
    % For the measurement noise
    eg(end)  = .5*dlogdetQ + .5*dyiQy;
    
    % Account for the log-transformed values
    %eg(end) = exp(w(end))*eg(end);

  end
  
  % Remove those gradient that was not calculated
  eg(isnan(eg)) = [];
  
  % Account for the log-transformed values
  eg(:) = exp(w(:)).*eg(:);
  
end