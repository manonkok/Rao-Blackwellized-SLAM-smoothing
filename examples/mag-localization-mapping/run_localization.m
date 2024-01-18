function run_localization(data,saveFrames,saveVideo)
%% run_localization - Run magentic terrain matching example
%
% See main.m for data loading.
%
% References:
%   [1] Arno Solin, Simo Särkkä, Juho Kannala, and Esa Rahtu (2016). 
%       Terrain navigation in the magnetic landscape: Particle filtering 
%       for indoor positioning. In Proceedings of the European Navigation 
%       Conference (ENC). Pages 1–9. IEEE.
%   [2] Manon Kok and Arno Solin (2018). Scalable magnetic field SLAM in 
%       3D using Gaussian process maps. In Proceedings of the International 
%       Conference on Information Fusion (FUSION). Pages 1353–1360. 
%       Cambridge, UK.
%
%% Parameters for particle filtering

  addpath(genpath('../../'))

  % Set options here
  showVisualization = true;
  %saveVideo = false;
  %saveFrames = false;
  showRobotOverlay = false;
  sameRndSeed = true;
  
  N_P = 1000; % Number of particles
  Q = blkdiag(diag(4^2*[0.01^2,0.01^2,0.01^2]), diag([1e-2 1e-2 1e-2]*pi/180).^2);
  params.Q = Q;
  theta = [500.0000    0.1178  384.6590    3.5859];
  params.theta = theta;
    
  % Define indices
  iPos = 1:3;                         % Indices pos in states
  iQuat = 4:7;                        % Indices quat in states

  % Always same seed
  if sameRndSeed
    rng(1,'twister')
  end
  
  
%% Data

  % Training data for magnetic field
  ind = data.s<3 | data.s==4;
  x = data.x(ind,:);
  y = data.y(ind,:);
  x = [x 0*x(:,1)];  
  
  % Thinning
  x = x(1:10:end,:);
  y = y(1:10:end,:);
  
  % Data along the observed path
  ind = data.s==3;
  xp = data.x(ind,:);
  yp = data.y(ind,:);
  xp = [xp 0*xp(:,1)];  
  
  % Thinning
  xp = xp(1:50:end,:);
  yp = yp(1:50:end,:);
  
  % Control signal and time steps
  dx = diff(xp);
  dt = 0.1;
  
  % Orientation (angle -> rotm -> quat)
  psi = atan2(dx(:,2),dx(:,1)); psi = [psi; psi(end)];
  N = numel(psi);
  R = [reshape(cos(psi),[1 1 N]), reshape(sin(psi),[1 1 N]), zeros(1,1,N) ; ...
       reshape(-sin(psi),[1 1 N]), reshape(cos(psi),[1 1 N]), zeros(1,1,N) ; ...
       zeros(1,1,N), zeros(1,1,N), ones(1,1,N)];
  %quat = rotm2quat(R); % using robottics toolbox 
  quat = rmat2quat(R)';

  % Increments (~control signal)
  initState0 = [xp(1,:)' ; quat(1,:)'];
  dPos = diff(xp);
  dQuat = squeeze(multiprod(qLeft(qInv(quat(1:end-1,:))), ...
            reshape(quat(2:end,:)',[4 1 length(quat)-1])))';
  u = [dPos dQuat];

  % Rotate magnetic field
  for i=1:size(yp,1)
    yp(i,:) = (R(:,:,i)*yp(i,:)')';
  end
  
  
%% Preliminaries on GP model

  % Linear states: basis functions
  nBasisFunctions = 1000;              
  d = 3;
  
  % Domain boundaries
  rangex = range(x);
  pm = 0.2*min(rangex(rangex>0));
  LL = [min(x)-pm; max(x)+pm];
    
  % Test points
  x1t = linspace(LL(1,1),LL(2,1),100);
  x2t = linspace(LL(1,2),LL(2,2),100);
  [X1t,X2t]=meshgrid(x1t,x2t);
  xt = [X1t(:) X2t(:) 0*X1t(:)];  
  
  % Shift origin
  centerpoint = mean(LL,1);
  x  = bsxfun(@minus,x,mean(LL,1));
  xt = bsxfun(@minus,xt,mean(LL,1));
  xp  = bsxfun(@minus,xp,mean(LL,1));  
  LL = range(LL)/2;
    
  % Eigenbasis
  [eigenval,~,eigenfun_dx,NN] = domain_cartesian_dx(nBasisFunctions,d,LL);

  % The eigenvalues
  lambda = eigenval(NN);

  % Hyperparameters 
  linSigma2   = theta(1); % Linear model scale parameter
  lengthScale = theta(2); % Length scale parameter 
  magnSigma2  = theta(3); % Magnitude parameter
  sigma2 = theta(4);

  % Evaluate the spectral density
  Sse = @(w,lengthScale,magnSigma2) ...
    magnSigma2*sqrt(2*pi)^d*lengthScale^d*exp(-w.^2*lengthScale^2/2);   
  S = @(w,linSigma2,lengthScale,magnSigma2) ...
    [linSigma2; linSigma2; linSigma2; Sse(w(1:end),lengthScale,magnSigma2)];
  k = S(sqrt(lambda),linSigma2,lengthScale,magnSigma2);

  % Evaluate Phi for the observations
  dPhix = eigenfun_dx(NN,x,1);
  dPhiy = eigenfun_dx(NN,x,2);
  dPhiz = eigenfun_dx(NN,x,3);
  
  % Account for the linear covariance function
  dPhix  = [ones(size(x,1),1) zeros(size(x,1),2) dPhix];
  dPhiy  = [zeros(size(x,1),1) ones(size(x,1),1) zeros(size(x,1),1) dPhiy];
  dPhiz  = [zeros(size(x,1),2) ones(size(x,1),1) dPhiz];
  
  % Pre-calculate basis functions
  Phi    = [dPhix; dPhiy; dPhiz];
  PhiPhi = Phi'*Phi;        % O(nm^2)
  Phiy   = Phi'*y(:);
  
  % Pre-calculate predictions for weighting
  L = chol(PhiPhi + diag(sigma2./k),'lower');
  foo = (L'\(L\Phiy));
  

%% Initial state

 initState = initState0 * ones(1,N_P);
 rx = range(x);
 mx = min(x);
 initState(1,:) = mx(1) + rx(1)*rand(1,N_P);
 initState(2,:) = mx(2) + rx(2)*rand(1,N_P);
 initState(iQuat,:) = initState(iQuat,:);
  
 
%% Domain for visualization
  
  nt = length(xt);
  dPhixt = eigenfun_dx(NN,xt,1);
  dPhiyt = eigenfun_dx(NN,xt,2);
  dPhizt = eigenfun_dx(NN,xt,3);
  dPhixt = [ones(nt,1), zeros(nt,2), dPhixt];
  dPhiyt = [zeros(nt,1), ones(nt,1), zeros(nt,1), dPhiyt];
  dPhizt = [zeros(nt,2), ones(nt,1), dPhizt];

  
%% For homography overlay visualization

  % Corners
  corner = [1 1; -1 1; -1 -1; 1 -1];

  % Markers
  markers = ...
  [ 4.7741   -3.9192    0.1101
    4.9283    1.2333    0.0896
   -1.2583    1.0878   -0.0092
   -1.2033   -2.7739    0.0136]; 
  
  % Markers on image
  im_markers = [1314,130; 279,213; 139,1017; 1798,749];

  % Grid
  z1 = [linspace(-1,1,32) nan];
  z2 = [linspace(-1,1,32) nan];
  [G1,G2] = meshgrid(z1,z2);
  Z = [G1(:) G2(:)];
  
  % Mapping
  [A,C] = homography_estimation(markers(:,1:2),corner);

  % Transform
  G = bsxfun(@rdivide,A*[Z ones(size(Z,1),1)]',C*[Z ones(size(Z,1),1)]')';
    
  % Mapping
  [A,C] = homography_estimation(im_markers,markers(:,1:2));
    
  % Show grid
  Y = bsxfun(@rdivide,A*[G(:,1:2) ones(size(G,1),1)]',C*[G(:,1:2) ones(size(G,1),1)]')';
  G1 = reshape(Y(:,1),size(G1));
  G2 = reshape(Y(:,2),size(G2));
  
  
%% Initial state and covariance, measurement noise and process noise covariance

  % Initial state
  x0_nonLin = initState;

  % Measurement noise
  ny = 3;                      % Measurement dimension
  R = sigma2 * eye(3);


%% Run filter

  % Prepare the new file
  if saveVideo == true
    vidObj = VideoWriter('robot-pf.mp4','MPEG-4');
    vidObj.FrameRate = 10;
    open(vidObj);
  end
  
  % Run filter
  [traj_max,traj_mean] = particleFilterLocalization(@dynModel,@measModel,u,yp, ...
    x0_nonLin,Q,R,N_P,dt,@makePlotsFilter_dense3D_magfield);
  
  % Close video
  if saveVideo == true
    close(vidObj);
  end
  
  
%% Supporting functions
function w = measModel(yt,xn)
    %keyboard
    
    N_pred = size(xn,2);
    dPhix = eigenfun_dx(NN,xn(iPos,:)',1);
    dPhiy = eigenfun_dx(NN,xn(iPos,:)',2);
    dPhiz = eigenfun_dx(NN,xn(iPos,:)',3);
    dPhix = [ones(N_pred,1), zeros(N_pred,2), dPhix];
    dPhiy = [zeros(N_pred,1), ones(N_pred,1), zeros(N_pred,1), dPhiy];
    dPhiz = [zeros(N_pred,2), ones(N_pred,1), dPhiz];
    Rnb = quat2rmat(xn(iQuat,:)');
    dy = zeros(N_pred,ny,size(dPhix,2));
    for i = 1:N_pred
        dy(i,:,:) = Rnb(:,:,i)' * ...
            [dPhix(i,:) ; dPhiy(i,:) ; dPhiz(i,:)];
    end
    
  % Solve GP with optimized hyperparameters and 
  % return predictive mean and variance
  dEft = [dPhix*foo dPhiy*foo dPhiz*foo];
  dVarft = sigma2*[sum((dPhixt/L').^2,2) ...
                   sum((dPhiyt/L').^2,2) ...
                   sum((dPhizt/L').^2,2)];  
  % Rotations
  Rnb = quat2rmat(xn(iQuat,:)');
              
  % Return weights
  w = nan(1,N_pred);
  for i=1:N_pred
    w(i) = sum(normpdf(yt,(Rnb(:,:,i)'*dEft(i,:)')',sqrt(dVarft(i,:)+sigma2)));  % + sigma2
  end
end

function [xpred, dQuat] = dynModel(xn,dx,dt,Q)
    % Predict through dynamic model. Also optionally output dQuat for
    % generating odometry data
    xpred_pos = xn(iPos) + dx(iPos)' + sqrt(dt * Q(iPos,iPos))*randn(3,1);
    xpred_quat = qLeft(qRight(xn(iQuat)) * dx(iQuat)') * ...
            expq(sqrt(dt * Q(4:6,4:6))*randn(3,1));
    xpred = [xpred_pos ; xpred_quat]; 
end

function makePlotsFilter_dense3D_magfield(xn,traj_max,yhattraj,xn_traj,traj_mean)

  % Skip if no visualization to show
  if ~showVisualization
      fprintf('Running step %i/%i\n',sum(~isnan(traj_max(1,:))),size(traj_max,2))
      return; 
  end   
    
    %{
    % Simple visualization
    figure(2); cla; hold on
      plot(x(:,1),x(:,2),'.g')    
      plot(xn(1,:),xn(2,:),'.k')
      t = sum(~isnan(traj_max(1,:)));
      %plot(traj_max(1,1:t),traj_max(2,1:t),'k')
      plot(xp(1:t,1),xp(1:t,2),'-k')
      %plot(traj_max(1,1:t),traj_max(2,1:t),'-r')
      %plot(traj_mean(1,1:t),traj_mean(2,1:t),'-b')
      axis equal
      xlim([min(xt(:,1)) max(xt(:,1))])
      ylim([min(xt(:,2)) max(xt(:,2))]) 
      box on
      drawnow
      
    %}
    
  % Frame index
  t = sum(~isnan(traj_max(1,:)));

  % Load background
  if showRobotOverlay
    [bg,img,mask] = mask_for_frame(t);
  else
    bg = imread('frame.png');
  end
  
  % Transform particle cloud
  Z = xn(1:2,:)'+centerpoint(1:2);
  X = bsxfun(@rdivide,A*[Z ones(size(Z,1),1)]',C*[Z ones(size(Z,1),1)]')';

  % Transform ground-truth
  t = sum(~isnan(traj_max(1,:)));
  Z = xp(1:t,1:2)+centerpoint(1:2);
  path = bsxfun(@rdivide,A*[Z ones(size(Z,1),1)]',C*[Z ones(size(Z,1),1)]')';
  
  % Transform path mean
  Z = traj_mean(1:2,1:t)'+centerpoint(1:2);
  path_mean = bsxfun(@rdivide,A*[Z ones(size(Z,1),1)]',C*[Z ones(size(Z,1),1)]')';

  % Transform path max
  Z = traj_max(1:2,1:t)'+centerpoint(1:2);
  path_max = bsxfun(@rdivide,A*[Z ones(size(Z,1),1)]',C*[Z ones(size(Z,1),1)]')';
    
  % Create figure
  figure(1); cla; hold on

  % Show background
  imagesc(bg)
  
  % Show grid
  plot(G1,G2,'-',G1',G2','-','Color',[.7 .7 .7],'LineWidth',.5)

  % Show ground-truth path
  plot(path(:,1),path(:,2),'--k','LineWidth',1)
 
  % Show mean
  plot(path_mean(:,1),path_mean(:,2),'-b','LineWidth',1)

  % Show max
  %plot(path_max(:,1),path_max(:,2),'--r','LineWidth',1)
    
  % Plot particles
  scatter(X(:,1),X(:,2),45,'filled','MarkerFaceAlpha',.1,'cdata',[0 0 1]);
  
  % Axis
  axis ij image off
  xlim([0 1920]), ylim([0 1080])

  % Show robot on top
  if showRobotOverlay
    h=imagesc(img);
    mask = .8*double(mask);
    %mask(1) = 1;
    set(h,'AlphaData',mask,'AlphaDataMapping','scaled')
  end
  
  % BG Color
  set(gcf,'Color','w')
  
  drawnow
  pause(.1)
 
  % Save as video
  if saveVideo == true  
    currFrame = getframe(gca);
    writeVideo(vidObj,currFrame);
  end
  
  % Save frame images every 10th frame
  if saveFrames == true && rem(t,10) == 1
    frame = getframe(gca);
    imwrite(frame.cdata,sprintf('robot-frame-%03d.png',t)) 
  end
      
end

end