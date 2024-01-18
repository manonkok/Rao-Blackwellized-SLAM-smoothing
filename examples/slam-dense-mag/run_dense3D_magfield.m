function [savedData,groundTruth] = run_dense3D_magfield(params)
%% RUN_DENSE3D_MAGFIELD - Runs a dense 3D magnetic field experiments
%
% Syntax:
%   [savedData, groundTruth] = run_dense3D_magfield(params)
%
% In:
%   params - Struct containing various parameters for the simulation
%
% Out:
%   savedData    - Data saved from the simulation
%   groundTruth  - Struct containing ground truth data for the simulation
%
% Description:
%   This function runs dense 3D magnetic field SLAM based on the 
%   provided parameters. It generates simulation data, including odometry 
%   and sensor measurements, and can run with default settings or custom 
%   parameters provided in 'params'. The function also allows for 
%   visualization and uses a specific random seed if run separately. It's 
%   suited for trajectories like 'bean_6D' and others, with specific noise 
%   and Gaussian Process (GP) settings.
%
% See also:
%   generateData_dense, dynModel
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

% Parameters should either be input via params struct or default parameters
% are used 
if nargin < 1
    % If file run separately then we want to always use the same seed and
    % make plots
    sameRndSeed = 1;
    makePlots = 1;
    params.visualiseResults = 1;
    N_K = 3; % Number of smoother iterations
    % Default is bean trajectory in 6D ... 
    params.trajType = 'bean_6D';
    % ... with some specific process noise and GP settings
    % Q = blkdiag(diag(5^2*[0.05^2,0.05^2,0.01^2]), diag([0.01 0.01 0.12]*pi/180).^2);            
    % Q = blkdiag(diag(30^2*[0.05^2,0.05^2,0.01^2]), diag([0.01 0.01 5*1]*pi/180).^2);
    % Q = blkdiag(diag(2^2*[0.05^2,0.05^2,0.01^2]), diag([0.01 0.01 2*0.12]*pi/180).^2);
    Q = blkdiag(diag(8^2*[0.05^2,0.05^2,0.01^2]), diag([0.01 0.01 0.25]*pi/180).^2);
    params.Q = Q;
    theta = [628.5196 ; 1.3286 ; 162.0930 ; 5.9900];
    % theta = [628.5196 ; 1.3286 ; 162.0930 ; 3];
%     theta = [628.5196 ; 0.5 ; 162.0930 ; 3];
    params.theta = theta;
    params.magDisti = zeros(1,3);
else
    % If parameters input via params struct ... 
    % ... don't use the same seed and don't visualise any intermediate
    % results
    sameRndSeed = 0;
    makePlots = 0;
    params.visualiseResults = 0;
    % Extract info from params struct
    Q = params.Q;
    theta = params.theta;
    N_K = params.N_K;
end

%% Generate data and define details about model
if sameRndSeed
    rng(1,'twister')
end

% Define indices and measurement dimension
iPos = 1:3;     % Indices pos in states
iQuat = 4:7;    % Indices quat in states
ny = 3;         % Measurement dimension

% Define sampling time
dt = 0.01;
params.dt = dt;

params.makePlots = makePlots;
[dx, initState, y, groundTruth] = generateData_dense(params,@dynModel);

y = y + params.magDisti;

%% Preliminaries on GP model
% nBasisFunctions = 256; % # basis functions
nBasisFunctions = 512; % # basis functions
groundTruth.nBasisFunctionsEst = nBasisFunctions;
d = 3; % The dimensionality of the inputs
groundTruth.d = d;
% Eigenbasis
LL = groundTruth.LL;
[eigenval,~,eigenfun_dx,NN] = domain_cartesian_dx(nBasisFunctions,d,LL);
lambda = eigenval(NN); % Eigenvalues
if max(NN(:,3))<2
    warning('Too few basis functions / misspecified domain size')
end
% Hyperparameters 
linSigma2   = theta(1); % Linear model scale parameter
lengthScale = theta(2); % Length scale parameter 
magnSigma2  = theta(3); % Magnitude parameter
sigma2 = theta(4);      % Noise variance parameter

% Evaluate the spectral density
Sse = @(w,lengthScale,magnSigma2) ...
    magnSigma2*sqrt(2*pi)^d*lengthScale^d*exp(-w.^2*lengthScale^2/2);   
S = @(w,linSigma2,lengthScale,magnSigma2) ...
    [linSigma2; linSigma2; linSigma2; Sse(w(1:end),lengthScale,magnSigma2)];
k = S(sqrt(lambda),linSigma2,lengthScale,magnSigma2);

%% Domain for visualization
x1t = groundTruth.fullMap_x1t;
x2t = groundTruth.fullMap_x2t;
[X1t,X2t] = meshgrid(x1t,x2t);
xt = [X1t(:) X2t(:) 0*X1t(:)];
nt = length(xt);
dPhixt = eigenfun_dx(NN,xt,1);
dPhiyt = eigenfun_dx(NN,xt,2);
dPhizt = eigenfun_dx(NN,xt,3);
dPhixt = [ones(nt,1), zeros(nt,2), dPhixt];
dPhiyt = [zeros(nt,1), ones(nt,1), zeros(nt,1), dPhiyt];
dPhizt = [zeros(nt,2), ones(nt,1), dPhizt];

%% Initial state and covariance, measurement noise covariance
% Initial linear and nonlinear state
x0_lin = [0; 0; 0; zeros(nBasisFunctions,1)];
x0_nonLin = initState;

% Initial covariance of the map
P0_lin = diag(k);

% Measurement noise covariance that was used to simulate data
R = sigma2 * eye(ny); 

%% Run filter
N_P = 100; % Number of particles
sparseFeatures = 0; % Run filter with dense features

runParticleFilter = 1;
if runParticleFilter
% Plotting settings
visualiseResults = 0;
% Run filter
[traj_max,traj_mean,xl_max,xl_mean,P_max,P_mean,traj_sample_iwmax,~] = ...
    particleFilter(@dynModel,@measModel,dx,y,...
    x0_nonLin,x0_lin,P0_lin,Q,R,N_P,dt,sparseFeatures,@makePlotsFilter_dense3D_magfield);

%% Compute RMSEs
% Compute position error for both maximum weight and weighted mean particles
% groundTruth.pos: Groundtruth position
% traj_max: Trajectory of maximum weight particles
% traj_mean: Trajectory of weighted mean particles
% traj_max_trans: transformed trajectory traj_max using the transformation
% from procrusted problem
% traj_mean_trans: transformed trajectory traj_mean using the transformation
% from procrusted problem
[~, traj_max_trans, ~] = procrustes(groundTruth.pos',traj_max(iPos,:)');
[~, traj_mean_trans, ~] = procrustes(groundTruth.pos',traj_mean(iPos,:)');

% Compute RMS position error represented in 3D
rmse_pos_max = rms(groundTruth.pos' - traj_max_trans);
rmse_pos_mean = rms(groundTruth.pos' - traj_mean_trans);

% Compute RMS orientation error
error_ori_max = zeros(length(traj_max),3);
error_ori_mean = zeros(length(traj_mean),3);
for ii = 1:length(traj_max)
    q_traj_max_trans = traj_max(iQuat,ii);
    q_traj_mean_trans = traj_mean(iQuat,ii);
    error_ori_max(ii,:) = ...
        quat2euler( qLeft(q_traj_max_trans) * ...
        qInv(groundTruth.quat(ii,:)'));
    error_ori_mean(ii,:) = ...
        quat2euler(qLeft( q_traj_mean_trans ) * ...
        qInv(groundTruth.quat(ii,:)'));
end
rmse_ori_max = rms(error_ori_max);
rmse_ori_mean = rms(error_ori_mean);

rmseFilter = [rmse_pos_max,rmse_pos_mean,...
  rmse_ori_max,rmse_ori_mean];

outputstr = '%.4f %.4f %.4f : ';
outputstr = repmat(outputstr,1,4);
fprintf(['RMSE filter : ' outputstr ' \n'], rmseFilter)

if nargout >= 1
    savedData = {};
    savedData.traj_max = traj_max;
    savedData.traj_mean = traj_mean;
    savedData.traj_sample_iwmax = traj_sample_iwmax;
    savedData.xl_max = xl_max;
    savedData.P_max = P_max;
    savedData.xl_mean = xl_mean;
    savedData.P_mean = P_mean;
    savedData.rmseFilter = rmseFilter;
end

end

%% Run smoother
runSmoother = 1;
if runSmoother
dynResNorm = @(xnk,xni,dx,dt,Q) [xnk(iPos) - xni(iPos) - dx(1:3)' ; ...
    logq( qLeft( qLeft(qInv(dx(4:7))) * qInv(xni(iQuat)) ) * xnk(iQuat) )]'/chol(dt * Q,'lower');
                
if N_K ~= 0
    if makePlots
        [XNK,XLK,PK] = particleSmoother(@dynModel,@measModel,dynResNorm,dx,y,...
            x0_nonLin,x0_lin,P0_lin,Q,R,N_P,N_K,dt,sparseFeatures,@makePlotsSmoother_dense3D_magfield);
    else
        [XNK,XLK,PK] = particleSmoother(@dynModel,@measModel,dynResNorm,dx,y,...
            x0_nonLin,x0_lin,P0_lin,Q,R,N_P,N_K,dt,sparseFeatures,[]);
    end
end

%% Compute RMSEs
rmseSmoother = zeros(N_K,6);
for k = 1:N_K
    [~, traj_smooth_trans, ~] = procrustes(groundTruth.pos',XNK(iPos,:,k)');
    rmse_pos = rms(groundTruth.pos' - traj_smooth_trans);
    
    % Compute RMS orientation error
    N_T = size(XNK,2);
    error_ori = zeros(N_T,3);
    for ii = 1:N_T
        q_traj_trans = XNK(iQuat,ii,k);
        error_ori(ii,:) = ...
            quat2euler( qLeft(q_traj_trans) * ...
            qInv(groundTruth.quat(ii,:)'));
    end
    rmse_ori = rms(error_ori);
    rmseSmoother(k,:) = [rmse_pos, rmse_ori];

    % Report
    outputstr = '%.4f %.4f %.4f : ';
    outputstr = repmat(outputstr,1,2);
    fprintf(['RMSE smoother k=%02i : ' outputstr '\n'], k, rmseSmoother(k,:))
end

if nargout >= 1 && N_K~=0
    savedData.XNK = XNK;
    savedData.XLK = XLK;
    savedData.PK = PK;
    savedData.rmseSmoother = rmseSmoother;
end
end

%% Run EKF
x0_ekf = [x0_nonLin(iPos) ; zeros(3,1) ; x0_lin];
q0_ekf = x0_nonLin(iQuat);
P0_ekf = blkdiag(zeros(6),P0_lin);
[xf_ekf,qnb_ekf,Pf_ekf] = ekf_dense(@dynModel_ekf,@measModel_ekf,dx,y,...
            x0_ekf,q0_ekf,P0_ekf,Q,R,dt);

[~, traj_ekf_trans, ~] = procrustes(groundTruth.pos',xf_ekf(1:3,:)');
rmse_pos_ekf = rms(groundTruth.pos' - traj_ekf_trans);

if nargout >= 1
    savedData.xf_ekf = xf_ekf;
    savedData.Pf_ekf = Pf_ekf;
    savedData.qnb_ekf = qnb_ekf;
    savedData.rmse_pos_ekf = rmse_pos_ekf;
end

%% Supporting functions
function dy = measModel(xn)
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
end

function [yhat,dy] = measModel_ekf(x,q)
    dPhix = eigenfun_dx(NN,x(iPos)',1);
    dPhiy = eigenfun_dx(NN,x(iPos)',2);
    dPhiz = eigenfun_dx(NN,x(iPos)',3);
    dPhix = [1, 0, 0, dPhix];
    dPhiy = [0, 1, 0, dPhiy];
    dPhiz = [0, 0, 1, dPhiz];
    Rnb = quat2rmat(q);
    dPhi = [dPhix ; dPhiy ; dPhiz];
    yhat = Rnb' * dPhi * x(7:end);
    dy = zeros(ny,size(dPhix,2)+6);
    J=reshape(...
        JacobianPhi3D(x(iPos),nBasisFunctions,LL(1,1),LL(2,1),LL(1,2),LL(2,2),LL(1,3),LL(2,3),NN),...
        9,nBasisFunctions)*x(10:end);
    J = reshape(J,3,3);
    dy(:,1:3) = Rnb' * J;
    dy(:,4:6) = Rnb' * mcross(dPhi * x(7:end));
    dy(:,7:end) = Rnb' * dPhi;
end

function [xpred, dQuat] = dynModel(xn,dx,dt,Q)
    % Predict through dynamic model. Also optionally output dQuat for
    % generating odometry data
    xpred_pos = xn(iPos) + dx(iPos)' + chol(dt * Q(iPos,iPos),'lower')*randn(3,1);
    dQuat = qLeft(dx(iQuat)') * expq(chol(dt * Q(4:6,4:6),'lower')*randn(3,1));
    xpred_quat = qLeft(xn(iQuat)) * dQuat;
    xpred = [xpred_pos ; xpred_quat]; 
end

function [xpred, qpred, F, G] = dynModel_ekf(x,q,dx)
    xpred = x;
    xpred(iPos) = x(iPos) + dx(iPos)';
    qpred = qLeft(q) * dx(iQuat)';
    F = eye(size(x,1));
    G = [blkdiag(eye(3),quat2rmat(qpred)) ; zeros(size(x,1)-6,6)];
end

function makePlotsFilter_dense3D_magfield(...
        xn,xl_max,~,traj_max,yhattraj,~,~,~,~)
    % Predicted GP for visualisation
    dEft = reshape([dPhixt*xl_max;...
             dPhiyt*xl_max; ...
             dPhizt*xl_max],[nt 3]);
    
    t = sum(~isnan(traj_max(1,:)));
    % Plot
    if visualiseResults
        figure(3);
        subplot(121); cla; hold on
          imagesc(x1t,x2t,reshape(sqrt(sum(dEft.^2,2)),size(X1t)));
          colorbar
          caxis([0 max(sqrt(sum(y.^2,2)))+10])
          axis equal
          xlim([min(xt(:,1)) max(xt(:,1))])
          ylim([min(xt(:,2)) max(xt(:,2))])
          title('Estimated magnetic field map')
        subplot(122); cla; hold on  
            caxis([0 max(sqrt(sum(y.^2,2)))+10])
           scatter(traj_max(1,1:t),traj_max(2,1:t),100,sqrt(yhattraj(1,1:t).^2 + ...
               yhattraj(2,1:t).^2 + yhattraj(3,1:t).^2)','filled')
           plot(traj_max(1,1:t),traj_max(2,1:t),'k')
           scatter(xn(1,:),xn(2,:),'k')
           axis equal
           xlim([min(xt(:,1)) max(xt(:,1))])
           ylim([min(xt(:,2)) max(xt(:,2))])
    end
end


function makePlotsSmoother_dense3D_magfield(xnk,xlk,k,XNK,XLK,PK)
    % Predict magnetic field map on sampled trajectory ...
    N_T = size(xnk,2);
    dPhix = eigenfun_dx(NN,xnk(iPos,:)',1);
    dPhiy = eigenfun_dx(NN,xnk(iPos,:)',2);
    dPhiz = eigenfun_dx(NN,xnk(iPos,:)',3);
    dPhix = [ones(N_T,1), zeros(N_T,2), dPhix];
    dPhiy = [zeros(N_T,1), ones(N_T,1), zeros(N_T,1), dPhiy];
    dPhiz = [zeros(N_T,2), ones(N_T,1), dPhiz];
    dEft_pred = reshape([dPhix*xlk;...
             dPhiy*xlk; ...
             dPhiz*xlk],[N_T 3]);
    
    % ... and on whole grid of testing points
    dEft = reshape([dPhixt*xlk;...
             dPhiyt*xlk; ...
             dPhizt*xlk],[nt 3]);

    figure(4);
    subplot(121); cla; hold on
      imagesc(x1t,x2t,reshape(sqrt(sum(dEft.^2,2)),size(X1t)));
      colorbar
      caxis([0 max(sqrt(sum(y.^2,2)))+10])
      axis equal
      xlim([min(xt(:,1)) max(xt(:,1))])
      ylim([min(xt(:,2)) max(xt(:,2))])
      title('Estimated magnetic field map')
    subplot(122); cla; hold on  
       caxis([0 max(sqrt(sum(y.^2,2)))+10])
       scatter(xnk(1,:),xnk(2,:),100,sqrt(dEft_pred(:,1).^2 + ...
           dEft_pred(:,2).^2 + dEft_pred(:,3).^2)','filled')
       plot(xnk(1,:),xnk(2,:),'k')
       
       axis equal
       xlim([min(xt(:,1)) max(xt(:,1))])
       ylim([min(xt(:,2)) max(xt(:,2))])
       
       if k <= 9
           figure(5); subplot(3,3,k); cla; hold on
           title(['It. ' num2str(k) ])
           caxis([0 max(sqrt(sum(y.^2,2)))+10])
           scatter(xnk(1,:),xnk(2,:),100,sqrt(dEft_pred(:,1).^2 + ...
               dEft_pred(:,2).^2 + dEft_pred(:,3).^2)','filled')
           plot(xnk(1,:),xnk(2,:),'k')

           axis equal
           xlim([min(xt(:,1)) max(xt(:,1))])
           ylim([min(xt(:,2)) max(xt(:,2))])
       elseif k > 9 && k <= 18
           figure(6); subplot(3,3,k-9); cla; hold on
           title(['It. ' num2str(k) ])
           caxis([0 max(sqrt(sum(y.^2,2)))+10])
           scatter(xnk(1,:),xnk(2,:),100,sqrt(dEft_pred(:,1).^2 + ...
               dEft_pred(:,2).^2 + dEft_pred(:,3).^2)','filled')
           plot(xnk(1,:),xnk(2,:),'k')

           axis equal
           xlim([min(xt(:,1)) max(xt(:,1))])
           ylim([min(xt(:,2)) max(xt(:,2))])
       elseif k > 18 && k <= 27
           figure(7); subplot(3,3,k-18); cla; hold on
           title(['It. ' num2str(k) ])
           caxis([0 max(sqrt(sum(y.^2,2)))+10])
           scatter(xnk(1,:),xnk(2,:),100,sqrt(dEft_pred(:,1).^2 + ...
               dEft_pred(:,2).^2 + dEft_pred(:,3).^2)','filled')
           plot(xnk(1,:),xnk(2,:),'k')

           axis equal
           xlim([min(xt(:,1)) max(xt(:,1))])
           ylim([min(xt(:,2)) max(xt(:,2))])
       end
    end


end