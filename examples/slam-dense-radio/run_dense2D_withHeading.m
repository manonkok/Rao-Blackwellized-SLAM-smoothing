function [rmseFilter,rmseSmoother,savedData,groundTruthSeq] = run_dense2D_withHeading(params)
%% RUN_DENSE2D_WITHHEADING - Runs a dense 2D simulation with heading information
%
% Syntax:
%   [rmseFilter, rmseSmoother, savedData, groundTruthSeq] = run_dense2D_withHeading(params)
%
% In:
%   params - Struct containing various parameters for the simulation
%
% Out:
%   rmseFilter      - Root Mean Square Error (RMSE) of the filter
%   rmseSmoother    - RMSE of the smoother
%   savedData       - Data saved from the simulation
%   groundTruthSeq  - Sequence of ground truth data for each Monte Carlo iteration
%
% Description:
%   This function runs a dense 2D simulation with heading information. It 
%   generates data based on dynamic models and process noise covariances for 
%   different trajectory types such as 'line_3D' and 'square_3D'. The 
%   function runs both a filter and smoother, calculating the RMSE for each. 
%   
% See also:
%   generateData_dense
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

% Parameters should either be input via params struct or standard
% parameters are used 
if nargin < 1
    % Default settings
    params = [];
    % Trajectory type
    trajType = 'line_3D';
    params.trajType = trajType;
    
    % GP hyperparameters
    params.theta = [0.25 ; 2 ; 0.01]; 
    theta = params.theta;
    
    % Number of smoother iterations
    N_K = 20;
    
    % If file run separately then we want to always use the same seed and
    % make plots
    makePlots = 1;
    params.makePlots = makePlots;
    sameRndSeed = 1;
    nMC = 1;
else
    % If parameters input via params struct ... 
    % ... don't use the same seed and don't visualise any intermediate
    % results
    theta = params.theta; 
    N_K = params.N_K;
    trajType = params.trajType;
    % If we run this code in a loop we don't always want to get the same results
    sameRndSeed = 0;
    makePlots = params.makePlots;
    nMC = params.nMC;
end

%% Make dynamic models and process noise covariances
switch trajType
    case 'line_3D'
        T = 1;
        params.T = T;
        params.dt = T;
        N = 32;
        params.N = N;
        Q = 1E-6 * ones(N,1);
        Q(N/2) = 0.3^2;
        Q = reshape(Q,[1 1 N]);
        params.Q = Q;
        dynModel = @(xn,dx,dt,Q) [xn(1:2) + [cos(xn(3)), -sin(xn(3)) ; ...
                sin(xn(3)), cos(xn(3))]' * dx(1:2)' ; xn(3) + dx(3) + chol(dt * Q,'lower') * randn];  
        dynResNorm = @(xnk,xni,dx,dt,Q) (xnk(3) - xni(3) - dx(3))'/chol(dt * Q,'lower');
 
    case 'square_3D'
        T = 1;
        params.T = T;
        params.dt = T;
        N = 48;
        params.N = N;
        Q = 1E-6 * ones(N,1);
        Q(N/4 + N/4*(0:2)) = 0.1^2;        
        Q = reshape(Q,[1 1 N]);
        params.Q = Q;
        dynModel = @(xn,dx,dt,Q) [xn(1:2) + [cos(xn(3)), -sin(xn(3)) ; ...
            sin(xn(3)), cos(xn(3))]' * dx(1:2)' ; xn(3) + dx(3) + chol(dt * Q,'lower') * randn];
        dynResNorm = @(xnk,xni,dx,dt,Q) (xnk(3) - xni(3) - dx(3))'/chol(dt * Q,'lower');    
end

%% Generate data and define details about model
if sameRndSeed
  rng(1,'twister') 
end

visualiseResults = 1;
params.visualiseResults = visualiseResults;

[dx, initState, y, groundTruth] = generateData_dense(params,dynModel);

% Define indices
iPos = 1:2;
 
%% Preliminaries on GP model
nBasisFunctions = 128;              % Linear states: basis functions
groundTruth.nBasisFunctionsEst = nBasisFunctions;
d = 2; % The dimensionality of the inputs
% Domain determined in pre-processing: x \in [-L,L]
LL = groundTruth.LL;
% Eigenbasis
[eigenval,eigenfun,~,NN] = domain_cartesian_dx(nBasisFunctions,d,LL);
lambda = eigenval(NN); % Eigenvalues
if max(NN(:,1))<2 || max(NN(:,2))<2 % Check if there are at least 2 basis functions in each dimension
    warning('Too few basis functions / misspecified domain size')
end

% Hyperparameters 
lengthScale = theta(1); % Length scale parameter 
magnSigma2  = theta(2); % Magnitude parameter
sigma2 = theta(3);      % Noise variance parameter

% Evaluate the spectral density
S = @(w,lengthScale,magnSigma2) ...
    magnSigma2*sqrt(2*pi)^d*lengthScale^d*exp(-w.^2*lengthScale^2/2);   
k = S(sqrt(lambda),lengthScale,magnSigma2);

%% Domain for visualization
x1t = groundTruth.fullMap_x1t;
x2t = groundTruth.fullMap_x2t;
[X1t,X2t] = meshgrid(x1t,x2t);
xt = [X1t(:) X2t(:)];
Phit = eigenfun(NN,xt);

%% Initial state and covariance, measurement noise
% Initial linear and nonlinear states
x0_nonLin = initState;
x0_lin = zeros(nBasisFunctions,1);

% Initial covariance of the map
P0_lin = diag(k);

% Measurement noise and process noise
R = sigma2 * eye(size(y,2));

%% Run Monte Carlo simulations for filter and smoother using same map but 
% different realisation of the measurement noise and the odometry noise
savedData = cell(1,nMC);
groundTruthSeq = cell(nMC,1);
groundTruthSeq{1} = groundTruth;
rmseFilter = zeros(nMC,3);
rmseSmoother = zeros(N_K,2,nMC);
for iMC = 1:nMC
    % Make new odometry
    if iMC ~= 1
        fieldData.f = groundTruth.f;
        fieldData.y = y;
        [dx, ~, y, groundTruth] = generateData_dense(params,dynModel,fieldData);
        groundTruthSeq{iMC} = groundTruth;
    end
    
    % Run filter
    N_P = 100; % Number of particles
    sparseFeatures = 0;

    measModel = @(xn) eigenfun(NN,xn(iPos,:)');
    if makePlots
        plotFunctionFilter = @makePlotsFilter_dense2D_withHeading;
    else
        plotFunctionFilter = [];
    end

    % Run filter
    [traj_max,traj_mean,xl_max,xl_mean,P_max,P_mean,~,xn_traj] = particleFilter(dynModel,measModel,dx,y,...
            x0_nonLin,x0_lin,P0_lin,Q,R,N_P,T,sparseFeatures,plotFunctionFilter);

    % Fill output struct
    if nargout >= 3
        savedData{iMC}.traj_max = traj_max;
        savedData{iMC}.traj_mean = traj_mean;
        savedData{iMC}.xl_max = xl_max;
        savedData{iMC}.P_max = P_max;
        savedData{iMC}.xl_mean = xl_mean;
        savedData{iMC}.P_mean = P_mean;
        if iMC == nMC
            savedData{iMC}.xn_traj = xn_traj;
        end
    end

    %% Run smoother
    if N_K ~= 0
        if makePlots
            plotFunctionSmoother = @makePlotsSmoother_dense2D_withHeading;
        else
            plotFunctionSmoother = [];
        end

        [XNK,XLK,PK] = particleSmoother(dynModel,measModel,dynResNorm,dx,y,...
                x0_nonLin,x0_lin,P0_lin,Q,R,N_P,N_K,T,sparseFeatures,plotFunctionSmoother);
    end
    
    if nargout >= 3 && N_K~=0
        savedData{iMC}.XNK = XNK;
        savedData{iMC}.XLK = XLK;
        savedData{iMC}.PK = PK;
    end
end

%% Plotting function
function makePlotsFilter_dense2D_withHeading(xn,xl_max,~,traj_max,yhattraj,~,~,~,~)
    
    % Predicted GP for visualisation
    Eft = Phit*xl_max;
    
    t = sum(~isnan(traj_max(1,:)));
    % Plot
    figure(3);
    subplot(121); cla; hold on
      imagesc(x1t,x2t,reshape(Eft,size(X1t)));
      colorbar
      caxis([min(groundTruth.f) max(groundTruth.f)])
      axis equal
      xlim([min(xt(:,1)) max(xt(:,1))])
      ylim([min(xt(:,2)) max(xt(:,2))])
      title('Estimated map')
    subplot(122); cla; hold on  
       caxis([min(groundTruth.f) max(groundTruth.f)])
       scatter(traj_max(1,1:t),traj_max(2,1:t),100,yhattraj(1:t)','filled')
       plot(traj_max(1,1:t),traj_max(2,1:t),'k')
       scatter(xn(1,:),xn(2,:),'k')
       axis equal
       xlim([min(xt(:,1)) max(xt(:,1))])
       ylim([min(xt(:,2)) max(xt(:,2))])
end

%%
    function makePlotsSmoother_dense2D_withHeading(xnk,xlk,k,XNK,XLK,PK)
    % Predict magnetic field map on sampled trajectory ...
    Phi = eigenfun(NN,xnk(1:2,:)');
    Eft_pred_traj = Phi*xlk;
    
    % ... and on whole grid of testing points
    Eft = Phit*xlk; 
    if params.visualiseResults
        figure(3);
            subplot(121); cla; hold on
                imagesc(x1t,x2t,reshape(Eft,size(X1t)));
                colorbar
                caxis([min(groundTruth.f) max(groundTruth.f)])
                axis equal
                xlim([min(xt(:,1)) max(xt(:,1))])
                ylim([min(xt(:,2)) max(xt(:,2))])
                title('Estimated magnetic field map')
            subplot(122); cla; hold on  
                caxis([min(groundTruth.f) max(groundTruth.f)])
                scatter(xnk(1,:),xnk(2,:),100,Eft_pred_traj,'filled')
                plot(xnk(1,:),xnk(2,:),'k')
                axis equal
                xlim([min(xt(:,1)) max(xt(:,1))])
                ylim([min(xt(:,2)) max(xt(:,2))])

        if k <= 9
            figure(4); subplot(3,3,k); cla; hold on
            title(['It. ' num2str(k) ])
            caxis([min(groundTruth.f) max(groundTruth.f)])
            scatter(xnk(1,:),xnk(2,:),100,Eft_pred_traj,'filled')
            plot(xnk(1,:),xnk(2,:),'k')
            axis equal
            xlim([min(xt(:,1))-2 max(xt(:,1))+2])
            ylim([min(xt(:,2)) max(xt(:,2))])
        elseif k > 9 && k <= 18
            figure(5); subplot(3,3,k-9); cla; hold on
            title(['It. ' num2str(k) ])
            caxis([min(groundTruth.f) max(groundTruth.f)])
            scatter(xnk(1,:),xnk(2,:),100,Eft_pred_traj,'filled')
            plot(xnk(1,:),xnk(2,:),'k')
            axis equal
            xlim([min(xt(:,1))-2 max(xt(:,1))+2])
            ylim([min(xt(:,2)) max(xt(:,2))])
        elseif k > 18 && k <= 27
            figure(6); subplot(3,3,k-18); cla; hold on
            title(['It. ' num2str(k)])
            caxis([min(groundTruth.f) max(groundTruth.f)])
            scatter(xnk(1,:),xnk(2,:),100,Eft_pred_traj,'filled')
            plot(xnk(1,:),xnk(2,:),'k')
            axis equal
            xlim([min(xt(:,1))-2 max(xt(:,1))+2])
            ylim([min(xt(:,2)) max(xt(:,2))])
        end
    end
end

end