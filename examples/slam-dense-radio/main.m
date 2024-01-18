%% DENSE (~RADIO) SLAM EXAMPLE
%
% This example demonstrates the use of the filtering vs. smoothing
% approaches presented in the paper (see Sec. 5.1). The model estimates a
% Gaussian process map for a dense signal field with one output
% dimensionality. This can correspond, e.g., to dense radio SLAM where the
% map/field is the RSSI signal strength.
% 
% The switches below are for running and plotting the results in the
% Figs. 3 and 4 in the paper. 
%
% Copyright:
%   2023-   Manon Kok and Arno Solin


%% Make simulation results for radio SLAM

  clear, clc, close all
  rng(1,'twister') % To make sure we always get the same results when re-running
  addpath(genpath('../../.'))

  % Default settings
  params = [];
  params.theta = [0.25 ; 2 ; 0.01];
  params.trajType = 'line_3D';
  params.N_K = 50;
  params.nMC = 100; % Change this to 1 if you just quickly want to test
  params.makePlots = 0;

  % Settings for running (choose here what to run)
  runAlgs = false;
  runAlgsDegeneracy = false;
  makePlotsLine = true;
  makePlotsDegeneracy = true;
  savePlots = false;

%% Make plots 2D line

  if runAlgs
    params.nLL = 4;
    [rmseFilter,rmseSmoother,savedData,groundTruth] = ...
        run_dense2D_withHeading(params);
    save('results/resultsSimRadioSLAM','rmseFilter','rmseSmoother','savedData','groundTruth')
  end
  if makePlotsLine
    load('results/resultsSimRadioSLAM')

    nMC = length(savedData);

    % Plotting options
    color = [0, 93, 141]/255;
    cmap = parula(64);
    clim = [-1 1]*max(abs([min(groundTruth{1}.fullMap_f) max(groundTruth{1}.fullMap_f)]));
    alim = []; 
    lims = [-.7 .7 -2 2]; % xlim / ylim
    
    % Make map based on true values f, which were constant over the
    % simulations
    m = groundTruth{1}.nBasisFunctionsEst; % Basis functions
    d = 2;
    LL = groundTruth{1}.LL; % Domain
    x1t = groundTruth{1}.fullMap_x1t;
    x2t = groundTruth{1}.fullMap_x2t;
    [X1t,X2t] = meshgrid(x1t,x2t);
    xt = [X1t(:) X2t(:)];
    [~,eigenfun,~,NN] = domain_cartesian_dx(m,d,LL);
    Phit = eigenfun(NN,xt);
    
    % Make plot odometry + true map
    figure(1), clf, hold all
      %title('Odometry')
      set(gcf,'color','w')
      imagesc(x1t,x2t,reshape(groundTruth{1}.fullMap_f,size(X1t)));

      % Trajectories
      for i = 1:nMC
        plot(groundTruth{i}.odometry(:,1),groundTruth{i}.odometry(:,2),'-', ...
            'color',color,'linewidth',1)
      end
      box on
      %colorbar
      caxis(clim)
      axis equal on
      xlim(lims(1:2))
      ylim(lims(3:4))
      set(gca,'YTickLabel',[]);
      set(gca,'XTickLabel',[]);
    
    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'line-odometry.png') 
    end      
      
    % Make plot maximum weight particle
    figure(2), clf, hold all
      set(gcf,'color','w')
      % Map from one of the maximum weight particles
      Eft = Phit*savedData{1}.xl_max;
      Varft = sum((Phit*chol(savedData{1}.P_max,'lower')).^2,2);    
      imagescalpha(x1t,x2t,reshape(Eft,size(X1t)),reshape(sqrt(Varft),size(X1t)),alim);
      
      % Trajectories
      for i = 1:nMC
        plot(savedData{i}.traj_max(1,:),savedData{i}.traj_max(2,:),'-', ...
            'color',color,'linewidth',1)
      end
      box on
      %colorbar
      axis equal on
      xlim(lims(1:2))
      ylim(lims(3:4))
      set(gca,'YTickLabel',[]);
      set(gca,'XTickLabel',[]);
      
    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'line-filter-max.png') 
    end
    
    % Make plot weighted mean particles
    figure(3), clf, hold all
      set(gcf,'color','w')
      % Map from one of the weighted mean particles
      Eft = Phit*savedData{1}.xl_mean;
      Varft = sum((Phit*chol(savedData{1}.P_mean,'lower')).^2,2);    
      imagescalpha(x1t,x2t,reshape(Eft,size(X1t)),reshape(sqrt(Varft),size(X1t)),alim);
      
      % Trajectories
      for i = 1:nMC
        plot(savedData{i}.traj_mean(1,:),savedData{i}.traj_mean(2,:),'-', ...
            'color',color,'linewidth',1)
      end
      box on
      %colorbar
      axis equal on
      xlim(lims(1:2))
      ylim(lims(3:4))
      set(gca,'YTickLabel',[]);
      set(gca,'XTickLabel',[]);

    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'line-filter-mean.png') 
    end

    % Make plot sampled trajectories PS
    figure(4), clf, hold all
      set(gcf,'color','w')
      %title('Samples smoother, after 10 iterations')
      % Map from one of the trajectories from one of the iterations
      Eft = Phit*savedData{1}.XLK(:,end);
      Varft = sum((Phit*chol(savedData{1}.PK(:,:,end),'lower')).^2,2);    
      imagescalpha(x1t,x2t,reshape(Eft,size(X1t)),reshape(sqrt(Varft),size(X1t)),alim);
          
      % Trajectories
      N_K = size(savedData{1}.XNK,3);
      for i = 1:nMC
        for k = 1:N_K
            plot(squeeze(savedData{i}.XNK(1,:,k)),...
                squeeze(savedData{i}.XNK(2,:,k)),'-', ...
                'color',color,'linewidth',1)
        end
      end
      
      box on
      %colorbar
      axis equal on
      xlim(lims(1:2))
      ylim(lims(3:4))
      set(gca,'YTickLabel',[]);
      set(gca,'XTickLabel',[]);

    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'line-smoother.png') 
    end
  end

%% Make plots path degeneracy

  if runAlgsDegeneracy
    rng(2,'twister') % To make sure we always get the same results when re-running

    % Default settings
    params = [];
    params.theta = [0.25 ; 2 ; 0.01]; 
    params.nLL = 4;
    params.trajType = 'square_3D';
    params.N_K = 50;
    params.nMC = 1;
    params.makePlots = 1;
    params.makePlotsDegeneracy = 1;
    params.visualiseResults = 0;

    [rmseFilter,rmseSmoother,savedData,groundTruth] = ...
        run_dense2D_withHeading(params);

    save('results/resultsSimRadioSLAM_degen','rmseFilter','rmseSmoother','savedData','groundTruth')
  end
  if makePlotsDegeneracy
    load('results/resultsSimRadioSLAM_degen')
    
    % Plotting options
    color = [0, 93, 141]/255;
    cmap = parula(64);
    clim = [-1 1]*max(abs([min(groundTruth{1}.fullMap_f) max(groundTruth{1}.fullMap_f)]));
    alim = [0.05 1.39]; % Alpha limits
    m = groundTruth{1}.nBasisFunctionsEst; % Basis functions
    d = 2;
    pos = groundTruth{1}.pos; % Since ground truth position always the same
    LL = groundTruth{1}.LL; % Domain
    x1t = groundTruth{1}.fullMap_x1t;
    x2t = groundTruth{1}.fullMap_x2t;
    [X1t,X2t] = meshgrid(x1t,x2t);
    xt = [X1t(:) X2t(:)];
    [eigenval,eigenfun,~,NN] = domain_cartesian_dx(m,d,LL);
    Phit = eigenfun(NN,xt);
    lims = [[min(xt(:,1))+.3 max(xt(:,1))-.3] ...
            [min(xt(:,2))+.3 max(xt(:,2))-.3]];
    
    % Plot ground-truth
    figure(6); clf; hold on
        set(gcf,'color','w')    
        imagesc(x1t,x2t,reshape(groundTruth{1}.fullMap_f,size(X1t)));
        
        axis equal square on
        xlim(lims(1:2))
        ylim(lims(3:4))
        box on
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        plot([-1 -1 1 1 -1],[-1 1 1 -1 -1],'-','color',color,'linewidth',1)
        set(gca,'Layer','top')
        colorbar
            
    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'degeneracy-gt.png') 
    end
  
    % Plot degenerate filtering solution
    figure(7), clf, hold all
        set(gcf,'color','w')    
        % Compute map of highest weight particle, which is exactly the one
        % of which we plot the trajectory
        Eft = Phit*savedData{1}.xl_max;
        Varft = sum((Phit*chol(savedData{1}.P_max,'lower')).^2,2);
        
        % Visualize magnetic field with uncertainties
        imagescalpha(x1t,x2t,reshape(Eft,size(X1t)),reshape(sqrt(Varft),size(X1t)),alim);
             
        % Trajectories
        for iParticle = 1:size(savedData{end}.xn_traj,2)
            plot(squeeze(savedData{end}.xn_traj(1,iParticle,:))',...
                squeeze(savedData{end}.xn_traj(2,iParticle,:))','-', ...
                'color',color,'linewidth',1)
        end
                
        axis equal square on
        xlim(lims(1:2))
        ylim(lims(3:4))
        box on
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(gca,'Layer','top')        
    
    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'degeneracy-filter.png') 
    end
    
    % Plot non-degenerate smoother solution  
    figure(8), clf, hold all
        set(gcf,'color','w')
        % Compute map of one of the trajectories
        Eft = Phit*savedData{1}.XLK(:,end);
        Varft = sum((Phit*chol(savedData{1}.PK(:,:,end),'lower')).^2,2);
        h=imagescalpha(x1t,x2t,reshape(Eft,size(X1t)),reshape(sqrt(Varft),size(X1t)),alim);
        
        for iSampleTraj = 1:size(savedData{1}.XNK,3)
            plot(squeeze(savedData{1}.XNK(1,:,iSampleTraj))',...
                squeeze(savedData{1}.XNK(2,:,iSampleTraj))','-', ...
                'color',color,'linewidth',1)
        end
        axis xy equal square on
        xlim(lims(1:2))
        ylim(lims(3:4))
        box on
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(gca,'Layer','top')
        colorbar
        
    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'degeneracy-smoother.png') 
    end
  end
