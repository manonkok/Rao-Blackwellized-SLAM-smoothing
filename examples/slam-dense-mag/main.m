%% DENSE MAGNETIC SLAM EXAMPLE
%
% This example demonstrates the use of the filtering vs. smoothing
% approaches presented in the paper (see Sec. 5.2). The model estimates a
% Gaussian process map for a dense magnetic field with three output
% dimension.
% 
% The switches below are for running and plotting the results in Fig. 5
% in the paper. 
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

%% Make plots for dense toy example

  clear, clc, close all
  rng(1,'twister') % To make things deterministic so we always get the same results
  addpath(genpath('../../.'))

  % Default settings
  params = [];
  params.Q = blkdiag(diag(10^2*[0.05^2,0.05^2,0.01^2]), diag([0.01 0.01 0.3]*pi/180).^2);
  params.theta = [650 ; 1.2 ; 200 ; 10];

  params.trajType = 'bean_6D';
  N_K = 10;
  params.N_K = N_K;
  params.makePlots = 0; % For plotting during filter and smoother

  % Settings for running (choose what to run)
  runComparisons = false;
  makePlots = true; % For plotting in this script
  savePlots = false;

%% Run code to generate results for different disturbances
if runComparisons
    nSim = 20;
    rmses_ekf_3D = zeros(nSim,3);
    rmses_pf_3D = zeros(nSim,12);
    rmses_ps_3D = zeros(nSim,N_K,3);
    magDist = [zeros(1,3) ; 0, 1, 0 ; 0, 5, 0 ; 0 , 10, 0];
    for iMagDist = 1:size(magDist,1)
        params.magDisti = magDist(iMagDist,:);
        for i = 1:nSim
            [savedData{i},groundTruth{i}] = ...
                run_dense3D_magfield(params);
            rmses_ekf_3D(i,:) = savedData{i}.rmse_pos_ekf;
            rmses_pf_3D(i,:) = savedData{i}.rmseFilter;
            rmses_ps_3D(i,:,:) = savedData{i}.rmseSmoother(:,1:3);
        end
        totalRMSE = [sqrt(mean(mean(rmses_ekf_3D.^2))),...
            sqrt(mean(mean(rmses_pf_3D(:,1:3).^2))),...
            sqrt(mean(mean(rmses_pf_3D(:,4:6).^2))), ...
            sqrt(mean(mean(mean(rmses_ps_3D(:,:,:).^2))))];
        rmses_ekf_pf = [sqrt(mean(rmses_ekf_3D.^2,2)),sqrt(mean(rmses_pf_3D(:,1:3).^2,2)),...
            sqrt(mean(rmses_pf_3D(:,4:6).^2,2))];
        rmses_ps = [sqrt(mean(rmses_ps_3D.^2,3))];
        filename = ['results/resultsSimMagSLAM_dist' num2str(magDist(iMagDist,2))];
        save(filename,'savedData','groundTruth','rmses_ekf_pf','rmses_ps','magDist')
    end
end

%% Make plots
if makePlots
    magDist = [zeros(1,3) ; 0, 1, 0 ; 0, 5, 0 ; 0 , 10, 0];
    % Load results 
    dist0 = load('results/resultsSimMagSLAM_dist0.mat');
    dist1 = load('results/resultsSimMagSLAM_dist1.mat');
    dist5 = load('results/resultsSimMagSLAM_dist5.mat');
    dist10 = load('results/resultsSimMagSLAM_dist10.mat');

    % Extract data for plotting; extract indices 1 and 3 since only
    % weighted mean for PF plotted
    X = [dist0.rmses_ekf_pf(:,[1,3]),dist0.rmses_ps(:,end), ...
                dist1.rmses_ekf_pf(:,[1,3]),dist1.rmses_ps(:,end), ...
                dist5.rmses_ekf_pf(:,[1,3]),dist5.rmses_ps(:,end), ...
                dist10.rmses_ekf_pf(:,[1,3]),dist10.rmses_ps(:,end)];

    figure(1);clf; hold on
        lims = [0 0.3];
        labels = {};
        
        for i = 1:size(magDist,1)
          labels = [labels {'EKF','PF','PS'}]; %#ok
        end
        
        h = fill([.5 3.5 3.5 .5 .5],[lims(1) lims(1) lims(2) lims(2) lims(1)],1);
        set(h,'FaceColor',[.9 .9 .9],'EdgeColor',[.95 .95 .95])
        
        h = fill([.5 3.5 3.5 .5 .5]+6,[lims(1) lims(1) lims(2) lims(2) lims(1)],1);
        set(h,'FaceColor',[.9 .9 .9],'EdgeColor',[.95 .95 .95])
        
        h=boxplot(X,'colors','rrb','boxstyle','outline');
        
        mycolor = [0.00000,0.36471,0.55294];
        
        % Customize boxplot
        for i=1:numel(h)
          if strcmpi(get(h(i),'LineStyle'),'--')
            set(h(i),'LineStyle','-')
          end
          if strcmpi(get(h(i),'tag'),'Box')
            foo = get(h(i),'XData');
            mid = (max(foo)+min(foo))/2;
            set(h(i),'XData',mid+.5*(foo-mid))
          end
          if strcmpi(get(h(i),'tag'),'Box') && all(get(h(i),'Color')==[0 0 1])
            hp=patch(get(h(i),'XData'),get(h(i),'YData'),1);
            set(hp,'FaceColor',mycolor)
            set(hp,'EdgeColor',mycolor)
          end
          set(h(i),'Color',mycolor)
        end
        
        % Set labels
        for i=1:size(magDist,1)
          text(2+(i-1)*3,lims(end)*0.95, ...
              sprintf('$o = %.0f$',magDist(i,2)), ...
              'HorizontalAlignment','center')
        end
        
        ylabel('RMSE in final SLAM path')
        ylim(lims)
        box on
        xlim([.5 size(X,2)+.5])
        set(gca,'XTickLabel',labels)
        set(gca,'Layer','top')

    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'boxplot-mag.png') 
    end 

    figure(2), clf, hold all
        set(gcf,'color','w')
        color = [0, 93, 141]/255;
        cmap = parula(64);
        
        m = dist0.groundTruth{1}.nBasisFunctionsEst;
        d = dist0.groundTruth{1}.d;
        LL = dist0.groundTruth{1}.LL;
        pos = dist0.groundTruth{1}.pos; % Since ground truth position always the same
        x1t = dist0.groundTruth{1}.fullMap_x1t;
        x2t = dist0.groundTruth{1}.fullMap_x2t;
        [X1t,X2t] = meshgrid(x1t,x2t);
        xt = [X1t(:) X2t(:) 0*X1t(:)];
        nt = length(xt);
        [eigenval,~,eigenfun_dx,NN] = domain_cartesian_dx(m,d,LL);
        dPhixt = eigenfun_dx(NN,xt,1);
        dPhiyt = eigenfun_dx(NN,xt,2);
        dPhizt = eigenfun_dx(NN,xt,3);
        dPhixt = [ones(nt,1), zeros(nt,2), dPhixt];
        dPhiyt = [zeros(nt,1), ones(nt,1), zeros(nt,1), dPhiyt];
        dPhizt = [zeros(nt,2), ones(nt,1), dPhizt];
        %     
        cmax = max(sqrt(sum(dist0.groundTruth{1}.fullMap_df.^2,2)));
        cmin = min(sqrt(sum(dist0.groundTruth{1}.fullMap_df.^2,2)));
        h=imagesc(x1t,x2t,reshape(sqrt(sum(dist0.groundTruth{1}.fullMap_df.^2,2)),size(X1t)));
        plot(dist0.groundTruth{1}.pos(1,:)',dist0.groundTruth{1}.pos(2,:)','color',color,'linewidth',1)
        colorbar
        caxis([cmin cmax])
        axis equal
        xlim([min(xt(:,1)) max(xt(:,1))])
        ylim([min(xt(:,2)) max(xt(:,2))])
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);

    if savePlots
      frame = getframe(gca);
      imwrite(frame.cdata,'mag-path-field.png') 
    end 
end