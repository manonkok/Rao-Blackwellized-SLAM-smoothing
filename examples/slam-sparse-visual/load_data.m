function [Y,u,map,groundTruth,initPos,initTheta,f,fp,fw,nonnanswaps] = ...
  load_data(n_shuffle,posVar,posBias,angleVar)
%% LOAD_DATA - Load and preprocess simulation data for experiments
%
% Syntax:
%   [Y,u,map,groundTruth,initPos,initTheta,f,fp,fw,nonnanswaps] = 
%       load_data(n_shuffle,posVar,posBias,angleVar)
%
% In:
%   n_shuffle   - Number of times to shuffle the data (default: 0)
%   posVar      - Variance of position noise (default: 0.04^2)
%   posBias     - Bias in position measurements (default: 0.01)
%   angleVar    - Variance of angle noise (default: (.001^2)^2)
%
% Out:
%   Y           - Noisy observations
%   u           - Control inputs with noise and drift
%   map         - Map of the environment (if available)
%   groundTruth - Ground truth of the robot's path
%   initPos     - Initial position of the robot
%   initTheta   - Initial orientation of the robot
%   f           - Focal length of the camera model
%   fp          - Principal point of the camera model
%   fw          - Image width of the camera model
%   nonnanswaps - Number of valid swaps in data (if applicable)
%
% Description:
%   Load and preprocess simulation data for robotic experiments. The function 
%   sets up a camera model, loads a predefined path and simulated observations, 
%   and adds noise and bias to the control inputs and observations. 
%   It can also shuffle the data if required.
%
% Examples:
%   [Y,u,map,gt,ip,it,f,fp,fw,nns] = load_data();
%   [Y,u,~,~,~,~,~,~,~,~] = load_data(10, 0.05^2, 0.02, (.002^2)^2);
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

%% Defaults

  if nargin<1 || isempty(n_shuffle)
    n_shuffle = 0;
  end
  if nargin<2 || isempty(posVar)
    posVar = 0.04^2;
  end
  if nargin<3 || isempty(posBias)
    posBias = .01;
  end
  if nargin<4 || isempty(angleVar)
    angleVar = (.001^2)^2;
  end
  
%% Set camera model

  % Camera matrix
  f = 1.5; % Focal length
  fp = 0;  % Principal point
  fw = 1;  % Image width
  fov = 2*atan2(fw,f);


%% Load the path and simulated observations

  load('curve-x2.mat')
  

%% Observations

  %rng(0,'twister')

  initPos = p(:,1);
  initTheta = th(1);
  
  dPos = diff(p,[],2);
  dTheta = diff(unwrap(th))';
  T = 1; % = dt
  
  u = [dPos; dTheta]';
  
  % Add noise + drift
  u(:,1:2) = u(:,1:2)+sqrt(posVar)*randn(size(u(:,1:2))) + posBias;
  u(:,3) = u(:,3)+sqrt(angleVar)*randn(size(u(:,3)));
  
  % Add noise
  Y = Yclean + .01*randn(size(Y));

  groundTruth = [p; th(:)'];
  
  % Naive reconstruction
  %figure(1); clf; hold on
  
    % Background
    %image(xi,yi,(255-.3*alpha))
    %colormap(gray), caxis([0 255])
    
    %foo = cumsum(u);
    %plot(initPos(1)+foo(:,1),initPos(2)+foo(:,2));
  
    %plot_map([initPos; initTheta; map(:)],f,fp,fw,cmap);
    
    %axis off
    %set(gcf,'color','w')
    
    %frame = getframe(gca); 
    %imwrite(frame.cdata,'bean-odometry.png')
    
%% Shuffle data observations (corrupt them)    

  nonnanswaps = 0;

  if n_shuffle>0
    
    % Shuffle some of the observations
    t_rand = sort(randi(size(Y,2),1,n_shuffle));
    for i=t_rand
        i_swap = randi(size(Y,1)/2-1);
        foo = Y(:,i);
        tmp = foo(i_swap);
        foo(i_swap) = foo(i_swap+1);
        foo(i_swap+1) = tmp;
        Y(:,i) = foo(:);
        %fprintf('t=%i - swapped %i and %i: %.2f <-> %.2f\n',i,i_swap,i_swap+1, ...
        %  foo(i_swap), foo(i_swap+1))
        if ~isnan(foo(i_swap)) || ~isnan(foo(i_swap+1)), nonnanswaps=nonnanswaps+1; end
    end
    
    %nonnanswaps / sum(~isnan(Y(:)))
    
  end