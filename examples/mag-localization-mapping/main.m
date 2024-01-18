%% Magnetic localization / mapping examples (not SLAM)
%
% Compose the parts required for the figure in the paper that shows
% magnetic GP maps as well as runs a particle filter for magnetic terrain
% matching. Simplified example that does not do any re-initialization,
% smart initialization, etc.
%
% References:
%   [1] Arno Solin, Manon Kok, Niklas Wahlström, Thomas B. Schön, and 
%       Simo Sarkka (2018). Modeling and interpolation of the ambient 
%       magnetic field by Gaussian processes. IEEE Transactions on 
%       Robotics (T-RO), 34(4):1112–1127.
%   [1] Arno Solin, Simo Sarkka, Juho Kannala, and Esa Rahtu (2016). 
%       Terrain navigation in the magnetic landscape: Particle filtering 
%       for indoor positioning. In Proceedings of the European Navigation 
%       Conference (ENC). Pages 1–9. IEEE.
%   [2] Manon Kok and Arno Solin (2018). Scalable magnetic field SLAM in 
%       3D using Gaussian process maps. In Proceedings of the International 
%       Conference on Information Fusion (FUSION). Pages 1353–1360. 
%       Cambridge, UK.
%
%% Load data

  clear, close all, clc

  % Check if data is downloaded
  if isfolder('magnetic-data')==false
    error('Follow instructions to download the magnetic data first.')
  end

  % Data folder cloned from the public GitHub data repository
  datapath = './magnetic-data';
  
  % Sensor (invensense/nexus/trivisio)
  sensor = 'invensense';
  
  % Allocate space
  x = [];
  y = [];
  t = [];
  s = [];
  
  % Loop through
  for i=1:9
    
      % Read in files
      xi = load(fullfile(datapath,'data',sensor, ...
        sprintf('%i-loc.csv',i)));
      yi = load(fullfile(datapath,'data',sensor, ...
        sprintf('%i-mag.csv',i)));
      ti = load(fullfile(datapath,'data',sensor, ...
        sprintf('%i-time.csv',i)));
      
      % Concatenate
      x = [x; xi]; %#ok
      y = [y; yi]; %#ok
      t = [t; ti]; %#ok
      s = [s; i+0*ti]; %#ok
     
  end
  
  % Visualize all the separate paths
  figure(1); clf
  for i=1:9
    subplot(3,3,i)
      ii = (s==i);
      scatter(x(ii,1),x(ii,2),6,y(ii,1))
      axis equal 
  end
  
  % Visualize the combined data
  figure(2); clf; hold on
  for i=1:3
      ii = (s==i);
      scatter(x(ii,1),x(ii,2),6,y(ii,1))
      axis equal 
  end
  
  data.x = x;
  data.y = y;
  data.t = t;
  data.s = s;

  
%% GP model for the magnetic field

  addpath ../../tools
  
  % Number of basis functions
  m = 1000;

  % Dimensionality
  d = 3;
  
  % Training data
  ind = data.s<4;
  x = data.x(ind,:);
  y = data.y(ind,:);
  x = [x 0*x(:,1)];  
  
  % Domain boundaries
  rangex = range(x);
  pm = 0.2*min(rangex(rangex>0));
  LL = [min(x)-pm; max(x)+pm];
  
  % Test points
  xi = linspace(LL(1,1),LL(2,1),100);
  yi = linspace(LL(1,2),LL(2,2),100);
  [XI,YI]=meshgrid(xi,yi);
  xt = [XI(:) YI(:) 0*YI(:)];
  
  % Hyperparameters
  theta = [500.0000 0.36 26 1];
  opt = [false true true true];
  
  % Run GP learning and inference
  [Eft,dEft,Varft,theta,lik,dVarft] = gp_scalar_potential_fast(x,y,xt,m,LL,theta,opt);


%% Set up homography overlay

  % Corners
  corner = [1 1; -1 1; -1 -1; 1 -1];

  % Markers
  markers = ...
  [ 4.7741   -3.9192    0.1101
    4.9283    1.2333    0.0896
   -1.2583    1.0878   -0.0092
   -1.2033   -2.7739    0.0136]; 
  
  % Grid
  z1 = [linspace(-1,1,32) nan];
  z2 = [linspace(-1,1,32) nan];
  [G1,G2] = meshgrid(z1,z2);
  Z = [G1(:) G2(:)];
  
  % Mapping
  [A,C] = homography_estimation(markers(:,1:2),corner);

  % Transform
  G = bsxfun(@rdivide,A*[Z ones(size(Z,1),1)]',C*[Z ones(size(Z,1),1)]')';
  
  % Discretize
  z1 = linspace(-1.3,1.3,64);
  z2 = linspace(-1.3,1.3,64);
  [Z1,Z2] = meshgrid(z1,z2);
  Z = [Z1(:) Z2(:)];

  % Transform
  X = bsxfun(@rdivide,A*[Z ones(size(Z,1),1)]',C*[Z ones(size(Z,1),1)]')';

  % Test points
  X1 = reshape(X(:,1),size(Z1));
  X2 = reshape(X(:,2),size(Z2));
  xt = X(:,1:2);
  
  
%% Visualization transformation

  % Markers on image
  im_markers = [1314,130; 279,213; 139,1017; 1798,749];
  
  % Mapping
  [A,C] = homography_estimation(im_markers,markers(:,1:2));
    
  % Show grid
  Y = bsxfun(@rdivide,A*[G(:,1:2) ones(size(G,1),1)]',C*[G(:,1:2) ones(size(G,1),1)]')';
  G1 = reshape(Y(:,1),size(G1));
  G2 = reshape(Y(:,2),size(G2));

  % Transform test coordinates
  Y = bsxfun(@rdivide,A*[X(:,1:2) ones(size(X,1),1)]',C*[X(:,1:2) ones(size(X,1),1)]')';
  Y1 = reshape(Y(:,1),size(Z1));
  Y2 = reshape(Y(:,2),size(Z2));
  

%% Predict overlay
  
  % Training data for magnetic field
  ind = data.s<3 | data.s==4;
  x = data.x(ind,:);
  y = data.y(ind,:);
  x = [x 0*x(:,1)];  
  
  % Thinning
  x = x(1:10:end,:);
  y = y(1:10:end,:);
  
  % Mapping path
  path = bsxfun(@rdivide,A*[x(:,1:2) ones(size(x,1),1)]',C*[x(:,1:2) ones(size(x,1),1)]')';
  
  % Localization path
  ind = data.s==3;
  xl = data.x(ind,:);
  path_loc = bsxfun(@rdivide,A*[xl(:,1:2) ones(size(xl,1),1)]',C*[xl(:,1:2) ones(size(xl,1),1)]')';
  
  % GP prediction
  [Eft,dEft,Varft,theta,lik,dVarft] = ...
      gp_scalar_potential_fast(x,y,[xt 0*xt(:,1)],m,LL,theta,false);


%% Save plots of GP magnetic field maps

  savePlot = false;

  % The magnitude field
  figure(1); clf; hold on
    [clim]=plotgpmap(dEft,dVarft,G1,G2,Y1,Y2,[nan nan nan],y,'frame.png',1:3);
    plot(path(:,1),path(:,2),'--','LineWidth',.5,'color',[.5 .5 .5])
  
  drawnow
  
  % Save figure
  if savePlot
    frame = getframe(gca);
    imwrite(frame.cdata,sprintf('robot-map-mag.png'))
  end
  
  % Components
  for j=1:3
    figure(1+j); clf
      plotgpmap(dEft,dVarft,G1,G2,Y1,Y2,[nan nan nan],y,[],j);
      caxis([-60 30])
      drawnow
      
    % Save figure
    xyz = {'x','y','z'};
    if savePlot
      frame = getframe(gca);
      imwrite(frame.cdata,sprintf('robot-map-%s.png',xyz{j})) 
    end  
  end

  
%% Run the localization example that uses a magnetic map

  % One of the paths is used as a test path 
  % (ground truth visualized by a dashed line)

  % Run localization script (see options in script)
  run_localization(data,false,true)

  