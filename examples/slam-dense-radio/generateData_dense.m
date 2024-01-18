function [dx, initState, y, groundTruth] = generateData_dense(params,dynModel,fieldData)
%% GENERATEDATA_DENSE - Generates dense data for simulation experiments
%
% Syntax:
%   [dx, initState, y, groundTruth] = generateData_dense(params, dynModel, fieldData)
%
% In:
%   params    - Struct containing various parameters for data generation
%   dynModel  - Function handle for the dynamic model
%   fieldData - Field data for the simulation (optional)
%
% Out:
%   dx         - Simulated odometry data
%   initState  - Initial state for the simulation
%   y          - Simulated sensor measurements
%   groundTruth - Struct containing ground truth data for position, 
%                 orientation, odometry, and other relevant information
%
% Description:
%   This function generates simulated data for dense experiments, including
%   odometry data, sensor measurements, and ground truth information. The 
%   function uses provided parameters and dynamic models, and can simulate 
%   data for various trajectory types such as 2D circles, 3D squares, and 
%   others. It can also generate field data if not provided.
%
% Copyright:
%   2023-   Manon Kok and Arno Solin

%% Generate odometry data 
% Parameters should either be input via params struct or standard inputs
% are used 
if nargin >= 1   
    % If file run separately then extract info
    makePlots = params.makePlots;
    trajType = params.trajType;
    theta = params.theta;
else
    % Else use some standard settings
    makePlots = 1;
    trajType = 'square_3D';
    theta = [0.5 ; 1 ; 0.01];
    
    params.dt = 1;

    N = 128;
    params.N = N;
    Q = 1E-6 * ones(N,1);
    Q(N/8 + N/8*(0:6)) = 0.2^2;
    Q = reshape(Q,[1 1 N]);
    params.Q = Q;
    dynModel = @(xn,dx,dt,Q) [xn(1:2) + [cos(xn(3)), -sin(xn(3)) ; ...
                sin(xn(3)), cos(xn(3))]' * dx(1:2)' ; xn(3) + dx(3) + sqrt(Q) * randn]; 
end

% If we have a 6D trajectory, we simulate a 3D curl-free field; if the
% trajectory is 2D (planar) or 3D (planar + heading), we only simulate a 1D
% field from a squared exponential
if strcmp(trajType(end-1:end),'6D')
    field3D = 1;
else 
    field3D = 0;
end

%% Simulate groundtruth data and noiseless odometry for different cases
groundTruth = [];
switch trajType
    case 'circle_2D'
        disp('Generating 2D data on circle')
        radius = 2; % Radius of circular trajectory
        nLaps = 3; % Number of laps
        dpsi = 5; % Obtain odometry sample every 5 degrees of a circle
        psi = (0:dpsi:360 * nLaps - dpsi)*pi/180;
        N = length(psi);

        pos = [radius * cos(psi); radius * sin(psi)]; % True position
        groundTruth.pos = pos; % Save ground truth data
        initState = pos(:,1); % Initial state
        dx = diff(pos'); % Resulting odometry
    case 'bean_2D'
        disp('Generating bean-shaped 2D data')
        nLaps = 3; % Number of laps
        nDataPointsPerLap = 63; % Number of points per lap
        a = 5;
        psi = linspace(0,pi,nDataPointsPerLap)';
        r = a*sin(psi).^3 + a*cos(psi).^3;
        u = r.*cos(psi)-.3;
        v = r.*sin(psi)-.3;
%         yu = diff(u);
%         yv = diff(v);
%         th = atan2(yv,yu);
        pos = [u'; v'];
        % Center positions around zero
        pos = pos - mean([min(pos,[],2) , max(pos,[],2)],2); 
        % Make multiple laps
        pos = [pos , repmat(pos(:,2:end),[1 nLaps-1])]; % True position
        
        groundTruth.pos = pos; % Save ground truth data
        initState = pos(:,1); % Initial state
        dx = diff(pos'); % Resulting odometry
        N = length(pos);
    case 'square_3D'
        disp('Generating 3D data on square')    
        N = params.N;
        pos = [zeros(1,N/4), linspace(0,2,N/4), 2*ones(1,N/4), linspace(2,0,N/4); ...
            linspace(0,2,N/4), 2*ones(1,N/4), linspace(2,0,N/4), zeros(1,N/4)];
        % Center positions around zero
        pos = pos - mean(pos,2);
        
        groundTruth.pos = pos; % Save ground truth data
        initState = [pos(:,1);0]; % Initial state
        dx = [diff(pos'),zeros(N-1,1)]; % Resulting odometry
    case {'line_3D','line_2D','line_3D_withPos'}
        disp('Generating data on line')  
        % Simulate position data
        if isfield(params,'N')
            N = params.N;
        else 
            N = 32;
        end
        pos = [zeros(1,N); ...
            [linspace(0,3,N/2) , ...
            linspace(3,0,N/2)]];
        % Center positions around zero
        pos = pos - mean(pos,2);
        
        groundTruth.pos = pos; % Save ground truth data
        initState = pos(:,1); % Initial position
        dx = diff(pos'); % Resulting odometry
        if strcmp(trajType,'line_3D') || strcmp(trajType,'line_3D_withPos')
            initState = [initState ; 0];
            dx = [dx, zeros(N-1,1)];
        end
    case 'line_6D'
        disp('Generating 6D data on line')
        N = 32;
        pos = [zeros(1,N); ...
            [linspace(0,3,N/2) , ...
            linspace(3,0,N/2)] ; ...
            zeros(1,N)];
        pos = pos - mean(pos,2);
        quat = [[ones(N/2,1),zeros(N/2,3)] ; ...
            [zeros(N/2,3),-ones(N/2,1)]];

        % Ground truth
        groundTruth.pos = pos;
        groundTruth.quat = quat;

        % Odometry measurements
        initState = [pos(:,1) ; quat(1,:)'];
        dQuat = squeeze(multiprod(qLeft(qInv(quat(1:end-1,:))), ...
            reshape(quat(2:end,:)',[4 1 length(quat)-1])))';

        dPos = diff(pos');
        dx = [dPos , dQuat];
        N = length(pos);
    case 'circle_6D'
        disp('Generating 6D data on circle')
        dtheta = 5;
        radius = 2;
        nLaps = 2;
        psi = (0:dtheta:360-dtheta);
        psi = repmat(psi,1,nLaps);
        N = length(psi);
        
        % True position and orientation
        pos = [radius * cos(psi*pi/180); radius * sin(psi*pi/180); zeros(1,length(psi))];
        R = [reshape(cos(psi*pi/180),[1 1 N]), reshape(sin(psi*pi/180),[1 1 N]), zeros(1,1,N) ; ...
                reshape(-sin(psi*pi/180),[1 1 N]), reshape(cos(psi*pi/180),[1 1 N]), zeros(1,1,N) ; ...
                zeros(1,1,N), zeros(1,1,N), ones(1,1,N)];
        quat = rotm2quat(R);

        groundTruth.pos = pos;
        groundTruth.quat = quat;

        % Odometry measurements
        initState = [pos(:,1) ; quat(1,:)'];
        dPos = diff(pos');
        dQuat = squeeze(multiprod(qLeft(qInv(quat(1:end-1,:))), ...
            reshape(quat(2:end,:)',[4 1 length(quat)-1])))';
        dx = [dPos, dQuat];
    case 'bean_6D'
        disp('Generating bean-shaped 6D data')
        % Simulate position data
        nLaps = 3;
        nDataPointsPerLap = 64;
        a = 15;
        psi = linspace(0,nLaps*pi,nLaps * nDataPointsPerLap);
        r = a*sin(psi).^3 + a*cos(psi).^3;
        u = r.*cos(psi)-.3;
        v = r.*sin(psi)-.3;
        yu = diff(u);
        yv = diff(v);
        th = atan2(yv,yu); th = [th, th(end)];
        pos = [u; v ; zeros(1,length(u))];
        N = length(pos);
        R = [reshape(cos(th),[1 1 N]), reshape(sin(th),[1 1 N]), zeros(1,1,N) ; ...
                reshape(-sin(th),[1 1 N]), reshape(cos(th),[1 1 N]), zeros(1,1,N) ; ...
                zeros(1,1,N), zeros(1,1,N), ones(1,1,N)];
        %quat = rotm2quat(R);
        quat = rmat2quat(R)';
        % Center positions around zero
        pos = pos - mean([min(pos,[],2) , max(pos,[],2)],2); 
        
        % Save ground truth data
        groundTruth.pos = pos; 
        groundTruth.quat = quat;
        
        % Odometry measurements
        initState = [pos(:,1) ; quat(1,:)']; % Initial state
        dPos = diff(pos'); % Resulting odometry
        dQuat = squeeze(multiprod(qLeft(qInv(quat(1:end-1,:))), ...
            reshape(quat(2:end,:)',[4 1 length(quat)-1])))';
        dx = [dPos, dQuat];
end

%% Simulate field data in navigation frame
m = 2000; % Simulate with many basis functions
% Parameter to define how many length scales the domain should extend
% outside of movement area
if isfield(params,'nLL')
    nLL = params.nLL; % The parameter can be specified
else
    nLL = 2; % or is by default set to 2
end
if field3D
    lengthScale = theta(2); % Length scale parameter 
    % Make 3D domain sufficiently larger around trajectory
    LL = [min(pos(1,:)) - nLL*lengthScale , max(pos(1,:)) + nLL*lengthScale ; ...
        min(pos(2,:)) - nLL*lengthScale, max(pos(2,:)) + nLL*lengthScale; ...
        -nLL*lengthScale, nLL*lengthScale]';
    
    % Make test points
    x1t = linspace(LL(1,1),LL(2,1),100);
    x2t = linspace(LL(1,2),LL(2,2),100);
    [X1t,X2t] = meshgrid(x1t,x2t);
    xt = [X1t(:) X2t(:) 0*X1t(:)];
    
    % Simulate function and measurements both at locations visited as well
    % as at testing points
    [f,df,yn] = gp_rnd_scalar_potential_fast([pos' ; xt],m,LL,theta);
    
    % Save complete map and grid
    groundTruth.fullMap_f = f(length(pos)+1:end); 
    groundTruth.fullMap_df = df(length(pos)+1:end,:); 
    groundTruth.fullMap_x1t = x1t;
    groundTruth.fullMap_x2t = x2t;
    % And save ground truth function values along trajectory
    groundTruth.df = df(1:length(pos),:);
    
    % We only need the measurements along the trajectory
    yn = yn(1:length(pos),:); 
    % Rotate magnetometer data to body frame
    y = zeros(N,3);
    for i = 1:N
        Ri = quat2rmat(quat(i,:));
        y(i,:) = Ri' * yn(i,:)';
    end
else
    lengthScale = theta(1); % Length scale parameter 
    % Make 2D domain sufficiently larger around trajectory
    LL = [min(pos(1,:)) - nLL*lengthScale , max(pos(1,:)) + nLL*lengthScale ; ...
        min(pos(2,:)) - nLL*lengthScale, max(pos(2,:)) + nLL*lengthScale]';
    % The third input argument is field data; We use this option if we
    % don't want to simulate a new field but new measurements and new
    % odometry but keep the same field
    if nargin < 3
        % Make test points
        x1t = linspace(LL(1,1),LL(2,1),100);
        x2t = linspace(LL(1,2),LL(2,2),100);
        [X1t,X2t] = meshgrid(x1t,x2t);
        xt = [X1t(:) X2t(:)];
        
        % Simulate function and measurements both at locations visited as
        % well as at testing points
        [f,y] = gp_rnd_SE1D_fast([pos' ; xt],m,LL,theta);
        
        % Save complete map and grid
        groundTruth.fullMap_f = f(length(pos)+1:end); % Save complete map and grid
        groundTruth.fullMap_x1t = x1t;
        groundTruth.fullMap_x2t = x2t;

        % We only need the measurements along the trajectory
        y = y(1:length(pos)); 
    else
        f = fieldData.f; % Keep same map
        y = f + sqrt(theta(3))*randn(size(f)); % But add new measurement noise
    end
    % Save ground truth function values along trajectory
    groundTruth.f = f(1:length(pos));
end
groundTruth.LL = LL;
groundTruth.nBasisFunctionsSim = m;

%% Add measurement noise to odometry by using the dynamic model to put it 
% on the way that is assume in the model
Q = params.Q;
if size(Q,3) == 1 % Allow for both time-varying and constant Q
      Q = repmat(Q,[1 1 N]);
end
% Run dynamic model forward 
dt = params.dt;
if field3D
    x = zeros(N,length(initState));
    dQuat = zeros(N-1,4);
    x(1,:) = initState;
    for i = 2:N
        [x(i,:),dQuat(i-1,:)] = dynModel(x(i-1,:)',dx(i-1,:),dt,Q(:,:,i-1));
    end
    dx = [diff(x(:,1:3)),dQuat];
else
    x = zeros(N,length(initState));
    x(1,:) = initState;
    for i = 2:N
        x(i,:) = dynModel(x(i-1,:)',dx(i-1,:),dt,Q(:,:,i-1));
    end

    switch trajType
        case {'line_3D','square_3D'}
            dx = [dx(:,1:2), diff(x(:,3))];
        case {'line_2D','bean_2D','circle_2D'}
            dx = diff(x);
    end
end
groundTruth.odometry = x;
groundTruth.Q = Q;

%% Visualize ground truth and odometry
if makePlots && params.visualiseResults
  figure(1); cla; hold on
      if field3D
          imagesc(x1t,x2t,reshape(sqrt(sum(df(length(pos)+1:end,:).^2,2)),size(X1t)));
          colorbar
          caxis([0 max(sqrt(sum(yn.^2,2)))+10])
      else
          imagesc(x1t,x2t,reshape(f(length(pos)+1:end),size(X1t)));
          colorbar
          caxis([min(f) max(f)])
      end
      plot(groundTruth.pos(1,:),groundTruth.pos(2,:),'k','LineWidth',2)
      plot(groundTruth.odometry(:,1),groundTruth.odometry(:,2),'r','LineWidth',2)
      axis equal
      xlim([min(xt(:,1)) max(xt(:,1))])
      ylim([min(xt(:,2)) max(xt(:,2))])
      title('True map, true trajectory, odometry')
end

end