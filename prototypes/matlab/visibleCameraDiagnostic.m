function ST = visibleCameraDiagnostic(filename,camPos,hitPointPos,camParams,optics)
% camParams = [fieldOfViewCone, numAngleViews, numRotAngles, samplingPoints, toroidalSlices]
% camPos = [poloidal direction in degrees, toroidal direction in degrees]
% hitPointPos = [poloidal direction in degrees, toroidal direction in degrees]
close all

ST = struct;
ST.filename = filename; % name of file containing the relevant data
% ST.camPos = camPos;
% ST.hitPointPos = hitPointPos;
% ST.camParams = camParams;
% ST.optics = optics;


ST.data = loadData(ST);

ST.stats = statistics(ST,{'k','pitch','angularMomentum'});


% EXAMPLES

% ST.SP = fastCamera(ST);
end

function data = loadData(ST)
disp('Loading data...')

% Fundamental constants
Eo = 3E6; % MeV
c = 2.9979E8; % Speed of light, in m/s
qe = 1.602176E-19; % Electron charge, in Coulombs
me = 9.109382E-31; % Electron mass, in kg.
vo = sqrt( 1 - (me*c^2/(Eo*qe))^2 );
gamma = 1/sqrt(1 - vo^2);
K = 8.987E9; % In Nm^2/C^2
% Fundamental constants

NR = 60;
NZ = 60;

data = struct;

pitch = [10];
RZ = [];
k = [];
for pp=1:numel(pitch)
    data.raw = load([ST.filename '_po_' num2str(pitch(pp)) '.mat']);
    RZ = vertcat(RZ,data.raw.RZ{:});
    k = vertcat(k,data.raw.k{:});
end

% Here we calculate the radiated power in a poloidal plane using the
% curvature.

% RZ = vertcat(data.raw.RZ{:});

Rmin = min(RZ(:,1));
Rmax = max(RZ(:,1));

Zmin = min(RZ(:,2));
Zmax = max(RZ(:,2));

% k = horzcat(data.raw.k{:});

DR = (Rmax - Rmin)/NR;
DZ = (Zmax - Zmin)/NZ;

R_vertices = Rmin:DR:Rmax;
Z_vertices = Zmin:DZ:Zmax;

R_nodes = R_vertices(1:end-1) + 0.5*DR;
Z_nodes = Z_vertices(1:end-1) + 0.5*DZ;

polPlane = zeros(NZ,NR);
count = zeros(NZ,NR);

for ii = 1:numel(k)
    try
        ind_R = floor( (RZ(ii,1) - Rmin)/DR ) + 1;
        ind_Z = floor( (RZ(ii,2) + abs(Zmin))/DZ ) + 1;
        
        if ind_R == NR + 1
            ind_R = ind_R - 1;
        end
        
        if ind_Z == NZ + 1
            ind_Z = ind_Z - 1;
        end

        polPlane(ind_Z,ind_R) = polPlane(ind_Z,ind_R) + k(ii)^2;
        count(ind_Z,ind_R) = count(ind_Z,ind_R) + 1;
    catch
        disp(['Index in R: ' num2str(ind_R)])
        disp(['Index in Z: ' num2str(ind_Z)])
    end
    
end

polPlane = (2*K*qe^2*gamma^4*c/3)*polPlane;
polPlane = polPlane./count;

% polPlane(~isfinite(polPlane)) = min(min(polPlane));

surfc(R_nodes,Z_nodes,polPlane,'LineStyle','none')
colormap(jet)
box on; grid on; axis square
% shading interp
xlabel('R','Interpreter','latex','FontSize',16)
ylabel('Z','Interpreter','latex','FontSize',16)

data.polPlane = polPlane;
data.R_nodes = R_nodes;
data.Z_nodes = Z_nodes;

disp('Data loaded!')
end

function stats = statistics(ST,listOfParams)
disp('Starting statistical analysis')

% Fundamental constants
Eo = 3E6; % MeV
c = 2.9979E8; % Speed of light, in m/s
qe = 1.602176E-19; % Electron charge, in Coulombs
me = 9.109382E-31; % Electron mass, in kg.
vo = sqrt( 1 - (me*c^2/(Eo*qe))^2 );
gamma = 1/sqrt(1 - vo^2);
K = 8.987E9; % In Nm^2/C^2
% Fundamental constants

cadence = 2E4;
DT = 1E-2;
numTracers = numel(ST.data.raw.timeSeries);

numBins = 50;
stats = struct;

for pp =1:numel(listOfParams)
    numSnapshots = numel(ST.data.raw.timeSeries{1}.(listOfParams{pp}));
    param = zeros(numTracers,numSnapshots);
    
    if strcmp(listOfParams{pp},'k')
        for ii=1:numTracers
            param(ii,:) = (2*K*qe^2*gamma^4*c/3)*ST.data.raw.timeSeries{ii}.(listOfParams{pp}).^2;
        end
    elseif strcmp(listOfParams{pp},'pitch')
        for ii=1:numTracers
            param(ii,:) = 180*ST.data.raw.timeSeries{ii}.(listOfParams{pp})/pi;
        end
    else
        for ii=1:numTracers
            param(ii,:) = ST.data.raw.timeSeries{ii}.(listOfParams{pp});
        end
    end
    
    min_val = min(min(param));
    max_val = max(max(param));
    edges = linspace(min_val, max_val,numBins);
    
    stats.pdfs.(listOfParams{pp}).vals = zeros(numBins-1,numSnapshots);
    stats.pdfs.(listOfParams{pp}).edges = edges;
    
    Dbin = mean(diff(stats.pdfs.(listOfParams{pp}).edges));
    
    f = figure;
    hold on
    for it=1:numSnapshots
        [h,~] = histcounts(param(:,it),edges);
        h = h/(Dbin*sum(h));
        stats.pdfs.(listOfParams{pp}).vals(:,it) = h;
        plot3(edges(1:numBins-1),DT*cadence*it*ones(1,numBins-1),h,'k')
    end
    hold off
    xlim([min_val, max_val])
    ylim([1, DT*cadence*numSnapshots])
    box on; grid on
    xlabel([listOfParams{pp} '(t)'],'Interpreter','latex','FontSize',16)
    ylabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    zlabel('Number of $e^{-}$','Interpreter','latex','FontSize',16)
    
    if strcmp(listOfParams{pp},'k')
        rp = zeros(1,numSnapshots); % instantaneous radiated power
        rp_pdfs = zeros(1,numSnapshots); % instantaneous radiated power (using pdfs)
        
        Dbin = mean(diff(stats.pdfs.(listOfParams{pp}).edges));
        P = stats.pdfs.(listOfParams{pp}).edges(1:numBins-1) + 0.5*Dbin;
        
        for it=1:numSnapshots
            rp(it) = mean(param(:,it)); 
            rp_pdfs(it) = ...
                Dbin*sum( P.*stats.pdfs.(listOfParams{pp}).vals(:,it)' );
        end
        
        figure
        plot(DT*cadence*(1:1:numSnapshots),rp,'k',DT*cadence*(1:1:numSnapshots),rp_pdfs,'r')
        grid on; box on
        xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
        ylabel('Total radiated power []','Interpreter','latex','FontSize',16)
        
    end
end





disp('Statistics done!')
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Camera diagnostic

function LOS = calculateLinesOfSight(ST,posCam,posHitPoint,maxAngleView)
disp('Calculating lines of sight...')
LOS = struct;

% Define position of the camera
cam = struct;

% The camera is placed at the given poloidal and toroidal angles, right
% at the outer boundary of the simulation domain
if (posCam(1)>225) && (posCam(1)<300) || (posHitPoint(1)>225) && (posHitPoint(1)<300)
    error('Camera cannot be placed at the specified location.')
else
    polAngle_cam = posCam(1);% poloidal angle in degrees w.r.t. the outer midplane plasma.
    torAngle_cam = posCam(2);% toroidal angle in degrees.
end

polAngle_hitPoint = posHitPoint(1);
torAngle_hitPoint = posHitPoint(2);

% transform to radians
polAngle_cam = pi*polAngle_cam/180;
polAngle_hitPoint = pi*polAngle_hitPoint/180;

torAngle_cam = pi*torAngle_cam/180;
torAngle_hitPoint = pi*torAngle_hitPoint/180;
% transform to radians

ST = innerBoundary(ST); % figure handle of 3D torus

cam.position = zeros(1,3);% (x,y,z)
cam.hitPoint = zeros(1,3);% (x,y,z)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Camera and hit point position

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% r and poloidal angle (pol) defined with respect to R0 and Z0.
wall_r = sqrt((ST.geom.wall.R - ST.geom.R0).^2 + ...
    (ST.geom.wall.Z - ST.geom.Z0).^2);
wall_pol = atan2(ST.geom.wall.Z - ST.geom.Z0,ST.geom.wall.R - ST.geom.R0);
% r and theta (pol) defined with respect to R0 and Z0.

I = wall_pol < 0;
wall_pol(I) = wall_pol(I) + 2*pi;
clear I

[~,I] = min(wall_pol);
wall_pol = [wall_pol(I:end);wall_pol(1:I-1)];
wall_r = [wall_r(I:end);wall_r(1:I-1)];
R = [ST.geom.wall.R(I:end);ST.geom.wall.R(1:I-1)];
Z = [ST.geom.wall.Z(I:end);ST.geom.wall.Z(1:I-1)];
clear I

REP = [];
for ii=1:length(R)-1
    logic1 = R(ii) == R(ii+1);
    logic2 = Z(ii) == Z(ii+1);
    if logic1 && logic2
        REP = [REP;ii];
    end
end

wall_r(REP) = [];
wall_pol(REP) = [];

R(REP) = [];
Z(REP) = [];

% Flip if wall_r is clockwise
rotDir = mean(diff(wall_pol)./abs(diff(wall_pol))); % "> 0" anti-clockwise, "> 0" clockwise
if rotDir < 0
    wall_r = flip(wall_r);
    wall_pol = flip(wall_pol);
    R = flip(R);
    Z = flip(Z);
end
% Flip if wall_r is clockwise
[~,I] = min(wall_pol);
wall_pol = [wall_pol(I:end);wall_pol(1:I-1)];
wall_r = [wall_r(I:end);wall_r(1:I-1)];
R = [R(I:end);R(1:I-1)];
Z = [Z(I:end);Z(1:I-1)];
clear I

wall_r = [wall_r(end); wall_r; wall_r(1)];
wall_pol = [wall_pol(end)-2*pi; wall_pol; wall_pol(1)+2*pi];
R = [R(end); R; R(1)];
Z = [Z(end); Z; Z(1)];


[~,I] = min(abs(wall_pol - polAngle_cam));

cam.position(1) = R(I)*cos(torAngle_cam);
cam.position(2) = -R(I)*sin(torAngle_cam);
cam.position(3) = Z(I);
clear I

[~,I] = min(abs(wall_pol - polAngle_hitPoint));
cam.hitPoint(1) = R(I)*cos(torAngle_hitPoint);
cam.hitPoint(2) = -R(I)*sin(torAngle_hitPoint);
cam.hitPoint(3) = Z(I);
clear I

figure(ST.IB.figh)
hold on
plot3(cam.position(1),cam.position(2),cam.position(3),'ks',...
    'MarkerSize',10,'MarkerFaceColor',[0 0 0])
plot3(cam.hitPoint(1),cam.hitPoint(2),cam.hitPoint(3),'ko',...
    'MarkerSize',10,'MarkerFaceColor',[0 0 0])
hold off
% Camera and hit point position
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

LOS.main.camera = cam;

% main line of sight
mainLineOfSight = zeros(ST.camParams(4),3);% (x,y,z)
mainLineOfSight(:,1) = linspace(cam.position(1),cam.hitPoint(1),ST.camParams(4));
mainLineOfSight(:,2) = linspace(cam.position(2),cam.hitPoint(2),ST.camParams(4));
mainLineOfSight(:,3) = linspace(cam.position(3),cam.hitPoint(3),ST.camParams(4));
% main line of sight

LOS.main.sampling = struct;
LOS.main.sampling.coord = mainLineOfSight;
LOS.main.sampling.R = sqrt(LOS.main.sampling.coord(:,1).^2 + ...
    LOS.main.sampling.coord(:,2).^2);
LOS.main.sampling.Z = LOS.main.sampling.coord(:,3);

daw = (pi*maxAngleView/180)/ST.camParams(2);
angleView = daw:daw:ST.camParams(2)*daw;

dangle = 2*pi/ST.camParams(3);
rotAngle = 0:dangle:(ST.camParams(3)-1)*dangle;

LOS.angleView = angleView;
LOS.rotAngle = rotAngle;

sampling = cell(ST.camParams(2),ST.camParams(3));

% compute relevant quantities
Ls = zeros(ST.camParams(4),3); % (x,y,z)
Lxz = zeros(ST.camParams(4),3); % rotated main line of sight lying on the xz-plane
Lx = zeros(ST.camParams(4),3); % rotated main line of sight lying on the xz-plane

Ls(:,1) = LOS.main.sampling.coord(:,1) - LOS.main.camera.position(1);
Ls(:,2) = LOS.main.sampling.coord(:,2) - LOS.main.camera.position(2);
Ls(:,3) = LOS.main.sampling.coord(:,3) - LOS.main.camera.position(3);

anglxy= atan(Ls(end,2)/Ls(end,1));
Lxz(:,1) = cos(anglxy)*Ls(:,1) + sin(anglxy)*Ls(:,2);
Lxz(:,2) = -sin(anglxy)*Ls(:,1) + cos(anglxy)*Ls(:,2);
Lxz(:,3) = Ls(:,3);

anglxz = atan(Lxz(end,3)/Lxz(end,1));
Lx(:,1) = cos(anglxz)*Lxz(:,1) + sin(anglxz)*Lxz(:,3);
Lx(:,2) = Lxz(:,2);
Lx(:,3) = -sin(anglxz)*Lxz(:,1) + cos(anglxz)*Lxz(:,3);

for ii=1:ST.camParams(2)
    Ltmp = zeros(ST.camParams(4),3);
    Ltmp(:,1) = Lx(:,1)*cos(angleView(ii)) + Lx(:,3)*sin(angleView(ii));
    Ltmp(:,2) = Lx(:,2);
    Ltmp(:,3) = -Lx(:,1)*sin(angleView(ii)) + Lx(:,3)*cos(angleView(ii));
    
    for jj=1:ST.camParams(3)
        sampling{ii,jj} = struct;
        sampling{ii,jj}.coord = zeros(2*ST.camParams(4),3);
        sampling{ii,jj}.rotAngle = rotAngle(jj);
        sampling{ii,jj}.angleView = angleView(ii);
        
        L = zeros(ST.camParams(4),3);
        L(:,1) = Ltmp(:,1);
        L(:,2) = cos(rotAngle(jj))*Ltmp(:,2) + sin(rotAngle(jj))*Ltmp(:,3);
        L(:,3) = -sin(rotAngle(jj))*Ltmp(:,2) + cos(rotAngle(jj))*Ltmp(:,3);
        
        Lxztmp = zeros(ST.camParams(4),3);
        Lxztmp(:,1) = cos(anglxz)*L(:,1) - sin(anglxz)*L(:,3);
        Lxztmp(:,2) = L(:,2);
        Lxztmp(:,3) = sin(anglxz)*L(:,1) + cos(anglxz)*L(:,3);
        
        Lstmp = zeros(ST.camParams(4),3);
        Lstmp(:,1) = cos(anglxy)*Lxztmp(:,1) - sin(anglxy)*Lxztmp(:,2);
        Lstmp(:,2) = sin(anglxy)*Lxztmp(:,1) + cos(anglxy)*Lxztmp(:,2);
        Lstmp(:,3) = Lxztmp(:,3);
        
        sampling{ii,jj}.coord(:,1) = linspace(Lstmp(1,1),...
            2*Lstmp(end,1),2*ST.camParams(4));
        sampling{ii,jj}.coord(:,2) = linspace(Lstmp(1,2),...
            2*Lstmp(end,2),2*ST.camParams(4));
        sampling{ii,jj}.coord(:,3) = linspace(Lstmp(1,3),...
            2*Lstmp(end,3),2*ST.camParams(4));
        
        sampling{ii,jj}.coord(:,1) = sampling{ii,jj}.coord(:,1) + ...
            LOS.main.camera.position(1);
        sampling{ii,jj}.coord(:,2) = sampling{ii,jj}.coord(:,2) + ...
            LOS.main.camera.position(2);
        sampling{ii,jj}.coord(:,3) = sampling{ii,jj}.coord(:,3) + ...
            LOS.main.camera.position(3);
        
        sampling{ii,jj}.R = sqrt(sampling{ii,jj}.coord(:,1).^2 + ...
            sampling{ii,jj}.coord(:,2).^2);
        sampling{ii,jj}.Z = sampling{ii,jj}.coord(:,3);
        
    end
    clear Ltmp L Lxztmp Lstmp
end

LOS.numAngleViews = ST.camParams(2);
LOS.numRotAngles = ST.camParams(3);

% % Determine whether the sampling points are in the core region
inner_pol = [0,ST.IB.inner_pol,2*pi];
inner_r = [ST.IB.inner_r(1),ST.IB.inner_r,ST.IB.inner_r(end)];

% main line of sight
smpPt_r = sqrt((LOS.main.sampling.R - ST.geom.R0).^2 + ...
    LOS.main.sampling.Z.^2);
smpPt_pol = atan2(LOS.main.sampling.Z,LOS.main.sampling.R - ST.geom.R0);

I = smpPt_pol < 0;
smpPt_pol(I) = smpPt_pol(I) + 2*pi;
clear I

tmp_inner_r = interp1(inner_pol,inner_r,smpPt_pol,'linear');
aux = smpPt_r < tmp_inner_r;
LOS.main.sampling.R(aux) = [];
LOS.main.sampling.Z(aux) = [];
LOS.main.sampling.coord(aux,:) = [];

clear smpPt_pol smpPt_r tmp_inner_r aux
% main line of sight

for ii=1:LOS.numAngleViews
    for jj=1:LOS.numRotAngles
        smpPt_r = sqrt((sampling{ii,jj}.R - ST.geom.R0).^2 + ...
            sampling{ii,jj}.Z.^2);
        smpPt_pol = atan2(sampling{ii,jj}.Z,sampling{ii,jj}.R - ST.geom.R0);
        
        I = smpPt_pol < 0;
        smpPt_pol(I) = smpPt_pol(I) + 2*pi;
        clear I
        
        tmp_inner_r = interp1(inner_pol,inner_r,smpPt_pol,'linear');
        aux = smpPt_r < tmp_inner_r;
        sampling{ii,jj}.R(aux) = [];
        sampling{ii,jj}.Z(aux) = [];
        sampling{ii,jj}.coord(aux,:) = [];
        
        clear smpPt_pol smpPt_r tmp_inner_r aux
    end
end
% % Determine whether the sampling points are in the core region

% % % % % % % % % % % % MODIFY THIS PARAMETERS % % % % % % % % % % % % 
% Uncomment the plot below to choose the alternative coordinate systems
% Temporary figure
% figure
% plot3(R(2:end-1),Z(2:end-1),2:1:length(R)-1,'k-o','LineWidth',2)
% axis equal tight

switch ST.CASE
    case 'JET'
        Z_threshold = 0.75*min(Z); % JET
        R_threshold = mean(R); % JET
        
        RADD = [mean(R), 0.8*mean(R)]; % JET
        ZADD = [1.1*min(Z),min(Z)]; % JET
        
        R2 = 1.2*mean(R); % JET and ASDEX
        Z2 = min(Z); % JET and ASDEX
        
    case 'ASDEX'
        Z_threshold = 0.65*min(Z); % ASDEX
        R_threshold = mean(R); % ASDEX
        
        RADD = [0.9*mean(R), 0.5*mean(R)]; % ASDEX
        ZADD = [1.1*min(Z),0.88*min(Z)]; % ASDEX
        
        R2 = 1.2*mean(R); % JET and ASDEX
        Z2 = min(Z); % JET and ASDEX
        
    case 'WEST'
        Z_threshold = 0.7*min(Z); % WEST
        R_threshold = mean(R); % WEST
        
        RADD = [1.05*mean(R), 0.7*mean(R)]; % WEST
        ZADD = [1.2*min(Z),1.8*min(Z)]; % WEST
        
        R2 = 1.2*mean(R); % WEST
        Z2 = 0.95*min(Z); % WEST
        
    case 'TCV'
        Z_threshold = 2*min(Z); % TCV
        R_threshold = mean(R); % TCV
        
        RADD = []; % TCV
        ZADD = []; % TCV
        
        R2 = 1.2*mean(R); % TCV
        Z2 = 0.95*min(Z); % TCV
        
    otherwise
        error('Introduce a valid CASE');
end

% % % % % % % % % % % % MODIFY THIS PARAMETERS % % % % % % % % % % % % 

indZoneDivertor = find(Z < Z_threshold);
indZoneMainPlasma = find(Z >= Z_threshold);
indBreak = [find(diff(indZoneMainPlasma) ~= 1,1,'first'),...
    find(diff(indZoneMainPlasma) ~= 1,1,'first')+1];

% % % % % % % % % % % % MODIFY THIS PARAMETERS % % % % % % % % % % % % 
indZone = cell(1,length(RADD));
switch ST.CASE
    case 'JET'
        indZone{1} = 178:1:203; % JET
        indZone{2} = 159:1:203; % JET
        
    case 'ASDEX'
        indZone{1} = 37:1:45; % ASDEX
        indZone{2} = min(indZoneDivertor):1:36; % ASDEX
        
    case 'WEST'
        indZone{1} = 70:1:71; % WEST
        indZone{2} = min(indZoneDivertor)-1:1:70; % WEST
        
    case 'TCV'
        
    otherwise
        error('Introduce a valid CASE');
end
% % % % % % % % % % % % MODIFY THIS PARAMETERS % % % % % % % % % % % % 

% % % % % % % % % VARIABLES OF ALTERNATIVE COORD SYSTEMS % % % % % % %

wall_r_zone2 = sqrt( (R - R2).^2 + (Z - Z2).^2 );
wall_pol_zone2 = atan2(Z - Z2,R - R2);
I = wall_pol_zone2 < 0;
wall_pol_zone2(I) = wall_pol_zone2(I) + 2*pi;
clear I

wall_r_zone = zeros(length(R),length(RADD));
wall_pol_zone = zeros(length(R),length(RADD));

for ll=1:length(RADD)
    wall_r_zone(:,ll) = sqrt( (R - RADD(ll)).^2 + (Z - ZADD(ll)).^2 );
    wall_pol_zone(:,ll) = atan2(Z - ZADD(ll),R - RADD(ll));
    I = wall_pol_zone(:,ll) < 0;
    wall_pol_zone(I,ll) = wall_pol_zone(I,ll) + 2*pi;
    clear I
end
% % % % % % % % % VARIABLES OF ALTERNATIVE COORD SYSTEMS % % % % % % %

T=figure;
hold on
plot(R,Z,'k','LineWidth',2)
plot(ST.geom.R0,ST.geom.Z0,'k.','MarkerSize',12)
plot(R2,Z2,'g.','MarkerSize',12)
for ll=1:length(RADD)
    plot(RADD(ll),ZADD(ll),'g.','MarkerSize',12)
end
plot(linspace(min(R),max(R),10),Z_threshold*ones(1,10),'b')
hold off
axis equal tight
box on

for ii=1:LOS.numAngleViews
    disp(['Angle of view: ' num2str(ii)]);
    for jj=1:LOS.numRotAngles
        IND = [];
        INDADD = [];
        for pp=1:length(sampling{ii,jj}.R) % sampling points
            
            distToWall = sqrt((R - sampling{ii,jj}.R(pp)).^2 +...
                (Z - sampling{ii,jj}.Z(pp)).^2);
            [~,S] = sort(distToWall(2:end-1));
            S = S + 1;
            
            if isempty(find(indZoneDivertor == S(1),1,'first')) % NOT IN THE DIVERTOR
                
                smpPt_r = sqrt((sampling{ii,jj}.R(pp) - ST.geom.R0).^2 + ...
                    (sampling{ii,jj}.Z(pp) - ST.geom.Z0).^2);
                smpPt_pol = atan2(sampling{ii,jj}.Z(pp) - ST.geom.Z0,sampling{ii,jj}.R(pp) - ST.geom.R0);
                if smpPt_pol < 0
                    smpPt_pol = smpPt_pol + 2*pi;
                end
                
                [~,Spol] = sort(abs(smpPt_pol - wall_pol(2:end-1)));
                Spol = Spol + 1;
                
                closestWallPt = intersect(Spol,S,'stable');
                
                pivot = closestWallPt(1);
                refPolDist = smpPt_pol - wall_pol(pivot);

                try
                    if pivot <= indBreak(1)
                        upperLim = indBreak(1);
                        lowerLim = (min(indZoneMainPlasma)+1);
                    else
                        upperLim = (max(indZoneMainPlasma)-1);
                        lowerLim = indBreak(2);
                    end
                catch
                    upperLim = max(indZoneMainPlasma)-1;
                    lowerLim = min(indZoneMainPlasma)+1;
                end
                
                skipInterp = false;
                
                if (refPolDist == 0)
                    skipInterp = true;
                    wall_r_extr = wall_r(pivot);
                elseif (refPolDist > 0)
                    
                    for ll=pivot:upperLim
                        polDist = smpPt_pol - wall_pol(ll);
                        if (wall_pol(ll+1) - wall_pol(ll)) > polDist
                            index1_increasing = ll;
                            index2_increasing = ll+1;
                            break
                        end
                    end
                    
                    for ll=pivot:-1:lowerLim
                        polDist = smpPt_pol - wall_pol(ll);
                        if (wall_pol(ll-1) - wall_pol(ll)) > polDist
                            index1_decreasing = ll-1;
                            index2_decreasing = ll;
                            break
                        end
                    end
                    
                    var_bool1 = exist('index1_increasing','var') && exist('index2_increasing','var');
                    var_bool2 = exist('index1_decreasing','var') && exist('index2_decreasing','var');
                    
                    if var_bool1 && var_bool2
                        % distance to the line defined by the wall segments
                        m = diff(Z(index1_increasing:index2_increasing))/diff(R(index1_increasing:index2_increasing));
                        b = Z(index1_increasing) - m*R(index1_increasing);
                        dist1 = abs(m*sampling{ii,jj}.R(pp) - sampling{ii,jj}.Z(pp) +  b)/sqrt(m^2 + 1);
                        
                        m = diff(Z(index1_decreasing:index2_decreasing))/diff(R(index1_decreasing:index2_decreasing));
                        b = Z(index1_decreasing) - m*R(index1_decreasing);
                        dist2 = abs(m*sampling{ii,jj}.R(pp) - sampling{ii,jj}.Z(pp) +  b)/sqrt(m^2 + 1);
                        
                        if dist1 < dist2
                            INDEX1 = index1_increasing;
                            INDEX2 = index2_increasing;
                        else
                            INDEX1 = index1_decreasing;
                            INDEX2 = index2_decreasing;
                        end
                    elseif var_bool1
                        INDEX1 = index1_increasing;
                        INDEX2 = index2_increasing;
                    elseif var_bool2
                        INDEX1 = index1_decreasing;
                        INDEX2 = index2_decreasing;
                    end
                    
                    if ~var_bool2 && ~var_bool1
                        disp('EXCEPTION 1: FATAL')
                        disp(['At ii=' num2str(ii) ' jj=' num2str(jj) ' pp=' num2str(pp)]);
                        return
                    end
                    
                elseif (refPolDist < 0)
                    
                    for ll=pivot:-1:lowerLim
                        polDist = smpPt_pol - wall_pol(ll);
                        if (wall_pol(ll-1) - wall_pol(ll)) < polDist
                            index1_decreasing = ll-1;
                            index2_decreasing = ll;
                            break
                        end
                    end
                    
                    for ll=pivot:upperLim
                        polDist = smpPt_pol - wall_pol(ll);
                        if (wall_pol(ll+1) - wall_pol(ll)) < polDist
                            index1_increasing = ll;
                            index2_increasing = ll+1;
                            break
                        end
                    end
                    
                    var_bool1 = exist('index1_increasing','var') && exist('index2_increasing','var');
                    var_bool2 = exist('index1_decreasing','var') && exist('index2_decreasing','var');
                    
                    if var_bool1 && var_bool2
                        % distance to the line defined by the wall segments
                        m = diff(Z(index1_increasing:index2_increasing))/diff(R(index1_increasing:index2_increasing));
                        b = Z(index1_increasing) - m*R(index1_increasing);
                        dist1 = abs(m*sampling{ii,jj}.R(pp) - sampling{ii,jj}.Z(pp) +  b)/sqrt(m^2 + 1);
                        
                        m = diff(Z(index1_decreasing:index2_decreasing))/diff(R(index1_decreasing:index2_decreasing));
                        b = Z(index1_decreasing) - m*R(index1_decreasing);
                        dist2 = abs(m*sampling{ii,jj}.R(pp) - sampling{ii,jj}.Z(pp) +  b)/sqrt(m^2 + 1);
                        
                        if dist1 < dist2
                            INDEX1 = index1_increasing;
                            INDEX2 = index2_increasing;
                        else
                            INDEX1 = index1_decreasing;
                            INDEX2 = index2_decreasing;
                        end
                    elseif var_bool1
                        INDEX1 = index1_increasing;
                        INDEX2 = index2_increasing;
                    elseif var_bool2
                        INDEX1 = index1_decreasing;
                        INDEX2 = index2_decreasing;
                    end
                    
                    if ~var_bool2 && ~var_bool1
                        disp('EXCEPTION!')
                        disp(['At ii=' num2str(ii) ' jj=' num2str(jj) ' pp=' num2str(pp)]);
                    end
                    
                end
                
                clear index1_increasing index2_increasing
                clear index1_decreasing index2_decreasing
                
                if ~skipInterp
                    wall_r_extr = ...
                        interp1(wall_pol(INDEX1:INDEX2),wall_r(INDEX1:INDEX2),...
                        smpPt_pol,'linear');
                end
                
                if smpPt_r > wall_r_extr
                    IND = pp:1:length(sampling{ii,jj}.R);
                    break
                end
                
            elseif ~isempty(find(indZoneDivertor == S(1),1,'first')) % IN THE DIVERTOR
                
                smpPt_r = sqrt((sampling{ii,jj}.R(pp) - R2).^2 + ...
                    (sampling{ii,jj}.Z(pp) - Z2).^2);
                smpPt_pol = atan2(sampling{ii,jj}.Z(pp) - Z2,sampling{ii,jj}.R(pp) - R2);
                if smpPt_pol < 0
                    smpPt_pol = smpPt_pol + 2*pi;
                end
                
                for ll=S(1):max(indZoneDivertor)
                    if length(unique(wall_pol_zone2(ll:ll+1))) ~= 1 && length(unique(wall_r_zone2(ll:ll+1))) ~= 1
                        wall_r_extr = ...
                            interp1(wall_pol_zone2(ll:ll+1),wall_r_zone2(ll:ll+1),...
                            smpPt_pol,'linear');
                    else
                        continue
                    end
                    if ~isnan(wall_r_extr)
                        break
                    end
                end
                
                if isnan(wall_r_extr)
                    for ll=S(1):-1:min(indZoneDivertor)
                        if length(unique(wall_pol_zone2(ll-1:ll))) ~= 1 && length(unique(wall_r_zone2(ll-1:ll))) ~= 1
                            wall_r_extr = ...
                                interp1(wall_pol_zone2(ll-1:ll),wall_r_zone2(ll-1:ll),...
                                smpPt_pol,'linear');
                        else
                            continue
                        end
                        if ~isnan(wall_r_extr)
                            break
                        end
                    end
                end
                
                if smpPt_r < wall_r_extr
                    IND = pp:1:length(sampling{ii,jj}.R);
                    break
                end
                
                if isempty(IND)
                    for kk=1:length(RADD)
                        if ~isempty(find(indZone{kk} == S(1),1))
                            smpPt_r = sqrt((sampling{ii,jj}.R(pp) - RADD(kk)).^2 + ...
                                (sampling{ii,jj}.Z(pp) - ZADD(kk)).^2);
                            smpPt_pol = atan2(sampling{ii,jj}.Z(pp) - ZADD(kk),sampling{ii,jj}.R(pp) - RADD(kk));
                            if smpPt_pol < 0
                                smpPt_pol = smpPt_pol + 2*pi;
                            end
                            
                            for ll=S(1):max(indZone{kk})
                                if length(unique(wall_pol_zone(ll:ll+1,kk))) ~= 1 && length(unique(wall_r_zone(ll:ll+1,kk))) ~= 1
                                    wall_r_extr = ...
                                        interp1(wall_pol_zone(ll:ll+1,kk),wall_r_zone(ll:ll+1,kk),...
                                        smpPt_pol,'linear');
                                else
                                    continue
                                end
                                if ~isnan(wall_r_extr)
                                    break
                                end
                            end
                            
                            if isnan(wall_r_extr)
                                for ll=S(1):-1:min(indZone{kk})
                                    if length(unique(wall_pol_zone(ll-1:ll,kk))) ~= 1 && length(unique(wall_r_zone(ll-1:ll,kk))) ~= 1
                                        wall_r_extr = ...
                                            interp1(wall_pol_zone(ll-1:ll,kk),wall_r_zone(ll-1:ll,kk),...
                                            smpPt_pol,'linear');
                                    else
                                        continue
                                    end
                                    if ~isnan(wall_r_extr)
                                        break
                                    end
                                end
                            end
                            
                            if smpPt_r < wall_r_extr
                                INDADD = pp:1:length(sampling{ii,jj}.R);
                                break
                            end
                        end % if ~isempty(find(indZon{kk} == S(1),1))
                    end
                end % if isempty(IND)
                
            end % if isempty(find(indZoneDivertor == S(1),1)) % not in the divertor
            
            if isnan(wall_r_extr)
                INDADD = pp:1:length(sampling{ii,jj}.R);
            end
            
            if ~isempty(INDADD)
                break
            end
            
        end% sampling points
        
        sampling{ii,jj}.R(IND) = [];
        sampling{ii,jj}.Z(IND) = [];
        sampling{ii,jj}.coord(IND,:) = [];
        
        sampling{ii,jj}.R(INDADD) = [];
        sampling{ii,jj}.Z(INDADD) = [];
        sampling{ii,jj}.coord(INDADD,:) = [];
        clear smpPt_r smpPt_pol IND INDADD
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Determine whether the sampling points are within the simulation domain


figure(ST.IB.figh)
hold on
for ii=1:LOS.numAngleViews
    for jj=1:LOS.numRotAngles
        plot3(sampling{ii,jj}.coord(:,1),...
            sampling{ii,jj}.coord(:,2),...
            sampling{ii,jj}.coord(:,3),...
            'r.','MarkerSize',3)
    end
end
view([0 90])
xlabel('$x$-axis','FontSize',14,'Interpreter','latex')
ylabel('$y$-axis','FontSize',14,'Interpreter','latex')
zlabel('$z$-axis','FontSize',14,'Interpreter','latex')
hold off

figure(T)
hold on
for ii=1:LOS.numAngleViews
    for jj=1:LOS.numRotAngles
        plot(sampling{ii,jj}.R,sampling{ii,jj}.Z,...
            'r.','MarkerSize',4)
    end
end
hold off
axis equal tight
box on

LOS.sampling = sampling;

disp('Lines of sight calculated.')
end

function C = fastCamera(ST,paramToPlot)
disp('Initialising fast camera...')

CAM = struct;
% camera parameters
CAM.f = 40;% focal length
CAM.Rdiag = 5.0; % radius of the camera's diaphragm
% camera parameters

los = calculateLinesOfSight(ST,ST.camPos,ST.hitPointPos,ST.camParams(1));

% detector parameters
numLinesOfSight = los.numAngleViews*los.numRotAngles + 1;
CAM.vals = zeros(numLinesOfSight,3);% (x,y,value) per line of sight

nx = 1000; % 1 mega-pixel
ny = 1000; % 1 mega-pixel
% detector parameters

disp('Calculating scattered interpolant...')
if strcmp(paramToPlot,'density')
    R = reshape(ST.geom.knots.R,[numel(ST.geom.knots.R) 1]);
    Z = reshape(ST.geom.knots.Z,[numel(ST.geom.knots.Z) 1]);
    V = reshape(ST.data.(paramToPlot),[numel(ST.data.(paramToPlot)) 1]);
    F = scatteredInterpolant(R,Z,V);
else
    R = reshape(ST.geom.knots.R_tri,[numel(ST.geom.knots.R_tri) 1]);
    Z = reshape(ST.geom.knots.Z_tri,[numel(ST.geom.knots.Z_tri) 1]);
    V = reshape(ST.data.(paramToPlot),[numel(ST.data.(paramToPlot)) 1]);
    F = scatteredInterpolant(R,Z,V);
end
disp('Scattered interpolant CALCULATED.')

% value at the centre of the detector
dl = sqrt( min(diff(los.main.sampling.coord(:,1)).^2) + ...
    min(diff(los.main.sampling.coord(:,2)).^2) + ...
    min(diff(los.main.sampling.coord(:,3)).^2) );

CAM.vals(1,1) = 0;
CAM.vals(1,2) = 0;
tmp = F(los.main.sampling.R,los.main.sampling.Z);
r_tmp = sqrt( (los.main.sampling.coord(:,1) - los.main.camera.position(1)).^2 + ...
    (los.main.sampling.coord(:,2) - los.main.camera.position(2)).^2 + ...
    (los.main.sampling.coord(:,3) - los.main.camera.position(3)).^2 );
theta_tmp = atan(CAM.Rdiag./r_tmp);
if ST.optics
    CAM.val(1,3) = sum( 0.5*(1 - cos(theta_tmp(2:end))).*tmp(2:end))*dl;
else
    CAM.val(1,3) = sum( tmp )*dl;
end
clear tmp dl r_tmp theta_tmp

for ii=1:los.numAngleViews
    for jj=1:los.numRotAngles
        id = (ii-1)*los.numRotAngles + jj + 1;
        
        dl = sqrt( min(diff(los.sampling{ii,jj}.coord(:,1)).^2) + ...
            min(diff(los.sampling{ii,jj}.coord(:,2)).^2) + ...
            min(diff(los.sampling{ii,jj}.coord(:,3)).^2) );
        
        tmp = F(los.sampling{ii,jj}.R,los.sampling{ii,jj}.Z);
        
        h = -CAM.f*tan(los.sampling{ii,jj}.angleView);
        rotAngleSensor = los.sampling{ii,jj}.rotAngle - pi/2;
        
        CAM.vals(id,1) = h*sin(rotAngleSensor);
        CAM.vals(id,2) = h*cos(rotAngleSensor);
        
        %         Here we calculate the value of los.sampling{ii,jj}.vals captured
        %         by the camera. It is assumed that los.sampling{ii,jj} is an
        %         isotropic source.
        r_tmp = sqrt( (los.sampling{ii,jj}.coord(:,1) - los.main.camera.position(1)).^2 + ...
            (los.sampling{ii,jj}.coord(:,2) - los.main.camera.position(2)).^2 + ...
            (los.sampling{ii,jj}.coord(:,3) - los.main.camera.position(3)).^2 );
        theta_tmp = atan(CAM.Rdiag./r_tmp);
        if ST.optics
            CAM.vals(id,3) = sum( 0.5*(1 - cos(theta_tmp(2:end))).*tmp(2:end) )*dl;
        else
            CAM.vals(id,3) = sum( tmp )*dl;
        end
        clear r_tmp theta_tmp tmp dl rotAngleSensor
    end
end

PIC = scatteredInterpolant(CAM.vals(:,1),CAM.vals(:,2),CAM.vals(:,3));

xmin = min(CAM.vals(:,1));
xmax = max(CAM.vals(:,1));
ymin = min(CAM.vals(:,2));
ymax = max(CAM.vals(:,2));

tmpx = linspace(xmin,xmax,nx);
tmpy = linspace(ymin,ymax,ny);

[X,Y] = meshgrid(tmpx,tmpy);
S = PIC(X,Y);

sh=figure;
imagesc(tmpx,tmpy,fliplr(S'))
axis xy
axis([xmin xmax ymin ymax])
axis square
colormap(jet(512))
grid on

xlabel('x-axis','interpreter','latex','FontSize',14);
ylabel('y-axis','interpreter','latex','FontSize',14);

bound_xaxis = CAM.vals(end-los.numRotAngles+1:end,1);
bound_yaxis = CAM.vals(end-los.numRotAngles+1:end,2);

figure(sh)
hold on
plot(bound_xaxis,bound_yaxis,'k-','LineWidth',2)
hold off


disp('Fast camera: DONE!')
end