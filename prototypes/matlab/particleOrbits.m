function ST = particleOrbits(pathToBField,ND,res,timeStepParams,tracerParams,xo,vo_params,opt)
% Here ro and vo are the initial position and velocity of the particle.
% The components of the initial velocity must be entered as fractions of
% the speed of light in vacuum.
%
% EXAMPLE:
% The example below traces an alpha-particle in an analytical magnetic
% field. The parameters of the magnetic field are hard-coded in the
% functions 'analyticalB' and 'DiegosInvariant'.
% ST = particleOrbits('some_VMEC_file.dat','2D',[150,150],[1E4,1E-2,10],[2,7.2938E3],[6,0,1],[0.03,80]);
% ST = particleOrbits('','2D',[],[1E4,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.04,80]);
% ST = particleOrbits('CHEBYSHEV.dat','2D',[50,50],[1E4,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.04,80]);

narginchk(7,8);

close all

ST = struct;
% Script parameters
if strcmp(pathToBField,'')
    ST.analytical = true; % true = analytical B field; false = xpander fields
    ST.pathToBField = pathToBField; % Path to xpander fields
else
    ST.analytical = false; % true = analytical B field; false = xpander fields
    ST.pathToBField = pathToBField; % Path to xpander fields
    ST.res = res;
end

ST.ND = ND; % Dimensionality of xpander fields, 2D or 3D.
if nargin == 7
    ST.opt = true; % true = plot xpander fields, false = don't plot xpander fields
else
    ST.opt = opt;
end

% Plasma parameters and physical constants, all in SI units
ST.params = struct;
ST.params.kB = 1.38E-23; % Boltzmann constant
ST.params.mu0 = (4E-7)*pi; % Magnetic permeability
ST.params.ep = 8.854E-12;% Electric permitivity
ST.params.c = 2.9979E8; % Speed of light
ST.params.qe = 1.602176E-19; % Electron charge
ST.params.me = 9.109382E-31; % Electron mass

% Electric and magneticfield
ST.E = [0,0,0];
if ST.analytical
    ST.B = analyticalB([1,1,1],'initialize');
    ST.Bo = ST.B.Bo;
else
    ST.B = loadMagneticField(ST);
    ST.Bo = ST.B.Bo;
end

% Initial position and velocity of tracer, in SI units
ST.params.Xo = xo; % Initial position
ST.params.vo_params = vo_params; % vo_params = [velocity magnitude, pitch angle]
[ST.params.vo, ST.params.vpar, ST.params.vperp] = initializeVelocity(ST);

% ST.params.vo = 1E7*[-0.002449618114251   0.208217539711322   1.180942065102119];
% ST.params.vpar = -2.082319487310763E6;
% ST.params.vperp = -1.180942065102119E7;

% Particle's parameters
ST.params.q = tracerParams(1)*ST.params.qe; %alpha-particle
ST.params.m = tracerParams(2)*ST.params.me; % alpha-particle
ST.params.wc = sqrt(1 - sum(ST.params.vo.^2)/ST.params.c^2)*abs(ST.params.q)*ST.Bo/ST.params.m;
ST.params.rL = ST.params.vperp/ST.params.wc; % Larmor radius of particle

disp(['Cyclotron frequency ' num2str(ST.params.wc) ' Hz'])
disp(['Larmor radius ' num2str(ST.params.rL) ' m'])

% Normalisation parameters
ST.norm.q = abs(ST.params.q);
ST.norm.m = ST.params.m;
ST.norm.wc = ST.norm.q*ST.Bo/ST.norm.m;
ST.norm.l = ST.params.c/ST.norm.wc;
% Normalisation parameters

% Numerical parameters
ST.params.numIt = timeStepParams(1);
ST.params.dt = timeStepParams(2)*(2*pi/ST.params.wc);
ST.params.cadence = timeStepParams(3);
ST.params.inds = 1:ST.params.cadence:ST.params.numIt;

ST.time = zeros(ST.params.numIt,1);
for ii=1:ST.params.numIt
    ST.time(ii) = (ii-1)*ST.params.dt;
end

ST.PP = particlePusherLeapfrogMod(ST);
PoincarePlots(ST.PP);
% DiegosInvariants(ST);

% ST.PP = particlePusherMatlab(ST);
% PoincarePlots(ST.PP.X);
% DiegosInvariants(ST);

end

function [vo,vpar,vperp] = initializeVelocity(ST)
% ST.params.vo_params = [vo, pitch_angle]
v = ST.params.vo_params(1)*ST.params.c; % magnitude of initial velocity
pitchAngle = pi*ST.params.vo_params(2)/180;

vpar = v*cos(pitchAngle);
vperp = v*sin(pitchAngle);

[b,a] = unitVectors(ST);
vo = vpar*b + vperp*a;

disp(['Parallel velocity: ' num2str(vpar/ST.params.c)])
disp(['Perpendicular velocity: ' num2str(vperp/ST.params.c)])

end

function P = PoincarePlots(PP)
X = PP.X;

R = sqrt(X(:,1).^2 + X(:,2).^2);
Z = X(:,3);

if isfield(PP,'R')
    Rgc = sqrt(PP.R(:,1).^2 + PP.R(:,2).^2);
    Zgc = PP.R(:,3);
end

zeta = atan2(X(:,2),X(:,1));
locs = find(abs(diff(zeta)) > 6);

figure
plot(R,X(:,3),'k',R(locs),Z(locs),'r.','MarkerSize',15)
if isfield(PP,'R')
    hold on
    plot(Rgc,Zgc,'g','LineWidth',2)
    hold off
end
axis equal
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
title('Poincare plot','Interpreter','latex','FontSize',16)

end

function [b,a] = unitVectors(ST)
% initial condition of an electron drifting parallel to the local magnetic
% field.

% tmp = weightedInterpolation(ST,ST.params.Xo);

if ST.analytical
    B = analyticalB(ST.params.Xo);
    b = B/sqrt(sum(B.^2));
else
    B = interpMagField(ST,ST.params.Xo);
    b = B/sqrt(sum(B.^2));
end

az = -( b(1) + b(2) )/b(3);
if isfinite(az)
    a = [1,1,az]/( 2 + az^2 );
    a = a/sqrt(sum(a.^2));
else
    ay = -b(1)/b(2);
    a = [1,ay,0];
    a = a/sqrt(sum(a.^2));
end

end

function B = loadMagneticField(ST)
% All quantities in SI units
B = struct;

if isempty(ST.res) % magnetic field at scattered positions
    
    data = load(ST.pathToBField);
    if strcmp(ST.ND,'2D')
        B.R = data(:,1);
        B.phi = data(:,2);
        B.Z = data(:,3);
        
        B.BR = data(:,5);
        B.Bphi = data(:,6);
        B.BZ = data(:,7);
    elseif strcmp(ST.ND,'3D')
        B.R = data(:,1);
        B.phi = data(:,2);
        B.Z = data(:,3);
        
        B.BR = data(:,5);
        B.Bphi = data(:,6);
        B.BZ = data(:,7);
    end
    
    B.B = sqrt(B.BR.^2 + B.BZ.^2 + B.Bphi.^2);
    
    B.Bo = mean(mean(mean(B.B)));
    
    disp('NOTE: Magnetic field in non-uniform grid!')
    
    figure
    subplot(1,3,1)
    scatter3(B.R,B.Z,B.BR,'k.')
    axis square
    box on
    title('$B^R$ [T]','Interpreter','latex','FontSize',16)
    xlabel('R [m]','Interpreter','latex','FontSize',16)
    ylabel('Z [m]','Interpreter','latex','FontSize',16)
    
    subplot(1,3,2)
    scatter3(B.R,B.Z,B.Bphi,'k.')
    axis square
    box on
    title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
    xlabel('R [m]','Interpreter','latex','FontSize',16)
    ylabel('Z [m]','Interpreter','latex','FontSize',16)
    
    subplot(1,3,3)
    scatter3(B.R,B.Z,B.BZ,'k.')
    axis square
    box on
    title('$B^Z$ [T]','Interpreter','latex','FontSize',16)
    xlabel('R [m]','Interpreter','latex','FontSize',16)
    ylabel('Z [m]','Interpreter','latex','FontSize',16)
    
    B.SI = calculateChebyshevInterpolant(ST,B);
    
else % structured data for magnetic field
    
    if strcmp(ST.ND,'2D')
        data = load(ST.pathToBField);
        % X(R,Z,phi)
        B.NR = ST.res(1);
        B.NZ = ST.res(2);
        
        B.R = zeros(B.NR,B.NZ);
        B.Z = zeros(B.NR,B.NZ);
        
        B.B = zeros(B.NR,B.NZ); % magnitude
        B.BR = zeros(B.NR,B.NZ);
        B.BZ = zeros(B.NR,B.NZ);
        B.P = zeros(B.NR,B.NZ);
        
        if size(data,2) > 7
            
            for iz=1:B.NZ
                indi = (iz-1)*B.NR + 1;
                indf = iz*B.NR;
                B.R(:,iz) = data(indi:indf,1);
                B.Z(:,iz) = data(indi:indf,3);
                
                B.BR(:,iz) = data(indi:indf,4);
                B.Bphi(:,iz) = data(indi:indf,5);
                B.BZ(:,iz) = data(indi:indf,6);
                
                B.P(:,iz) = data(indi:indf,7);
            end
            
            disp('NOTE: Using a expanded magnetic field')
            
        else
            
            for iz=1:B.NZ
                indi = (iz-1)*B.NR + 1;
                indf = iz*B.NR;
                B.R(:,iz) = data(indi:indf,1);
                B.Z(:,iz) = data(indi:indf,3);
                
                B.BR(:,iz) = data(indi:indf,5);
                B.Bphi(:,iz) = data(indi:indf,6);
                B.BZ(:,iz) = data(indi:indf,7);
            end
            
            disp('NOTE: Using a VMEC magnetic field')
            
        end
        
        B.B = sqrt(B.BR.^2 + B.BZ.^2 + B.Bphi.^2);
        
        % Here the question what should be the characteristic magnetic field
        % used in the normalization and calculation of the time step.
        %     B.Bo = max(max(max(B.B)));
        B.Bo = mean(mean(mean(B.B)));
        
        sizegrid = 4;
        
        if ST.opt
            figure
            subplot(1,sizegrid,1)
            surfc(B.R,B.Z,B.BR,'LineStyle','none')
            view([0,90])
            axis equal
            box on
            colorbar
            title('$B_R$ [T]','Interpreter','latex','FontSize',16)
            xlabel('R [m]','Interpreter','latex','FontSize',16)
            ylabel('Z [m]','Interpreter','latex','FontSize',16)
            
            subplot(1,sizegrid,2)
            surfc(B.R,B.Z,B.BZ,'LineStyle','none')
            view([0,90])
            axis equal
            box on
            colorbar
            title('$B_Z$ [T]','Interpreter','latex','FontSize',16)
            xlabel('R [m]','Interpreter','latex','FontSize',16)
            ylabel('Z [m]','Interpreter','latex','FontSize',16)
            
            subplot(1,sizegrid,3)
            surfc(B.R,B.Z,B.Bphi,'LineStyle','none')
            view([0,90])
            axis equal
            box on
            colorbar
            title('$B_\phi$ [T]','Interpreter','latex','FontSize',16)
            xlabel('R [m]','Interpreter','latex','FontSize',16)
            ylabel('Z [m]','Interpreter','latex','FontSize',16)
            
            if size(data,2) > 7
                subplot(1,sizegrid,4)
                surfc(B.R,B.Z,B.P,'LineStyle','none')
                view([0,90])
                axis equal
                box on
                colorbar
                title('$P(R,Z)$','Interpreter','latex','FontSize',16)
                xlabel('R [m]','Interpreter','latex','FontSize',16)
                ylabel('Z [m]','Interpreter','latex','FontSize',16)
            else
                subplot(1,sizegrid,4)
                surfc(B.R,B.Z,B.B,'LineStyle','none')
                view([0,90])
                axis equal
                box on
                colorbar
                title('$P(R,Z)$','Interpreter','latex','FontSize',16)
                xlabel('R [m]','Interpreter','latex','FontSize',16)
                ylabel('Z [m]','Interpreter','latex','FontSize',16)
            end
            
            colormap(jet)
            
        end
        
    elseif strcmp(ST.ND,'3D')
        
        data = load(ST.pathToBField);
        % X(R,Z,phi)
        B.NR = ST.res(1);
        B.Nphi = ST.res(2);
        B.NZ = ST.res(3);
        
        B.R = zeros(B.NR,B.NZ,B.Nphi);
        B.Z = zeros(B.NR,B.NZ,B.Nphi);
        B.phi = zeros(B.NR,B.NZ,B.Nphi);
        
        phi = linspace(0,2*pi,B.Nphi);
        
        B.B = zeros(B.NR,B.NZ,B.Nphi); % magnitude
        B.BR = zeros(B.NR,B.NZ,B.Nphi);
        B.Bphi = zeros(B.NR,B.NZ,B.Nphi);
        B.BZ = zeros(B.NR,B.NZ,B.Nphi);
        B.P = zeros(B.NR,B.NZ,B.Nphi);
        
        for iphi = 1:B.Nphi;
            for iz=1:B.NZ
                indi = (iz-1)*B.NR + 1;
                indf = iz*B.NR;
                B.R(:,iz,iphi) = data(indi:indf,1);
                B.Z(:,iz,iphi) = data(indi:indf,3);
                
                B.BR(:,iz,iphi) = data(indi:indf,4);
                B.Bphi(:,iz,iphi) = data(indi:indf,5);
                B.BZ(:,iz,iphi) = data(indi:indf,6);
                
                B.P(:,iz,iphi) = data(indi:indf,7);
            end
            B.B(:,:,iphi) = sqrt(B.BR(:,:,iphi).^2 + B.BZ(:,:,iphi).^2 + B.Bphi(:,:,iphi).^2);
            B.phi(:,:,iphi) = phi(iphi);
        end
        
        B.Bo = mean(mean(mean(B.B)));
        
        figure
        subplot(1,3,1)
        surfc(B.R(:,:,1),B.Z(:,:,1),B.BR(:,:,1),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B_R$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        subplot(1,3,2)
        surfc(B.R(:,:,1),B.Z(:,:,1),B.BZ(:,:,1),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B_Z$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        subplot(1,3,3)
        surfc(B.R(:,:,1),B.Z(:,:,1),B.Bphi(:,:,1),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B_\phi$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        colormap(jet)
    else
        error('Use 2D or 3D');
    end
    
    B.SI = calculateChebyshevInterpolant(ST,B);
%     B.SI = calculatescatteredInterpolant(ST,B);
    
end % if isempty(ST.res)


end

function SI = calculateChebyshevInterpolant(ST,B)
% Cylindrical coordinates
% R = radius, phi = azimuthal angle, Z = z coordinate
% calculate interpolant of the field B
disp('Calculating chebfun2 interpolant...')
SI = struct;

if strcmp(ST.ND,'2D')
    Rmin = min(min(B.R));
    Rmax = max(max(B.R));
    Zmin = min(min(B.Z));
    Zmax = max(max(B.Z));
    
%     SI.BR = chebfun2(B.BR,[Rmin Rmax Zmin Zmax]);
    SI.BR = chebfun2(B.BR,[Zmin Zmax Rmin Rmax]);
    
    SI.BZ = chebfun2(B.BZ,[Zmin Zmax Rmin Rmax]);
    
    SI.Bphi = chebfun2(B.Bphi,[Zmin Zmax Rmin Rmax]);
    
    figure
    subplot(1,3,1)
    plot(SI.BR)
    subplot(1,3,2)
    plot(SI.Bphi)
    subplot(1,3,3)
    plot(SI.BZ)
    colormap(jet)
    
elseif strcmp(ST.ND,'3D')
    SI = calculatescatteredInterpolant(ST,B);
else
    error('Please, use 2D or 3D.');
end

disp('chebfun2 interpolant: done!')
end

function SI = calculatescatteredInterpolant(ST,B)
% Cylindrical coordinates
% R = radius, phi = azimuthal angle, Z = z coordinate
% calculate interpolant of the field B
disp('Calculating scattered interpolant...')
SI = struct;


if strcmp(ST.ND,'2D')
    R = reshape(B.R,[numel(B.R) 1]);
    Z = reshape(B.Z,[numel(B.Z) 1]);
    
    DATA = reshape(B.BR,[numel(B.BR) 1]);
    SI.BR = scatteredInterpolant(R,Z,DATA);
    clear DATA
    
    DATA = reshape(B.BZ,[numel(B.BZ) 1]);
    SI.BZ = scatteredInterpolant(R,Z,DATA);
    clear DATA
    
    DATA = reshape(B.Bphi,[numel(B.Bphi) 1]);
    SI.Bphi = scatteredInterpolant(R,Z,DATA);
    clear DATA
    
elseif strcmp(ST.ND,'3D')
    R = reshape(B.R,[numel(B.R) 1]);
    Z = reshape(B.Z,[numel(B.Z) 1]);
    phi = reshape(B.phi,[numel(B.phi) 1]);
    
%     [n, bin] = histc(R, unique(R));
%     mR = find(n > 1);
%     indR    = find(ismember(bin, mR));
%     
%     [n, bin] = histc(Z, unique(Z));
%     mZ = find(n > 1);
%     indZ    = find(ismember(bin, mZ));
%     
%     [n, bin] = histc(phi, unique(phi));
%     mphi = find(n > 1);
%     indphi    = find(ismember(bin, mphi));

    DATA = reshape(B.BR,[numel(B.BR) 1]);
    SI.BR = scatteredInterpolant(R,Z,phi,DATA);
    clear DATA
    
    DATA = reshape(B.BZ,[numel(B.BZ) 1]);
    SI.BZ = scatteredInterpolant(R,Z,phi,DATA);
    clear DATA
    
    DATA = reshape(B.Bphi,[numel(B.Bphi) 1]);
    SI.Bphi = scatteredInterpolant(R,Z,phi,DATA);
    clear DATA
    
else
    error('Use 2D or 3D');
end

disp('Scattered interpolant: done!')
end

function BF = interpMagField(ST,X)
% Cylindrical coordinates
% R = radius, phi = azimuthal angle, Z = z coordinate
if strcmp(ST.ND,'2D')
    R = sqrt(X(1)^2 + X(2)^2);
    Z = X(3);
    phi = atan2(X(2),X(1));
    if phi < 0
        phi = phi + 2*pi;
    end
    
    %     disp(['(R,phi,Z) = (' num2str(R) ',' num2str(phi) ',' num2str(Z) ')'])
    
    if isa(ST.B.SI.BR,'scatteredInterpolant')
        BR = ST.B.SI.BR(R,Z);
        Bphi = ST.B.SI.Bphi(R,Z);
        BZ = ST.B.SI.BZ(R,Z);
    else
        BR = ST.B.SI.BR(Z,R);
        Bphi = ST.B.SI.Bphi(Z,R);
        BZ = ST.B.SI.BZ(Z,R);
    end
    
    Bx = BR*cos(phi) - Bphi*sin(phi);
    By = BR*sin(phi) + Bphi*cos(phi);
    Bz = BZ;
    BF = [Bx,By,Bz];
    
elseif strcmp(ST.ND,'3D')
    
else
    error('Use 2D or 3D');
end

end

function B = analyticalB(X,opt)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

narginchk(1,2);

% Parameters of the analytical magnetic field
Bo = 5;
a = 2;% 0.6;% Minor radius in meters.
Ro = 6.0; % Major radius in meters.
qa = 10; % Safety factor at the separatrix (r=a)
co = 0.5; % Extra parameter
lamb = a/co;
Bpo = (a/Ro)*(Bo/qa)*(1-co^2)/co;
% Parameters of the analytical magnetic field

if nargin == 2
    if strcmp(opt,'initialize')
        B = struct;
        B.Bo = Bo;
        B.a = a;% 0.6;% Minor radius in meters.
        B.Ro = Ro; % Major radius in meters.
        B.qa = qa; % Safety factor at the separatrix (r=a)
        B.co = co; % Extra parameter
        B.lamb = lamb;
        B.Bpo = Bpo;
    elseif strcmp(opt,'normalize')
    end
else
    % Toroidal coordinates
    % r = radius, theta = poloidal angle, zeta = toroidal angle
    r = sqrt( (sqrt(X(1)^2 + X(2)^2) - Ro)^2 + X(3)^2 );
    theta = atan2(X(3),sqrt(X(1)^2 + X(2)^2) - Ro);
    if theta < 0
        theta = theta + 2*pi;
    end
    zeta = atan2(X(1),X(2));
    if zeta < 0
        zeta = zeta + 2*pi;
    end
    % Toroidal coordinates
    
    % Poloidal magnetic field
    Bp = Bpo*(r/lamb)/( 1 + (r/lamb)^2 );
    
    eta = r/Ro;
    Br = 1/( 1 + eta*cos(theta) );
    
    Bx = Br*( Bo*cos(zeta) - Bp*sin(theta)*sin(zeta) );
    By = -Br*( Bo*sin(zeta) + Bp*sin(theta)*cos(zeta) );
    Bz = Br*Bp*cos(theta);
    
    B = [Bx,By,Bz];
end

end

function PP = particlePusherLeapfrog(ST)
% Relativistic particle pusher.
% Vay, J.-L. 2008, PoP 15, 056701.

PP = struct;

PP.method = 'Leapfrog';

X = zeros(ST.params.numIt,3); % (it,ii), ii=X,y,z
v = zeros(ST.params.numIt,3);
u = zeros(ST.params.numIt,3);

% Normalization
X(1,:) = ST.params.Xo/ST.norm.l;
v(1,:) = ST.params.vo/ST.params.c;
q = ST.params.q/ST.norm.q;
m = ST.params.m/ST.norm.m;
if ST.analytical
    B = analyticalB(X(1,:)*ST.norm.l)/ST.Bo;
else
    B = interpMagField(ST,X(1,:)*ST.norm.l)/ST.Bo;
end
dt = ST.params.dt*ST.norm.wc;
E = ST.E/(ST.Bo*ST.params.c);

u(1,:) = v(1,:)/sqrt(1 - sum(v(1,:).^2));

% initial velocity
DT = 0.5*dt;
a = q*(DT)/m;

tau = 0.5*q*DT*B/m;
up = u(1,:) + a*(E + 0.5*cross(v(1,:),B));
gammap = sqrt(1 + sum(up.^2));
sigma = gammap^2 - sum(tau.^2);
us = sum(up.*tau);
gam = sqrt(0.5)*sqrt( sigma + sqrt(sigma^2 + 4*(sum(tau.^2) + sum(us.^2))) );
t = tau/gam;
s = 1/(1+sum(t.^2));

u(1,:) = s*(up + sum(up.*t)*t + cross(up,t));
v(1,:) = u(1,:)/sqrt(1 + sum(u(1,:).^2));
% initial velocity

a = q*dt/m;

for ii=2:ST.params.numIt
    X(ii,:) = X(ii-1,:) + dt*v(ii-1,:);
    if ST.analytical
        B = analyticalB(X(ii,:)*ST.norm.l)/ST.Bo;
    else
        B = interpMagField(ST,X(ii,:)*ST.norm.l)/ST.Bo;
    end
    
    tau = 0.5*q*dt*B/m;
    up = u(ii-1,:) + a*(E + 0.5*cross(v(ii-1,:),B));
    gammap = sqrt(1 + sum(up.^2));
    sigma = gammap^2 - sum(tau.^2);
    us = sum(up.*tau); % variable 'u^*' in paper
    gam = sqrt(0.5)*sqrt( sigma + sqrt(sigma^2 + 4*(sum(tau.^2) + sum(us.^2))) );
    t = tau/gam;
    s = 1/(1+sum(t.^2)); % variable 's' in paper
    
    u(ii,:) = s*(up + sum(up.*t)*t + cross(up,t));
    v(ii,:) = u(ii,:)/sqrt(1 + sum(u(ii,:).^2));
end

X = X*ST.norm.l;

PP.X = X;
PP.v = v;

PP.EK = 1./sqrt(1-sum(PP.v.^2,2));
PP.ERR = 100*(PP.EK(1) - PP.EK)./PP.EK(1);

figure
plot3(PP.X(ST.params.inds,1),PP.X(ST.params.inds,2),PP.X(ST.params.inds,3))
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)

time = ST.time/(2*pi/ST.params.wc);

figure
plot(time(ST.params.inds), PP.ERR(ST.params.inds))
box on
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Energy conservation [\%]','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)

PP.v = PP.v*ST.params.c;
end

function PP = particlePusherLeapfrogMod(ST)
% Relativistic particle pusher.
% Vay, J.-L. 2008, PoP 15, 056701.

PP = struct;

PP.method = 'Leapfrog';

X = zeros(ST.params.numIt,3); % (it,ii), ii=X,y,z
v = zeros(ST.params.numIt,3);
u = zeros(ST.params.numIt,3);
R = zeros(ST.params.numIt,3);

% Normalization
X(1,:) = ST.params.Xo/ST.norm.l;
v(1,:) = ST.params.vo/ST.params.c;
q = ST.params.q/ST.norm.q;
m = ST.params.m/ST.norm.m;
if ST.analytical
    B = analyticalB(X(1,:)*ST.norm.l)/ST.Bo;
else
    B = interpMagField(ST,X(1,:)*ST.norm.l)/ST.Bo;
end
dt = ST.params.dt*ST.norm.wc;
E = ST.E/(ST.Bo*ST.params.c);


% initial velocity at time level t = 0
u(1,:) = v(1,:)/sqrt(1 - sum(v(1,:).^2));
v(1,:) = u(1,:)/sqrt(1 + sum(u(1,:).^2));
R(1,:) = X(1,:) + m*cross(v(1,:),B)/(q*sum(B.^2));

% initial velocity
DT = 0.5*dt;
a = q*(DT)/m;

tau = 0.5*q*DT*B/m;
up = u(1,:) + a*(E + 0.5*cross(v(1,:),B));
gammap = sqrt(1 + sum(up.^2));
sigma = gammap^2 - sum(tau.^2);
us = sum(up.*tau);
gam = sqrt(0.5)*sqrt( sigma + sqrt(sigma^2 + 4*(sum(tau.^2) + sum(us.^2))) );
t = tau/gam;
s = 1/(1+sum(t.^2));

U = s*(up + sum(up.*t)*t + cross(up,t));
V = U/sqrt(1 + sum(U.^2));
% initial velocity

a = q*dt/m;

for ii=2:ST.params.numIt
    X(ii,:) = X(ii-1,:) + dt*V;
    if ST.analytical
        B = analyticalB(X(ii,:)*ST.norm.l)/ST.Bo;
    else
        B = interpMagField(ST,X(ii,:)*ST.norm.l)/ST.Bo;
    end
    
    tau = 0.5*q*dt*B/m;
    up = U + a*(E + 0.5*cross(V,B));
    gammap = sqrt(1 + sum(up.^2));
    sigma = gammap^2 - sum(tau.^2);
    us = sum(up.*tau); % variable 'u^*' in paper
    gam = sqrt(0.5)*sqrt( sigma + sqrt(sigma^2 + 4*(sum(tau.^2) + sum(us.^2))) );
    t = tau/gam;
    s = 1/(1+sum(t.^2)); % variable 's' in paper
    
    U = s*(up + sum(up.*t)*t + cross(up,t));
    V = U/sqrt(1 + sum(U.^2));
    
    u(ii,:) = U - 0.5*a*( E + cross(V,B) );
    v(ii,:) = u(ii,:)/sqrt(1 + sum(u(ii,:).^2));
    
    
    R(ii,:) = X(ii,:) + m*cross(v(ii,:),B)/(q*sum(B.^2));
    
    disp(['Iteration ' num2str(ii)])
end

X = X*ST.norm.l;
R = R*ST.norm.l;

PP.X = X;
PP.R = R;
PP.v = v;


PP.EK = 1./sqrt(1-sum(PP.v.^2,2));
PP.ERR = 100*(PP.EK(1) - PP.EK)./PP.EK(1);

figure
plot3(PP.X(ST.params.inds,1),PP.X(ST.params.inds,2),PP.X(ST.params.inds,3),'b')
hold on
plot3(PP.R(ST.params.inds,1),PP.R(ST.params.inds,2),PP.R(ST.params.inds,3),'r')
hold off
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)

time = ST.time/(2*pi/ST.params.wc);

figure
plot(time(ST.params.inds), PP.ERR(ST.params.inds))
box on
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Energy conservation [\%]','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)

PP.v = PP.v*ST.params.c;
end

function PP = particlePusherMatlab(ST)
% Function to advance in time the particle's trajectory using Matlab
% solvers for ODEs.

PP = struct;

PP.method = 'RK-Matlab';

X = zeros(ST.params.numIt,3); % (it,ii), ii=X,y,z
v = zeros(ST.params.numIt,3);
u = zeros(ST.params.numIt,3);

% Normalization
X(1,:) = ST.params.Xo/ST.norm.l;
v(1,:) = ST.params.vo/ST.params.c;
q = ST.params.q/ST.norm.q;
m = ST.params.m/ST.norm.m;
if ST.analytical
    B = analyticalB(X(1,:)*ST.norm.l)/ST.Bo;
else
    B = interpMagField(ST,X(1,:)*ST.norm.l)/ST.Bo;
end
E = ST.E/(ST.Bo*ST.params.c);

u(1,:) = v(1,:)/sqrt(1 - sum(v(1,:).^2));
v(1,:) = u(1,:)/sqrt(1 + sum(u(1,:).^2));

options = odeset('RelTol',1E-5,'AbsTol',1E-10);%,'Stats','on','OutputFcn',@odeplot)
% options = odeset('RelTol',1E-3,'AbsTol',1E-6);%,'Stats','on','OutputFcn',@odeplot)

tspan = ST.time*ST.norm.wc;
y0 = [X(1,:),u(1,:)];

if ST.analytical
    [~,y] = ode45(@(t,y) LorentzForce(t,y,q,m,ST.norm.l,ST.Bo,E),tspan,y0,options);  % NONSTIFF PROBLEM
else
    [~,y] = ode45(@(t,y) LorentzForceITER(t,y,q,m,ST.norm.l,ST.Bo,E,ST),tspan,y0,options);  % NONSTIFF PROBLEM
end
% [t,y] = ode15s(@(t,y) LorentzForce(t,y,q,m,ST.norm.l,ST.Bo,E),tspan,y0,options); % STIFF PROBLEM
% [t,y] = ode23t(@(t,y) LorentzForce(t,y,q,m,ST.norm.l,ST.Bo,E),tspan,y0,options);
% [t,y] = ode23tb(@(t,y) LorentzForce(t,y,q,m,ST.norm.l,ST.Bo,E),tspan,y0,options);
% [t,y] = ode23s(@(t,y) LorentzForce(t,y,q,m,ST.norm.l,ST.Bo,E),tspan,y0,options);

X = y(:,1:3);
u = y(:,4:6);
v(:,1) = u(:,1)./sqrt(1 + sum(u.^2,2));
v(:,2) = u(:,2)./sqrt(1 + sum(u.^2,2));
v(:,3) = u(:,3)./sqrt(1 + sum(u.^2,2));

X = X*ST.norm.l;

PP.X = X;
PP.v = v;

PP.EK = 1./sqrt(1-sum(PP.v.^2,2));
PP.ERR = 100*(PP.EK(1) - PP.EK)./PP.EK(1);

figure
plot3(PP.X(ST.params.inds,1),PP.X(ST.params.inds,2),PP.X(ST.params.inds,3))
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)

time = ST.time/(2*pi/ST.params.wc);

figure
plot(time(ST.params.inds), PP.ERR(ST.params.inds))
box on
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Energy conservation [\%]','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)

PP.v = PP.v*ST.params.c;

end

function dydt = LorentzForce(t,y,q,m,l,Bo,E)
% Lorentz Force
% y(1) = x;
% y(2) = y;
% y(3) = z;
% y(4) = ux;
% y(5) = uy;
% y(6) = uz;
B = analyticalB(y(1:3)*l)/Bo;

dydt = zeros(6,1);

gamma = sqrt( 1 + y(4)^2 + y(5)^2 + y(6)^2);

dydt(1) = y(4)/gamma; % x component of dx/dt = v
dydt(2) = y(5)/gamma; % y component of dx/dt = v
dydt(3) = y(6)/gamma; % y component of dx/dt = v

dydt(4) = (q/m)*(E(1) + (y(5)*B(3) - y(6)*B(2))/gamma );
dydt(5) = (q/m)*(E(2) + (y(6)*B(1) - y(4)*B(3))/gamma );
dydt(6) = (q/m)*(E(3) + (y(4)*B(2) - y(5)*B(1))/gamma );

end

function dydt = LorentzForceITER(t,y,q,m,l,Bo,E,ST)
% Lorentz Force
% y(1) = x;
% y(2) = y;
% y(3) = z;
% y(4) = ux;
% y(5) = uy;
% y(6) = uz;
B = interpMagFieldITER(ST,y(1:3)*l)/Bo;

dydt = zeros(6,1);

gamma = sqrt( 1 + y(4)^2 + y(5)^2 + y(6)^2);

dydt(1) = y(4)/gamma; % x component of dx/dt = v
dydt(2) = y(5)/gamma; % y component of dx/dt = v
dydt(3) = y(6)/gamma; % y component of dx/dt = v

dydt(4) = (q/m)*(E(1) + (y(5)*B(3) - y(6)*B(2))/gamma );
dydt(5) = (q/m)*(E(2) + (y(6)*B(1) - y(4)*B(3))/gamma );
dydt(6) = (q/m)*(E(3) + (y(4)*B(2) - y(5)*B(1))/gamma );

end

function BF = interpMagFieldITER(ST,X)
% Cylindrical coordinates
% R = radius, phi = azimuthal angle, Z = z coordinate
R = sqrt(X(1)^2 + X(2)^2);
Z = X(3);
phi = atan2(X(2),X(1));
if phi < 0
    phi = phi + 2*pi;
end

%     disp(['(R,phi,Z) = (' num2str(R) ',' num2str(phi) ',' num2str(Z) ')'])

BR = ST.B.SI.BR(R,Z);
Bphi = ST.B.SI.Bphi(R,Z);
BZ = ST.B.SI.BZ(R,Z);

Bx = BR*cos(phi) - Bphi*sin(phi);
By = BR*sin(phi) + Bphi*cos(phi);
Bz = BZ;
BF = [Bx,By,Bz];
end

function DI = DiegosInvariants(ST)
% All in Gaussian units
% Parameters of the analytical magnetic field
Bo = 1E4*ST.B.Bo;
a = 1E2*ST.B.a;% Minor radius in meters.
Ro = 1E2*ST.B.Ro; % Major radius in meters.
co = ST.B.co; % Extra parameter
lamb = a/co;
Bpo = 1E4*ST.B.Bpo;
% Parameters of the analytical magnetic field

c = 1E2*ST.params.c;
q = 3E9*ST.params.q;
m = 1E3*ST.params.m;
X = 1E2*ST.PP.X;
v = 1E2*ST.PP.v;

% Toroidal coordinates
% r = radius, theta = poloidal angle, phi = toroidal angle
r = sqrt( (sqrt(X(ST.params.inds,1).^2 + X(ST.params.inds,2).^2) - Ro).^2 + ...
    X(ST.params.inds,3).^2 );
theta = atan2(X(ST.params.inds,3),sqrt(X(ST.params.inds,1).^2 + X(ST.params.inds,2).^2) - Ro);
theta(theta<0) = theta(theta<0) + 2*pi;
zeta = atan2(X(ST.params.inds,1),X(ST.params.inds,2));
zeta(zeta<0) = zeta(zeta<0) + 2*pi;
% Toroidal coordinates

gamma = 1./sqrt(1 - sum(v(ST.params.inds,:).^2,2)/c^2);

eta = r./Ro;
psi = 0.5*lamb*Bpo*log(1 + r.^2/lamb^2);
wo = q*Bo./(m*c*gamma);

dzeta = ...
    (X(ST.params.inds,2).*v(ST.params.inds,1) - X(ST.params.inds,1).*v(ST.params.inds,2))./( sum(X(ST.params.inds,1:2).^2,2) );

DI = dzeta.*( 1 + eta.*cos(theta) ).^2 - wo.*psi/(Ro*Bo);

ERR = 100*(DI(1) - DI)./DI(1);

figure
subplot(2,1,1)
plot(ST.time(ST.params.inds),DI)
xlabel('Time $t$ [$\tau_c$]','Interpreter','latex','FontSize',16)
ylabel('Invariant','Interpreter','latex','FontSize',16)
title(ST.PP.method,'Interpreter','latex','FontSize',16)
subplot(2,1,2)
plot(ST.time(ST.params.inds),ERR)
xlabel('Time $t$ [$\tau_c$]','Interpreter','latex','FontSize',16)
ylabel('Relative error [\%]','Interpreter','latex','FontSize',16)
end

