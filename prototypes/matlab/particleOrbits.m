function ST = particleOrbits(pathToBField,fileType,ND,res,timeStepParams,tracerParams,xo,vo_params,opt)
% Here ro and vo are the initial position and velocity of the particle.
% The components of the initial velocity must be entered as fractions of
% the speed of light in vacuum.
%
% EXAMPLES:
% ST = particleOrbits('','','2D',[],[1E3,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% ST = particleOrbits('fields/ITER2D.dat','XPANDER','2D',[150,150],[1E3,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% ST = particleOrbits('fields/ITER3D.dat','XPANDER','3D',[150,100,150],[2E5,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% ST = particleOrbits('fields/ITER3D.dat','XPANDER','3D',[150,100,150],[1E5,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.1,70]);
% ST = particleOrbits('fields/CHEBYSHEV.dat','VMEC','2D',[60,60],[1E3,1.16E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% ST = particleOrbits('fields/VMEC2D.dat','VMEC','2D',[150,150],[1E3,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% ST = particleOrbits('fields/SIESTA.txt','SIESTA','3D',[100,99,149],[5E4,1E-2,1],[2,7.2938E3],[1.6,0,-0.5],[-0.03,0]);

narginchk(8,9);

% close all

ST = struct;
% Script parameters
if strcmp(pathToBField,'')
    ST.analytical = true; % true = analytical B field; false = xpander fields
    ST.pathToBField = pathToBField; % Path to xpander fields
else
    ST.analytical = false; % true = analytical B field; false = xpander fields
    ST.pathToBField = pathToBField; % Path to xpander fields
    ST.fileType = fileType;
    ST.res = res;
end

ST.ND = ND; % Dimensionality of xpander fields, 2D or 3D.
if nargin == 8
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
% ST.PP = particlePusherMatlab(ST);


PoincarePlots(ST);
% DiegosInvariants(ST);

munlock

end

% INITIALIZATION

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

function [b,a] = unitVectors(ST)
% initial condition of an electron drifting parallel to the local magnetic
% field.

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
B.fileType = ST.fileType;

data = load(ST.pathToBField);

switch ST.fileType
    case 'RAW' % This option changes quite frequently
        
        B.R = data(:,1);
        B.Z = data(:,3);
        
%         B.Ro = [B.R(1), B.Z(1)]; % This defines the position
        
        B.BR = data(:,4);
        B.Bphi = data(:,5); % minus sign
        B.BZ = data(:,6);
        
        B.B = sqrt(B.BR.^2 + B.BZ.^2 + B.Bphi.^2);
        
    case 'SIESTA'
        if strcmp(ST.ND,'2D')
            NR = ST.res(1);
            Nphi = ST.res(2);
            NZ = ST.res(3);
            
            R = data(:,1);
            Z = data(:,2);
            
            B.Ro = [R(1), Z(1)]; % This defines the position
            
            R = reshape(R,NR,Nphi,NZ);
            Z = reshape(Z,NR,Nphi,NZ);
            
            BR = data(:,4);
            BR = reshape(BR,NR,Nphi,NZ);
            Bphi = data(:,6);
            Bphi = reshape(Bphi,NR,Nphi,NZ);
            BZ = data(:,5);
            BZ = reshape(BZ,NR,Nphi,NZ);
            
            
            B.NR = NR;
            B.NZ = NZ;
            
            B.R = squeeze(R(:,1,:));
            B.Z = squeeze(Z(:,1,:));
            
            B.BR = squeeze(BR(:,1,:));
            B.Bphi = squeeze(Bphi(:,1,:));
            B.BZ = squeeze(BZ(:,1,:));
            
            B.B = sqrt(BR.^2 + BZ.^2 + Bphi.^2);
            
            B.Bo = mean(mean(mean(B.B)));
            
        elseif strcmp(ST.ND,'3D')
            B.NR = ST.res(1);
            B.Nphi = ST.res(2);
            B.NZ = ST.res(3);
            
            B.R = data(:,1);
            B.Z = data(:,2);
            B.phi = data(:,3);
            
            B.Ro = [B.R(1), B.Z(1)]; % This defines the position
            
            B.R = reshape(B.R,B.NR,B.Nphi,B.NZ);
            B.phi = reshape(B.phi,B.NR,B.Nphi,B.NZ);
            B.Z = reshape(B.Z,B.NR,B.Nphi,B.NZ);
            
            B.BR = data(:,4);
            B.BZ = data(:,5);
            B.Bphi = data(:,6);
            
            B.BR = reshape(B.BR,B.NR,B.Nphi,B.NZ);
            B.Bphi = reshape(B.Bphi,B.NR,B.Nphi,B.NZ);
            B.BZ = reshape(B.BZ,B.NR,B.Nphi,B.NZ);
            
            B.B = sqrt(B.BR.^2 + B.BZ.^2 + B.Bphi.^2);
            
            B.Bo = mean(mean(mean(B.B)));
        else
            error('Please use 2D or 3D fields');
        end
        
    case 'VMEC'
        
        if strcmp(ST.ND,'2D')
            B.NR = ST.res(1);
            B.NZ = ST.res(2);
            
            B.R = data(:,1);
            B.R = reshape(B.R,B.NR,B.NZ);
            B.Z = data(:,3);
            B.Z = reshape(B.Z,B.NR,B.NZ);
            
            B.BR = data(:,5);
            B.BR = reshape(B.BR,B.NR,B.NZ);
            B.Bphi = - data(:,6); % minus sign
            B.Bphi = reshape(B.Bphi,B.NR,B.NZ);
            B.BZ = data(:,7);
            B.BZ = reshape(B.BZ,B.NR,B.NZ);
            
        elseif strcmp(ST.ND,'3D')
            error('Not ready for using 3D fields of VMEC!')
        else
            error('Please use 2D or 3D fields');
        end
        
    case 'XPANDER'
        
        if strcmp(ST.ND,'2D')
            B.NR = ST.res(1);
            B.NZ = ST.res(2);
            
            B.R = data(:,1);
            B.R = reshape(B.R,B.NR,B.NZ);
            B.Z = data(:,3);
            B.Z = reshape(B.Z,B.NR,B.NZ);
            
            B.BR = data(:,4);
            B.BR = reshape(B.BR,B.NR,B.NZ);
            B.Bphi = data(:,5);
            B.Bphi = reshape(B.Bphi,B.NR,B.NZ);
            B.BZ = data(:,6);
            B.BZ = reshape(B.BZ,B.NR,B.NZ);
            
            B.P = data(:,7);
            B.P = reshape(B.P,B.NR,B.NZ);
        elseif strcmp(ST.ND,'3D')
            B.NR = ST.res(1);
            B.Nphi = ST.res(2);
            B.NZ = ST.res(3);
            
            B.R = zeros(B.NR,B.NZ,B.Nphi);
            B.Z = zeros(B.NR,B.NZ,B.Nphi);
            B.phi = zeros(B.NR,B.NZ,B.Nphi);
            
            phi = ((0:1:B.Nphi-1) + 0.5)*(2*pi/B.Nphi);
            
            B.B = zeros(B.NR,B.NZ,B.Nphi); % magnitude
            B.BR = zeros(B.NR,B.NZ,B.Nphi);
            B.Bphi = zeros(B.NR,B.NZ,B.Nphi);
            B.BZ = zeros(B.NR,B.NZ,B.Nphi);
            B.P = zeros(B.NR,B.NZ,B.Nphi);
            
            for iphi = 1:B.Nphi;
                for iz=1:B.NZ
                    indi = (iphi-1)*B.NR*B.NZ + (iz-1)*B.NR + 1;
                    indf = (iphi-1)*B.NR*B.NZ + iz*B.NR;
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
        else
            error('Please use 2D or 3D fields');
        end
        
    otherwise
        error('Unknown file format!')
end

B.B = sqrt(B.BR.^2 + B.BZ.^2 + B.Bphi.^2);

% Here the question what should be the characteristic magnetic field
% used in the normalization and calculation of the time step.
%     B.Bo = max(max(max(B.B)));
B.Bo = mean(mean(mean(B.B)));

% if ST.opt
%     plotLoadedMagneticField(B)
% end

% B.SI = calculateChebyshevInterpolant(ST,B);

B.SI = calculatescatteredInterpolant(ST,B);

B.R = [];
B.Z = [];
B.phi = [];

B.B = [];
B.BR = [];
B.Bphi = [];
B.BZ = [];
B.P = [];

end

function plotLoadedMagneticField(B)

switch B.fileType
    case 'RAW'
        S = 2*ones(numel(B.BR),1);
        
        figure
        subplot(1,4,1)      
        C = B.BR/max(abs(B.BR));
        scatter3(B.R,B.Z,B.BR,S,C)
        view([0,90])
        axis equal
        box on
        title('$B^R$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(1,4,2)
        C = B.Bphi/max(abs(B.Bphi));
        scatter3(B.R,B.Z,B.Bphi,S,C)
        view([0,90])
        axis equal
        box on
        title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(1,4,3)
        C = B.BZ/max(abs(B.BZ));
        scatter3(B.R,B.Z,B.BZ,S,C)
        view([0,90])
        axis equal
        box on
        title('$B^Z$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(1,4,4)
        C = B.B/max(abs(B.B));
        scatter3(B.R,B.Z,B.B,S,C)
        view([0,90])
        axis equal
        box on
        title('$B$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        colormap(jet)
    case 'SIESTA'
        figure
        subplot(2,2,1)
        surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.BR(:,1,:)),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B^R$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(2,2,2)
        surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.BZ(:,1,:)),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B^Z$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(2,2,3)
        surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.Bphi(:,1,:)),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        try
            subplot(2,2,4)
            surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.P(:,1,:)),'LineStyle','none')
            view([0,90])
            axis equal
            box on
            colorbar
            title('$P(R,Z)$','Interpreter','latex','FontSize',16)
            xlabel('R [m]','Interpreter','latex','FontSize',16)
            ylabel('Z [m]','Interpreter','latex','FontSize',16)
        catch
            subplot(2,2,4)
            surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.B(:,1,:)),'LineStyle','none')
            view([0,90])
            axis equal
            box on
            colorbar
            title('$B(R,Z)$','Interpreter','latex','FontSize',16)
            xlabel('R [m]','Interpreter','latex','FontSize',16)
            ylabel('Z [m]','Interpreter','latex','FontSize',16)
        end
        
        colormap(jet)
    otherwise
        figure
        subplot(2,2,1)
        surfc(B.R(:,:,1),B.Z(:,:,1),B.BR(:,:,1),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B^R$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(2,2,2)
        surfc(B.R(:,:,1),B.Z(:,:,1),B.BZ(:,:,1),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B^Z$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(2,2,3)
        surfc(B.R(:,:,1),B.Z(:,:,1),B.Bphi(:,:,1),'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        try
            subplot(2,2,4)
            surfc(B.R(:,:,1),B.Z(:,:,1),B.P(:,:,1),'LineStyle','none')
            view([0,90])
            axis equal
            box on
            colorbar
            title('$P(R,Z)$','Interpreter','latex','FontSize',16)
            xlabel('R [m]','Interpreter','latex','FontSize',16)
            ylabel('Z [m]','Interpreter','latex','FontSize',16)
        catch
            subplot(2,2,4)
            surfc(B.R(:,:,1),B.Z(:,:,1),B.B(:,:,1),'LineStyle','none')
            view([0,90])
            axis equal
            box on
            colorbar
            title('$B(R,Z)$','Interpreter','latex','FontSize',16)
            xlabel('R [m]','Interpreter','latex','FontSize',16)
            ylabel('Z [m]','Interpreter','latex','FontSize',16)
        end
        
        colormap(jet)
end

end

% INTERPOLANTS

function SI = calculateChebyshevInterpolant(ST,B)
% Cylindrical coordinates
% R = radius, phi = azimuthal angle, Z = z coordinate
% calculate interpolant of the field B
disp('Calculating chebfun2 interpolant...')
SI = struct;

if strcmp(ST.ND,'2D')   
    if mean(diff(B.R(:,1))) ~= 0 % non-uniform grid
        
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
        
        disp('chebfun2 interpolant: done!')
    elseif mean(diff(B.R(:,1))) == 0 % uniform grid
        SI = calculatescatteredInterpolant(ST,B);
    end
elseif strcmp(ST.ND,'3D')
    SI = calculatescatteredInterpolant(ST,B);
else
    error('Please, use 2D or 3D.');
end

end

function SI = calculatescatteredInterpolant(ST,B)
% Cylindrical coordinates
% R = radius, phi = azimuthal angle, Z = z coordinate
% calculate interpolant of the field B
disp('*** Uniform grid for the magnetic field detected! ***')
disp('Switching to scattered interpolant...')
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

    R = cat(3,B.R(:,:,end),B.R,B.R(:,:,1));
    Z = cat(3,B.Z(:,:,end),B.Z,B.Z(:,:,1));
    phi = cat(3,B.phi(:,:,end)-2*pi,B.phi,2*pi+B.phi(:,:,1));
    
    R = reshape(R,[numel(R) 1]);
    Z = reshape(Z,[numel(Z) 1]);
    phi = reshape(phi,[numel(phi) 1]);

    DATA = cat(3,B.BR(:,:,end),B.BR,B.BR(:,:,1));
    DATA = reshape(DATA,[numel(DATA) 1]);
%     SI.BR = scatteredInterpolant(R,Z,phi,DATA);
    SI.BR = scatteredInterpolant(R,Z,phi,DATA,'nearest','nearest');
    clear DATA
    
    DATA = cat(3,B.BZ(:,:,end),B.BZ,B.BZ(:,:,1));
    DATA = reshape(DATA,[numel(DATA) 1]);
%     SI.BZ = scatteredInterpolant(R,Z,phi,DATA);
    SI.BZ = scatteredInterpolant(R,Z,phi,DATA,'nearest','nearest');
    clear DATA
    
    DATA = cat(3,B.Bphi(:,:,end),B.Bphi,B.Bphi(:,:,1));
    DATA = reshape(DATA,[numel(DATA) 1]);
%     SI.Bphi = scatteredInterpolant(R,Z,phi,DATA);
    SI.Bphi = scatteredInterpolant(R,Z,phi,DATA,'nearest','nearest');
    clear DATA
    
else
    error('Use 2D or 3D');
end

disp('Scattered interpolant: done!')
end

% MODELS FOR MAGNETIC FIELD

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
    R = sqrt(X(1)^2 + X(2)^2);
    Z = X(3);
    phi = atan2(X(2),X(1));
    if phi < 0
        phi = phi + 2*pi;
    end
    
    BR = ST.B.SI.BR(R,Z,phi);
    Bphi = ST.B.SI.Bphi(R,Z,phi);
    BZ = ST.B.SI.BZ(R,Z,phi);
    
    Bx = BR*cos(phi) - Bphi*sin(phi);
    By = BR*sin(phi) + Bphi*cos(phi);
    Bz = BZ;
    BF = [Bx,By,Bz];
    
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

% LEAP-FROG PARTICLE PUSHER

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

X = zeros(ST.params.numIt,3); % (it,ii), ii=x,y,z
v = zeros(ST.params.numIt,3);
u = zeros(ST.params.numIt,3);
R = zeros(ST.params.numIt,3);

k = zeros(1,ST.params.numIt); % Curvature
T = zeros(1,ST.params.numIt); % Torsion

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
    
    % Curvature and torsion: 
    acc = q*cross(v(ii,:),B)/sqrt(1 + sum(u(ii,:).^2));
    aux = sum( cross(v(ii,:),acc).^2 );
    
    k(ii) = sqrt( aux )/sqrt( sum(v(ii,:).^2) )^3;
    
    dacc = q*cross(acc,B)/sqrt(1 + sum(u(ii,:).^2));
    
    T(ii) = det([v(ii,:); acc; dacc])/aux;
    
    %     disp(['Iteration ' num2str(ii)])
end

time = ST.time/(2*pi/ST.params.wc);

% Relative error in energy conservation
EK = 1./sqrt(1-sum(v.^2,2));
ERR = 100*(EK(1) - EK)./EK(1);
% Relative error in energy conservation

% Cylindrical coordinates
phi = atan2(X(:,2),X(:,1));
phi(phi < 0) = phi(phi < 0) + 2*pi;
% locs = [1;find(abs(diff(phi)) > 6);numel(phi)];

% Cylindrical coordinates

figure
subplot(2,1,1)
plot(time(ST.params.inds), ERR(ST.params.inds))
box on
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Energy conservation [\%]','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)
subplot(2,1,2)
plot(phi(ST.params.inds), ERR(ST.params.inds))
box on
xlabel('Azimuthal angle $\phi$ [rad]','Interpreter','latex','FontSize',16)
ylabel('Energy conservation [\%]','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)



% Return position and velocity with SI units
PP.X = X*ST.norm.l;
PP.R = R*ST.norm.l;
PP.v = v*ST.params.c;

PP.k = k/ST.norm.l; % curvature
PP.T = T/ST.norm.l; % curvature

figure
subplot(2,1,1)
plot(time(ST.params.inds), k(ST.params.inds))
box on
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Curvature $\kappa(t)$','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)
subplot(2,1,2)
plot(time(ST.params.inds), T(ST.params.inds),time(ST.params.inds),0*time(ST.params.inds),'k--')
box on
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Torsion $\tau(t)$','Interpreter','latex','FontSize',16)
title(PP.method,'Interpreter','latex','FontSize',16)

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

end

% MATLAB PARTICLE PUSHER

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

% POST-PROCESSING

function PoincarePlots(ST)
X = ST.PP.X;

R = sqrt(X(:,1).^2 + X(:,2).^2);
Z = X(:,3);

if isfield(ST.PP,'R')
    Rgc = sqrt(ST.PP.R(:,1).^2 + ST.PP.R(:,2).^2);
    Zgc = ST.PP.R(:,3);
end

zeta = atan2(X(:,2),X(:,1));
locs = find(abs(diff(zeta)) > 6);

figure
plot(R(locs),Z(locs),'r.','MarkerSize',15)
try 
    pol_angle = atan2(Z - ST.B.Ro(2),R - ST.B.Ro(1));
    locs = find(abs(diff(pol_angle)) > 6);
    hold on
    plot(R(locs(1):locs(2)),Z(locs(1):locs(2)),'k')
    if isfield(ST.PP,'R')
        plot(Rgc(locs(1):locs(2)),Zgc(locs(1):locs(2)),'g','LineWidth',2)
    end
    hold off
catch
    hold on
    plot(R,X(:,3),'k')
    if isfield(ST.PP,'R')
        plot(Rgc,Zgc,'g','LineWidth',2)
    end
    hold off
end
axis equal
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
title('Poincare plot','Interpreter','latex','FontSize',16)

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

