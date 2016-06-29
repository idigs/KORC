function ST = pOrbs(pathToBField,fileType,ND,res,timeStepParams,tracerParams,xo,vo,opt)
% Script to calculate particle orbits using either a relativistic modified
% leapfrog method or Matlab ODEs solvers.
%
% INPUT PARAMETERS' DESCRIPTION
% + pathToBField: absolute path to data files containing the magnetic
% field, for example, pathToBField='/home/Documents/jfit_165365_1400.txt'.
%
% + fileType: flag indicating the type of data file used in the simulation,
% different data files have different ways of presenting the data within
% the file. Its values can be: VMEC, XPANDER, SIESTA, or RAW (hard coded).
%
% + ND: dimensionality of the magnetic field, 2D for axisymmetric magnetic
% field, 3D for non-axisymmetric magnetic field.
%
% + res: a three-element vector specifying the dimensions along each coordinate
% of the magnetic field. For example, for a 3D XPANDER magnetic field 
% res = [NR,NPHI,NZ], where NR, NPHI and NZ are the number of nodes along
% the R, PHI, and Z axis (in cylindrical coordinates).
%
%  + timeStepParams: a three-element vector with parameters for controlling
%  the time integration of the particles' orbit. Here, 
% timeStepParams = [number_of_time_steps, time_step, cadence], where 
% 'time_step' is a time step in units of the particle gyro-period.
%
% + tracerParams: a two-element vector containing the charge of the
% particle in units of the absolute value of the electron charge and the 
% mass of the particle in terms of the electron mass, that is, 
% tracerParams=[charge,mass]. For example, for an electron tracerParams=[-1,-1]
%
% + xo: a three-dimensional vector in Cartesian coordinates specifying the
% initial position of the particle to be evolved in space and time.
%
% + vo_params: a two-dimensional vector containing the initial particle's
% velocity, in units of the speed of light, and the initial pitch angle.
% For example, for a particle initially streaming anti-parallel to the
% local magnetic field with a velocity of 0.5 the speed of light, 
% vo_params=[-0.5,0].
%
% + opt: (optional) logical flag specifying whether figures should be
% generated along the simulation. If opt=false, there will be no figures
% generated. All the information of the simulation will be stored in the
% matlab structure 'ST'.
%
% EXAMPLES:
% USING ANALYTICAL MAGNETIC FIELD
% ST = pOrbs('','','2D',[],[1E6,1E-2,10],[-1,1],[2.0,0,0],[50E6,170],true);
% USING TABULATED FIELDS OF THE ANALYTICAL MAGNETIC FIELD
% ST = pOrbs('fields/CHEBYSHEV.dat','VMEC','2D',[60,60],[1E3,1.16E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% USING XPAND FILES OF ITER FIELDS
% ST = pOrbs('fields/ITER2D.dat','XPANDER','2D',[150,150],[1E3,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% ST = pOrbs('fields/ITER3D.dat','XPANDER','3D',[150,100,150],[2E5,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% name = 'fields/xpand_iter3D_sc4_bet015_I87_hi_acc.dat';
% vo = 0.703270028800058;
% pitcho = 25.75;
% xo = 8.057345096062335;
% yo = -0.580890901345939;
% zo = 0.024044261000000;
% ST = pOrbs(name,'XPANDER','3D',[150,100,150],[5E6,1E-2,50],[-1,1],[xo,yo,zo],[vo,pitcho]);
% USING VMEC FILE OF ITER FIELDS
% ST = pOrbs('fields/VMEC2D.dat','VMEC','2D',[150,150],[1E3,1E-2,10],[2,7.2938E3],[6,0,-1],[-0.03,80]);
% USING SIESTA FILES
% ST = pOrbs('fields/SIESTA.txt','SIESTA','2D',[100,99,149],[5E5,1E-2,10],[-1,1],[1.6,0,-0.5],[0.99,0]);
% ST = pOrbs('fields/SIESTA_2.txt','SIESTA','3D',[100,99,149],[5E5,1E-2,10],[-1,1],[1.6,0,-0.5],[0.99,30]);
% name='fields/bfield_tracing/d3d_bfield_tracing-3E-3.dat';
% ST = pOrbs(name,'SIESTA','3D',[100,79,149],[1E7,1E-2,1],[-1,1],[1.5,0,-0.5],[0.99,30]);
% USING RAW FILES
% ST = pOrbs('jfit_165365_1400.txt','RAW','2D',[],[1E5,1E-2,10],[-1,1],[1.82,0,-0.4928],[0.99,70]);

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
ST.params.Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
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
vo(1) = sqrt(1 - (ST.params.me*ST.params.c^2/(vo(1)*ST.params.qe))^2);
ST.params.vo_params = vo; % vo_params = [velocity magnitude, pitch angle]
[ST.params.vo, ST.params.vpar, ST.params.vperp] = initializeVelocity(ST);

ST.coll.nimpurities = 1;
ST.coll.Te = 1.0*ST.params.qe; % Background electron temperature in eV
ST.coll.ne = 1.0E20; % Background electron density in 1/m^3
ST.coll.nH = ST.coll.ne;
ST.coll.Zo = [10.0]; % Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
ST.coll.Zj = [1.0]; % Average charge state of each impurity
ST.coll.nj = [5.0E20]; % Impurity densities
ST.coll.IZj = [21.5646]*ST.params.qe; % Ionization energy of impurity in eV

ST.coll.nef = ST.coll.ne + sum(ST.coll.Zj.*ST.coll.nj);
ST.coll.neb = ST.coll.ne + (ST.coll.Zo - ST.coll.Zj).*ST.coll.nj;
ST.coll.rD = sqrt( ST.params.ep*ST.coll.Te/(ST.coll.ne*ST.params.qe^2) );
ST.coll.re = ST.params.qe^2/( 4*pi*ST.params.ep*ST.params.me*ST.params.c^2 );
ST.coll.Ee_IZj = ST.params.me*ST.params.c^2/ST.coll.IZj;

% Particle's parameters
ST.params.q = tracerParams(1)*ST.params.qe; %alpha-particle
ST.params.m = tracerParams(2)*ST.params.me; % alpha-particle
ST.params.wc = sqrt(1 - sum(ST.params.vo.^2)/ST.params.c^2)*abs(ST.params.q)*ST.Bo/ST.params.m;
% ST.params.rL = ST.params.vperp/ST.params.wc; % Larmor radius of particle

if ST.opt
    disp(['Cyclotron frequency ' num2str(ST.params.wc) ' Hz'])
end

% Normalisation parameters
ST.norm.q = abs(ST.params.q);
ST.norm.m = ST.params.m;
ST.norm.wc = ST.norm.q*ST.Bo/ST.norm.m;
ST.norm.t = 1/ST.norm.wc;
ST.norm.l = ST.params.c*ST.norm.t;
ST.norm.E = ST.norm.m*ST.params.c^2;
ST.norm.n = 1/ST.norm.l^3;
% Normalisation parameters

% Numerical parameters
ST.params.numIt = timeStepParams(1);
ST.params.dt = timeStepParams(2)*(2*pi/ST.params.wc);
ST.params.cadence = timeStepParams(3);
ST.params.inds = 1:ST.params.cadence:ST.params.numIt;
ST.params.numSnapshots = numel(ST.params.inds);

ST.time = zeros(ST.params.numSnapshots,1);
for ii=1:ST.params.numSnapshots
    ST.time(ii) = (ii-1)*ST.params.dt*ST.params.cadence;
end

ST.PP = particlePusherLeapfrog(ST);

% Particle pusher using matlab ODEs solvers
% ST.PP = particlePusherMatlab(ST); % (Un)comment this line as required
% Particle pusher using matlab ODEs solvers

if ST.opt
    PoincarePlots(ST);
end
% ST.PP.angularMomentum = DiegosInvariants(ST);

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

if ST.opt
    disp(['Parallel velocity: ' num2str(vpar/ST.params.c)])
    disp(['Perpendicular velocity: ' num2str(vperp/ST.params.c)])
end
end

function [b,a] = unitVectors(ST)
% initial condition of an electron drifting parallel to the local magnetic
% field.

tol = 1E-10;

if ST.analytical
    B = analyticalB(ST.params.Xo);
    b = B/sqrt(sum(B.^2));
else
    B = interpMagField(ST,ST.params.Xo);
    b = B/sqrt(sum(B.^2));
end

rng('shuffle')

if all(b)
    c1 = b(1)^2 + b(2)^2;
    
    az = sqrt(c1)*rand;
    
    c2 = 2*b(1)*b(3)*az;
    c3 = (b(2)^2 + b(3)^2)*az^2 -b(2)^2;
    
    ax = (-c2 + sqrt(c2^2 - 4*c1*c3))/(2*c1);
    ay = sqrt(1 - ax^2 - az^2);
    
    a = [ax, ay, az];
    
    if abs(dot(a,b)) > tol      
        ax = (-c2 - sqrt(c2^2 - 4*c1*c3))/(2*c1);
        ay = sqrt(1 - ax^2 - az^2);
        
        a = [ax, ay, az];
    end
    
    if abs(dot(a,b)) > tol
        az = -az;
        
        c2 = 2*b(1)*b(3)*az;
        c3 = (b(2)^2 + b(3)^2)*az^2 -b(2)^2;
        
        ax = (-c2 + sqrt(c2^2 - 4*c1*c3))/(2*c1);
        ay = sqrt(1 - ax^2 - az^2);
        
        a = [ax, ay, az];
        
        if abs(dot(a,b)) > tol
            ax = (-c2 - sqrt(c2^2 - 4*c1*c3))/(2*c1);
            ay = sqrt(1 - ax^2 - az^2);
            
            a = [ax, ay, az];
        end
    end
elseif b(1) == 0
    ay = rand/sqrt( 1 + (b(2)/b(3))^2 );
    az = -ay*(b(2)/b(z));
    ax = sqrt( 1 - ay^2*( 1 + (b(2)/b(3))^2 ) );
    
    a = [ax, ay, az];
elseif b(2) == 0
    az = rand/sqrt( 1 + (b(3)/b(1))^2 );
    ax = -az*(b(3)/b(1));
    ay = sqrt( 1 - az^2*( 1 + (b(3)/b(1))^2 ) );
    
    a = [ax, ay, az];
elseif b(3) == 0
    ay = rand/sqrt( 1 + (b(2)/b(1))^2 );
    ax = -ay*(b(2)/b(1));
    az = sqrt( 1 - ay^2*( 1 + (b(2)/b(1))^2 ) );
    
    a = [ax, ay, az];
end

if ~isreal(a)
    error('Imaginary values in perpendicular vector to b');
end

if  abs(dot(a,b)) > tol
    error('a.b ~= 0')
end

% figure
% plot3([0,b(1)],[0,b(2)],[0,b(3)],'r',[0,a(1)],[0,a(2)],[0,a(3)],'b')

end

function b = unitParallelVector(ST,X)
% initial condition of an electron drifting parallel to the local magnetic
% field.

if ST.analytical
    B = analyticalB(X);
    b = B/sqrt(sum(B.^2));
else
    B = interpMagField(ST,X);
    b = B/sqrt(sum(B.^2));
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
        B.Bphi = - data(:,5); % minus sign
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
            BZ = data(:,5);
            Bphi = data(:,6);
            
            BR = reshape(BR,NR,Nphi,NZ);
            Bphi = reshape(Bphi,NR,Nphi,NZ);
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
            
            R = data(:,1);
            Z = data(:,2);
%             phi = data(:,3);
            phi = ((0:1:B.Nphi-1) + 0.5)*(2*pi/B.Nphi);
            
            B.Ro = [R(1), Z(1)]; % This defines the position
            
            R = reshape(R,B.NR,B.Nphi,B.NZ);
%             phi = reshape(phi,B.NR,B.Nphi,B.NZ);
            Z = reshape(Z,B.NR,B.Nphi,B.NZ);
            
            BR = data(:,4);
            BZ = data(:,5);
            Bphi = data(:,6);
            
            BR = reshape(BR,B.NR,B.Nphi,B.NZ);
            Bphi = reshape(Bphi,B.NR,B.Nphi,B.NZ);
            BZ = reshape(BZ,B.NR,B.Nphi,B.NZ);
            
            B.R = zeros(B.NR,B.NZ,B.Nphi);
            B.phi = zeros(B.NR,B.NZ,B.Nphi);
            B.Z = zeros(B.NR,B.NZ,B.Nphi);
            B.BR = zeros(B.NR,B.NZ,B.Nphi);
            B.Bphi = zeros(B.NR,B.NZ,B.Nphi);
            B.BZ = zeros(B.NR,B.NZ,B.Nphi);
            
            for iphi=1:B.Nphi
                B.R(:,:,iphi) = squeeze(R(:,iphi,:));
%                 B.phi(:,:,iphi) = squeeze(phi(:,iphi,:));
                B.phi(:,:,iphi) = phi(iphi);
                B.Z(:,:,iphi) = squeeze(Z(:,iphi,:));
                
                B.BR(:,:,iphi) = squeeze(BR(:,iphi,:));
                B.Bphi(:,:,iphi) = squeeze(Bphi(:,iphi,:));
                B.BZ(:,:,iphi) = squeeze(BZ(:,iphi,:));
            end
            
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
                    B.Bphi(:,iz,iphi) = - data(indi:indf,5);
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

if ST.opt
    plotLoadedMagneticField(B)
end

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
%     case 'SIESTA'
%         figure
%         subplot(2,2,1)
%         surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.BR(:,1,:)),'LineStyle','none')
%         view([0,90])
%         axis equal
%         box on
%         colorbar
%         title('$B^R$ [T]','Interpreter','latex','FontSize',16)
%         xlabel('R [m]','Interpreter','latex','FontSize',16)
%         ylabel('Z [m]','Interpreter','latex','FontSize',16)
%         
%         subplot(2,2,2)
%         surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.BZ(:,1,:)),'LineStyle','none')
%         view([0,90])
%         axis equal
%         box on
%         colorbar
%         title('$B^Z$ [T]','Interpreter','latex','FontSize',16)
%         xlabel('R [m]','Interpreter','latex','FontSize',16)
%         ylabel('Z [m]','Interpreter','latex','FontSize',16)
%         
%         subplot(2,2,3)
%         surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.Bphi(:,1,:)),'LineStyle','none')
%         view([0,90])
%         axis equal
%         box on
%         colorbar
%         title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
%         xlabel('R [m]','Interpreter','latex','FontSize',16)
%         ylabel('Z [m]','Interpreter','latex','FontSize',16)
%         
%         try
%             subplot(2,2,4)
%             surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.P(:,1,:)),'LineStyle','none')
%             view([0,90])
%             axis equal
%             box on
%             colorbar
%             title('$P(R,Z)$','Interpreter','latex','FontSize',16)
%             xlabel('R [m]','Interpreter','latex','FontSize',16)
%             ylabel('Z [m]','Interpreter','latex','FontSize',16)
%         catch
%             subplot(2,2,4)
%             surfc(squeeze(B.R(:,1,:)),squeeze(B.Z(:,1,:)),squeeze(B.B(:,1,:)),'LineStyle','none')
%             view([0,90])
%             axis equal
%             box on
%             colorbar
%             title('$B(R,Z)$','Interpreter','latex','FontSize',16)
%             xlabel('R [m]','Interpreter','latex','FontSize',16)
%             ylabel('Z [m]','Interpreter','latex','FontSize',16)
%         end
%         
%         colormap(jet)
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
if ST.opt
    disp('Calculating chebfun2 interpolant...')
end
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
        
        if ST.opt
            figure
            subplot(1,3,1)
            plot(SI.BR)
            subplot(1,3,2)
            plot(SI.Bphi)
            subplot(1,3,3)
            plot(SI.BZ)
            colormap(jet)
        end
        
        if ST.opt
            disp('chebfun2 interpolant: done!')
        end
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
if ST.opt
    disp('*** Uniform grid for the magnetic field detected! ***')
    disp('Switching to scattered interpolant...')
end
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
    SI.BR = scatteredInterpolant(R,Z,phi,DATA);
%     SI.BR = scatteredInterpolant(R,Z,phi,DATA,'nearest','nearest');
    clear DATA
    
    DATA = cat(3,B.BZ(:,:,end),B.BZ,B.BZ(:,:,1));
    DATA = reshape(DATA,[numel(DATA) 1]);
    SI.BZ = scatteredInterpolant(R,Z,phi,DATA);
%     SI.BZ = scatteredInterpolant(R,Z,phi,DATA,'nearest','nearest');
    clear DATA
    
    DATA = cat(3,B.Bphi(:,:,end),B.Bphi,B.Bphi(:,:,1));
    DATA = reshape(DATA,[numel(DATA) 1]);
    SI.Bphi = scatteredInterpolant(R,Z,phi,DATA);
%     SI.Bphi = scatteredInterpolant(R,Z,phi,DATA,'nearest','nearest');
    clear DATA
    
else
    error('Use 2D or 3D');
end

if ST.opt
    disp('Scattered interpolant: done!')
end
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
Bo = 2.19;
a = 0.6;% Minor radius in meters.
Ro = 1.695; % Major radius in meters.
qa = 3; % Safety factor at the separatrix (r=a)
co = 0.5; % Extra parameter
lamb = a/co;
Bpo = (a/Ro)*(Bo/qa)*(1+co^2)/co;
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
    % Minus sign = TEXTOR
    % Plus sign = default
    Bp = Bpo*(r/lamb)/( 1 + (r/lamb)^2 );
    
    eta = r/Ro;
    Br = 1/( 1 + eta*cos(theta) );
    
    Bx = Br*( Bo*cos(zeta) - Bp*sin(theta)*sin(zeta) );
    By = -Br*( Bo*sin(zeta) + Bp*sin(theta)*cos(zeta) );
    Bz = Br*Bp*cos(theta);
    
    B = [Bx,By,Bz];
end

end

function E = analyticalE(X)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

narginchk(1,2);

% Parameters of the analytical magnetic field
Eo = -0.0;
Ro = 1.695; % Major radius in meters.
% Parameters of the analytical magnetic field

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
eta = r/Ro;
Ezeta = Eo/( 1 + eta*cos(theta) );

Ex = Ezeta*cos(zeta);
Ey = -Ezeta*sin(zeta);
Ez = 0;

E = [Ex,Ey,Ez];

end

% LEAP-FROG PARTICLE PUSHER

function PP = particlePusherLeapfrog(ST)
% Relativistic particle pusher.
% Vay, J.-L. 2008, PoP 15, 056701.
if ST.opt
    disp('Initializing particle pusher...')
end

PP = struct;

PP.method = 'Leapfrog';

X = zeros(ST.params.numSnapshots,3); % (it,ii), ii=x,y,z
v = zeros(ST.params.numSnapshots,3);
u = zeros(ST.params.numSnapshots,3);
R = zeros(ST.params.numSnapshots,3);

POINCARE = struct;
POINCARE.R = [];
POINCARE.Z = [];
POINCARE.k = [];
POINCARE.T = [];
POINCARE.pitch = [];

k = zeros(1,ST.params.numSnapshots); % Curvature
T = zeros(1,ST.params.numSnapshots); % Torsion
vpar = zeros(1,ST.params.numSnapshots); % parallel velocity
vperp = zeros(1,ST.params.numSnapshots); % perpendicular velocity
mu = zeros(1,ST.params.numSnapshots); % instantaneous magnetic moment
EK = zeros(1,ST.params.numSnapshots); % kinetic energy

fL = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
f2 = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
f3 = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
fcolle = zeros(1,ST.params.numSnapshots); % Collision force with electrons
fcolli = zeros(1,ST.params.numSnapshots); % Collision force with ions

WL = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
WR = zeros(1,ST.params.numSnapshots); % Synchroton radiated power

% Normalization
X(1,:) = ST.params.Xo/ST.norm.l;
v(1,:) = ST.params.vo/ST.params.c;
q = ST.params.q/ST.norm.q;
m = ST.params.m/ST.norm.m;
if ST.analytical
    B = analyticalB(X(1,:)*ST.norm.l)/ST.Bo;
    E = analyticalE(X(1,:)*ST.norm.l)/(ST.Bo*ST.params.c);
else
    B = interpMagField(ST,X(1,:)*ST.norm.l)/ST.Bo;
    E = [0,0,0];
end
dt = ST.params.dt*ST.norm.wc;
ep = ST.params.ep*(ST.norm.m^2*ST.params.c^3/(ST.norm.q^3*ST.Bo));

ST.coll.Te = ST.coll.Te/ST.norm.E;
ST.coll.ne = ST.coll.ne/ST.norm.n;
ST.coll.nH = ST.coll.nH/ST.norm.n;
ST.coll.nj = ST.coll.nj/ST.norm.n;
ST.coll.IZj = ST.coll.IZj/ST.norm.E;

ST.coll.nef = ST.coll.nef/ST.norm.n;
ST.coll.neb = ST.coll.neb/ST.norm.n;
ST.coll.rD = ST.coll.rD/ST.norm.l;
ST.coll.re = ST.coll.re/ST.norm.l;
% Normalization

% initial velocity at time level t = 0
u(1,:) = v(1,:)/sqrt(1 - sum(v(1,:).^2));
gamma = sqrt(1 + sum(u(1,:).^2));
v(1,:) = u(1,:)/gamma;
R(1,:) = X(1,:) + gamma*m*cross(v(1,:),B)/(q*sum(B.^2));
EK(1) = gamma;

% % % % % % % % % % % % % % % % % % 
B_mag = sqrt(B*B');
b = B/B_mag;
vpar(1) = v(1,:)*b';
vperp(1) = sqrt( v(1,:)*v(1,:)' - vpar(1)^2 );
mu(1) = m*gamma^2*vperp(1)^2/(2*B_mag);

% Curvature and torsion
vmag = sqrt( sum(v(1,:).^2) );
aux =  cross(v(1,:),E) + v(1,:)*sum(v(1,:).*B) - B*vmag^2;
k(1) = abs(q)*sqrt( sum(aux.^2) )/(gamma*m*vmag^3);
% Curvature and torsion

% % % Leap-frog scheme for the radiation damping force % % %
F2 = ( q^3/(6*pi*ep*m^2) )*( (E*v(1,:)')*E + cross(E,B) +...
    cross(B,cross(B,v(1,:))) );
vec = E + cross(v(1,:),B);
F3 = ( gamma^2*q^3/(6*pi*ep*m^2) )*( (E*v(1,:)')^2 - vec*vec' )*v(1,:);

WL(1) = q*(E*v(1,:)');
WR(1) = ( q^4/(6*pi*ep*m^2) )*( E*E' + cross(v(1,:),B)*E' +...
    gamma^2*( (E*v(1,:)')^2 - vec*vec' ) );


tmp = (gamma - 1.0)*sqrt(gamma + 1.0);
Clog_ef = log(0.5*tmp*(ST.coll.rD/ST.coll.re)/gamma);
ae = ST.coll.nef*Clog_ef;
for ppi=1:ST.coll.nimpurities
    Clog_eb = log(tmp*ST.coll.Ee_IZj(ppi));
    ae = ae + ST.coll.neb(ppi)*Clog_eb;
end

tmp = (gamma^2 - 1.0)/gamma;
Clog_eH = log( tmp*(ST.coll.rD/ST.coll.re) );
ai = ST.coll.nH*Clog_eH;
for ppi=1:ST.coll.nimpurities
    Clog_eZj = log( ST.coll.rD/(ST.coll.Zj(ppi)*ST.coll.re*ST.coll.Ee_IZj(ppi)) );
    Clog_eZo = log(tmp*ST.coll.Ee_IZj(ppi));
    ai = ai + ST.coll.nj(ppi)*(Clog_eZj*ST.coll.Zj(ppi)^2 + Clog_eZo*ST.coll.Zo(ppi)^2);
end

tmp = gamma*(gamma + 1.0)/sqrt(u(1,:)*u(1,:)')^3;
Fcolle = -4.0*pi*ae*m*(ST.coll.re^2)*tmp*u(1,:);

tmp = gamma/sqrt(u(1,:)*u(1,:)')^3;
Fcolli = -4.0*pi*ai*m*(ST.coll.re^2)*tmp*u(1,:);
% Collisions
% % % Leap-frog scheme for the radiation damping force % % %

fL(1) = abs(q)*sqrt( vec*vec' );
f2(1) = abs(q)*sqrt( F2*F2' );
f3(1) = abs(q)*sqrt( F3*F3' );
fcolle(1) = abs(q)*sqrt( Fcolle*Fcolle' );
fcolli(1) = abs(q)*sqrt( Fcolli*Fcolli' );
% % % % % % % % % % % % % % % % % % 

% initial half-step for position
DT = 0.5*dt;
V = v(1,:);
U = u(1,:);
XX = X(1,:) + DT*V;
% initial half-step for position

a = q*dt/m;

for ii=2:ST.params.numSnapshots
    for jj=1:ST.params.cadence
        
        if ST.analytical
            B = analyticalB(XX*ST.norm.l)/ST.Bo;
            E = analyticalE(XX*ST.norm.l)/(ST.Bo*ST.params.c);
        else
            B = interpMagField(ST,XX*ST.norm.l)/ST.Bo;
        end
        
        U_R = U;
        U_L = U;
        
        % % % Leap-frog scheme for Lorentz force % % %
        U_hs = U_L + 0.5*a*(E + cross(V,B));
        
        % % % % % % % % % % % % % % %
        % Radiation losses operator
        vmag = sqrt( V*V' );
        aux =  cross(V,E) + V*(V*B') - B*vmag^2;
        curv = abs(q)*sqrt( aux*aux' )/(gamma*m*vmag^3);
        % Radiation losses operator
        % % % % % % % % % % % % % % %
        
        tau = 0.5*q*dt*B/m;
        up = U_hs + 0.5*a*E;
        gammap = sqrt(1 + up*up');
        sigma = gammap^2 - tau*tau';
        us = up*tau'; % variable 'u^*' in paper
        gamma = sqrt(0.5)*sqrt( sigma + sqrt(sigma^2 + 4*(tau*tau' + us^2)) );
        t = tau/gamma;
        s = 1/(1 + t*t'); % variable 's' in paper
        
        U_L = s*(up + (up*t')*t + cross(up,t));
        V = U_L/gamma;
        % % % Leap-frog scheme for Lorentz force % % %
        
        % % % Leap-frog scheme for the radiation damping force % % %
        U_eff = 0.5*(U_L + U);
        gamma_eff = sqrt(1 + U_eff*U_eff');
        V_eff = U_eff/gamma_eff;
        
        F2 = ( q^3/(6*pi*ep*m^2) )*( (E*V_eff')*E + cross(E,B) +...
            cross(B,cross(B,V_eff)) );
        vec = E + cross(V_eff,B);
        F3 = ( gamma_eff^2*q^3/(6*pi*ep*m^2) )*( (E*V_eff')^2 - vec*vec' )*V_eff;
        
        % Collisions
        tmp = (gamma_eff - 1.0)*sqrt(gamma_eff + 1.0);
        Clog_ef = log(0.5*tmp*(ST.coll.rD/ST.coll.re)/gamma_eff);
        ae = ST.coll.nef*Clog_ef;
        for ppi=1:ST.coll.nimpurities
            Clog_eb = log(tmp*ST.coll.Ee_IZj(ppi));
            ae = ae + ST.coll.neb(ppi)*Clog_eb;
        end
        
        
        tmp = (gamma_eff^2 - 1.0)/gamma_eff;
        Clog_eH = log( tmp*(ST.coll.rD/ST.coll.re) );
        ai = ST.coll.nH*Clog_eH;
        for ppi=1:ST.coll.nimpurities
            Clog_eZj = log( ST.coll.rD/(ST.coll.Zj(ppi)*ST.coll.re*ST.coll.Ee_IZj(ppi)) );
            Clog_eZo = log(tmp*ST.coll.Ee_IZj(ppi));
            ai = ai + ST.coll.nj(ppi)*(Clog_eZj*ST.coll.Zj(ppi)^2 + Clog_eZo*ST.coll.Zo(ppi)^2);
        end
        
        tmp = gamma_eff*(gamma_eff + 1.0)/sqrt(U_eff*U_eff')^3;
        Fcolle = -4.0*pi*ae*m*(ST.coll.re^2)*tmp*U_eff/q;
        
        tmp = gamma_eff/sqrt(U_eff*U_eff')^3;
        Fcolli = -4.0*pi*ai*m*(ST.coll.re^2)*tmp*U_eff/q;
        % Collisions
         
%         U_R = U_R + a*( F2 + F3 );
        U_R = U_R + a*( F2 + F3 + Fcolle + Fcolli);
        
        U = U_L + U_R - U;
        gamma = sqrt( 1 + U*U' );
        V = U/gamma;
        
%         U = U_L;
%         V = U_L/gamma;
        % % % Leap-frog scheme for the radiation damping force % % %fcoll
        
        zeta_previous = atan2(XX(2),XX(1));
        if zeta_previous < 0
            zeta_previous = zeta_previous + 2*pi;
        end
        XX = XX + dt*V;
        
        zeta_current = atan2(XX(2),XX(1));
        if zeta_current < 0
            zeta_current = zeta_current + 2*pi;
        end
        
        if abs(zeta_previous - zeta_current) > 6
            
            POINCARE.R = [POINCARE.R; sqrt(XX(1).^2 + XX(2).^2)];
            POINCARE.Z = [POINCARE.Z; XX(3)];
            
            b = B/sqrt(B*B');
            Vpar = V*b';
            Vperp = sqrt( V*V' - Vpar^2 );
            
            POINCARE.k = [POINCARE.k; curv];
            POINCARE.T = [POINCARE.T; 0.0];
            POINCARE.pitch = [POINCARE.pitch; atan2(Vperp,Vpar)];
        end
    end
    
    X(ii,:) = XX; % Time level n+1/2
    u(ii,:) = U; % Time level n+1
    v(ii,:) = V; % Time level n+1
    EK(ii) = gamma; % Time level n+1
    
    gamma_half_step = sqrt(1 + U_hs*U_hs');
    v_half_step = U_hs/gamma_half_step;
    
    R(ii,:) = X(ii,:) + gamma*m*cross(v_half_step,B)/(q*(B*B'));
    
    B_mag = sqrt(B*B');
    b = B/B_mag;
    vpar(ii) = v(ii,:)*b';
    vperp(ii) = sqrt( v(ii,:)*v(ii,:)' - vpar(ii)^2 );
    mu(ii) = m*gamma^2*vperp(ii)^2/(2*B_mag);
    
    k(ii) = curv;
    
    fL(ii) = abs(q)*sqrt( vec*vec' );
    f2(ii) = abs(q)*sqrt( F2*F2' );
    f3(ii) = abs(q)*sqrt( F3*F3' );
    fcolle(ii) = abs(q)*sqrt( Fcolle*Fcolle' );
    fcolli(ii) = abs(q)*sqrt( Fcolli*Fcolli' );
    
    F2 = ( q^3/(6*pi*ep*m^2) )*( (E*V_eff')*E + cross(E,B) +...
        cross(B,cross(B,V_eff)) );
    vec = E + cross(V_eff,B);
    F3 = ( gamma_eff^2*q^3/(6*pi*ep*m^2) )*( (E*V_eff')^2 - vec*vec' )*V_eff;
    
    WL(ii) = q*(E*V_eff');
    WR(ii) = ( q^4/(6*pi*ep*m^2) )*( E*E' + cross(V_eff,B)*E' +...
        gamma_eff^2*( (E*V_eff')^2 - vec*vec' ) );
end


if ST.opt
    time = ST.time;%/(2*pi/ST.params.wc);
    
    % Relative error in energy conservation
    ERR_EK = 100*(EK - EK(1))./EK(1);
    ERR_mu = 100*(mu - mu(1))./mu(1);
    % Relative error in energy conservation
    
    x = linspace(min(time),max(time),10);
    y = zeros(1,10);
    
    figure
    subplot(4,1,1)
    plot(time, EK)
    box on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$\gamma m c^2$ [$m_e c^2$]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(4,1,2)
    plot(time, ERR_EK,x,y,'k--')
    box on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('Error in $\gamma m c^2$ [\%]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(4,1,3)
    plot(time, mu)
    box on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$\mu$ [$c^2 m_o/B_o$]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(4,1,4)
    plot(time, ERR_mu,x,y,'k--')
    box on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('Error in $\mu$ [\%]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    
    
    figure
    subplot(3,1,1)
    plot(time, vpar)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$v_\parallel$ [c]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(3,1,2)
    plot(time, vperp)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$v_\perp$ [c]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    v_tot = sqrt(vpar.^2 +vperp.^2);
    subplot(3,1,3)
    plot(time, v_tot)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$v$ [c]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
end

% Return position and velocity with SI units
PP.X = X*ST.norm.l;
PP.R = R*ST.norm.l;
PP.v = v*ST.params.c;

PP.vpar = vpar;
PP.vperp = vperp;

PP.k = k/ST.norm.l; % curvature
PP.T = T/ST.norm.l; % curvature

PP.EK = EK;
PP.mu = mu;

POINCARE.R = POINCARE.R*ST.norm.l;
POINCARE.Z = POINCARE.Z*ST.norm.l;
POINCARE.k = POINCARE.k/ST.norm.l;
POINCARE.T = POINCARE.T/ST.norm.l;

PP.POINCARE = POINCARE;


if ST.opt
    figure
    subplot(2,1,1)
    plot(time, PP.k)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('Curvature $\kappa(t)$','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(2,1,2)
    plot(time, PP.T,time,0*time,'k--')
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('Torsion $\tau(t)$','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    
    figure
    plot3(PP.X(:,1),PP.X(:,2),PP.X(:,3),'b')
    hold on
    plot3(PP.R(:,1),PP.R(:,2),PP.R(:,3),'r')
    hold off
    axis equal
    box on
    xlabel('X','Interpreter','latex','FontSize',16)
    ylabel('Y','Interpreter','latex','FontSize',16)
    zlabel('Z','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    
    figure
    subplot(4,1,1)
    plot(time, fL)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|F_L|$','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(4,1,2)
    plot(time, f2)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|f_2|$','Interpreter','latex','FontSize',16)
    subplot(4,1,3)
    plot(time, f3)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|f_3|$','Interpreter','latex','FontSize',16)
    subplot(4,1,4)
    hold on
    plot(time, fcolle,'r')
    plot(time, fcolli,'b')
    hold off
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|f_{coll}|$','Interpreter','latex','FontSize',16)
    
    figure
    subplot(3,1,1)
    plot(time, WL)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|W_L|$','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(3,1,2)
    plot(time, WR)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|W_R|$','Interpreter','latex','FontSize',16)
    subplot(3,1,3)
    plot(time, abs(WR./WL))
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|W_R/W_L|$','Interpreter','latex','FontSize',16)
    disp(['Mean value: ' num2str(mean(WR./WL))])
end

if ST.opt
    disp('Particle pusher: done!')
end
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

options = odeset('RelTol',1E-6,'AbsTol',1E-10);%,'Stats','on','OutputFcn',@odeplot)
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

if ST.opt
    figure
    plot3(PP.X(ST.params.inds,1),PP.X(ST.params.inds,2),PP.X(ST.params.inds,3))
    axis equal
    xlabel('X','Interpreter','latex','FontSize',16)
    ylabel('Y','Interpreter','latex','FontSize',16)
    zlabel('Z','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    
    
    time = ST.time/(2*pi/ST.params.wc);
    
    figure
    plot(time, PP.ERR)
    box on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('Energy conservation [\%]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
end

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
zeta(zeta<0) = zeta(zeta<0) + 2*pi;
locs = find(abs(diff(zeta)) > 6);

figure
plot(R(locs),Z(locs),'r.','MarkerSize',15)
hold on
if isfield(ST.PP,'R')
    plot(Rgc(locs),Zgc(locs),'k.','MarkerSize',15)
end
hold off
% try 
%     pol_angle = atan2(Z - ST.B.Ro(2),R - ST.B.Ro(1));
% %     pol_angle(pol_angle<0) = pol_angle(pol_angle<0) + 2*pi;
%     locs = find(abs(diff(pol_angle)) > 6);
%     hold on
%     plot(R(locs(1):locs(2)),Z(locs(1):locs(2)),'k')
%     if isfield(ST.PP,'R')
%         plot(Rgc(locs(1):locs(2)),Zgc(locs(1):locs(2)),'g')
%     end
%     hold off
% catch
    hold on
    plot(R,X(:,3),'k')
    if isfield(ST.PP,'R')
        plot(Rgc,Zgc,'g','LineWidth',2)
    end
    hold off
% end
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
r = sqrt( (sqrt(X(:,1).^2 + X(:,2).^2) - Ro).^2 + X(:,3).^2 );
theta = atan2(X(:,3),sqrt(X(:,1).^2 + X(:,2).^2) - Ro);
theta(theta<0) = theta(theta<0) + 2*pi;
zeta = atan2(X(:,1),X(:,2));
zeta(zeta<0) = zeta(zeta<0) + 2*pi;
% Toroidal coordinates

gamma = 1./sqrt(1 - sum(v.^2,2)/c^2);

eta = r./Ro;
psi = 0.5*lamb*Bpo*log(1 + r.^2/lamb^2);
wo = q*Bo./(m*c*gamma);

dzeta = ...
    (X(:,2).*v(:,1) - X(:,1).*v(:,2))./( sum(X(:,1:2).^2,2) );

DI = dzeta.*( 1 + eta.*cos(theta) ).^2 - wo.*psi/(Ro*Bo);

ERR = 100*(DI(1) - DI)./DI(1);

time = ST.time/(2*pi/ST.norm.wc);

figure
subplot(2,1,1)
plot(ST.time(:),DI)
xlabel('Time $t$ [$\tau_c$]','Interpreter','latex','FontSize',16)
ylabel('Invariant','Interpreter','latex','FontSize',16)
title(ST.PP.method,'Interpreter','latex','FontSize',16)
subplot(2,1,2)
plot(ST.time(:),ERR)
xlabel('Time $t$ [$\tau_c$]','Interpreter','latex','FontSize',16)
ylabel('Relative error [\%]','Interpreter','latex','FontSize',16)
end

