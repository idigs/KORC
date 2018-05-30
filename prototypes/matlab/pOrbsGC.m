function ST = pOrbsGC(pathToBField,fileType,ND,res,timeStepParams,tracerParams,xo,vo,opt)
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
% ST = pOrbsGC('','','2D',[],[1E4,2.5E-2,100,500],[-1,1],[1.5,0.0,0.3],[15E6,10],false);
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
% ST = pOrbs(name,'SIESTA','3D',[100,79,149],[5E1,1E-2,20,1E-3],[-1,1],[1.5,0,-0.5],[40E6,10]);
% USING RAW FILES
% ST = pOrbs('jfit_165365_1400.mat','RAW','2D',[],[1E5,1E-2,10],[-1,1],[1.82,0,-0.4928],[0.99,70]);
% name = 'xpand_iter3D_sc4_bet015_I87_hi_acc';
% ST = pOrbs(name,'XPANDER','3D',[150,100,150],[1E4,1E-2,10,1E-3],[-1,1],[6.5,0,0],[1E6,5],true);

narginchk(8,9);

close all

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
ST.F = struct;
ST.F.init = true;
ST.F.norm = false;
ST.F = analyticalFields(ST,[1,1,1],true);
ST.Bo = ST.F.Bo;

% Particle's parameters and nitial position and velocity of tracer, in SI units
ST.params.q = tracerParams(1)*ST.params.qe; %alpha-particle
ST.params.m = tracerParams(2)*ST.params.me; % alpha-particle

ST.params.Xo = xo; % Initial position
ST.Eo = vo(1)*ST.params.qe + ST.params.m*ST.params.c^2;
vo(1) = sqrt(1 - (ST.params.m*ST.params.c^2/ST.Eo)^2);
ST.params.vo_params = vo; % vo_params = [velocity magnitude, pitch angle]
[ST.params.vo, ~, ~] = initializeVelocity(ST);

ST.params.wc = sqrt(1 - sum(ST.params.vo.^2)/ST.params.c^2)*abs(ST.params.q)*ST.Bo/ST.params.m;
ST.params.Tc = 2*pi/ST.params.wc;

ST.coll.nimpurities = 1;
ST.coll.Te = 1.0*ST.params.qe; % Background electron temperature in Joules
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

if ST.opt
    disp(['Cyclotron frequency ' num2str(ST.params.wc) ' Hz'])
end

% Normalisation parameters
ST.norm.v = ST.params.c;
ST.norm.q = abs(ST.params.q);
ST.norm.m = ST.params.m;
ST.norm.B = ST.Bo;
ST.norm.E = ST.norm.v*ST.norm.B;
ST.norm.wc = ST.norm.q*ST.norm.B/ST.norm.m;
ST.norm.t = 1/ST.norm.wc;
ST.norm.l = ST.norm.v*ST.norm.t;
ST.norm.Eo = ST.norm.m*ST.norm.v^2;
ST.norm.T = ST.norm.m*ST.norm.v^2;
ST.norm.p = ST.norm.m*ST.norm.v;
ST.norm.n = 1/ST.norm.l^3;
ST.norm.mu = ST.norm.m*ST.norm.v^2/ST.norm.B;
% Normalisation parameters

% Numerical parameters
ST.params.numIt = timeStepParams(1);
ST.params.dt = timeStepParams(2)*(2*pi/ST.params.wc);
ST.params.cadence = timeStepParams(3);
ST.params.subcycling = timeStepParams(4);
ST.params.inds = 1:ST.params.cadence:ST.params.numIt;
ST.params.numSnapshots = numel(ST.params.inds);

ST.time = zeros(ST.params.numSnapshots,1);
for ii=1:ST.params.numSnapshots
    ST.time(ii) = (ii-1)*ST.params.dt*ST.params.cadence;
end

ST.time = ST.time;%/(2*pi/ST.params.wc);

ST.cOp = initializeCollisionOperators(ST);

% Normalization of fields
ST.F.norm = true;
ST.F = analyticalFields(ST,[ST.F.Ro,0,0],false);
% Normalization of fields

ST.PP = particlePusherLeapfrog(ST,false,false);

ST.GC = particlePusherGC(ST,false,false);

% Particle pusher using matlab ODEs solvers
% ST.PP = particlePusherMatlab(ST); % (Un)comment this line as required
% Particle pusher using matlab ODEs solvers

% if ST.opt
%     PoincarePlots(ST);
% end

% ST.PP.invariant = invariants(ST);

munlock

end

% INITIALIZATION

function [vo,vpar,vperp] = initializeVelocity(ST,Xo,V,pitchAngle)
% ST.params.vo_params = [vo, pitch_angle]
x = [1,0,0];
y = [0,1,0];
z = [0,0,1];

if nargin == 1
    Xo = ST.params.Xo;
    V = ST.params.vo_params(1)*ST.params.c; % magnitude of initial velocity
    pitchAngle = pi*ST.params.vo_params(2)/180;
end

rng('shuffle')
angle = 2*pi*rand;

vpar = V*cos(pitchAngle);
vperp = V*sin(pitchAngle);

v1 = V*cos(pitchAngle);
v2 = V*sin(pitchAngle)*cos(angle);
v3 = V*sin(pitchAngle)*sin(angle);

[b1,b2,b3] = unitVectors(ST,Xo);

vo = [0,0,0];
vo(1) = v1*(b1*x') + v2*(b2*x') + v3*(b3*x');
vo(2) = v1*(b1*y') + v2*(b2*y') + v3*(b3*y');
vo(3) = v1*(b1*z') + v2*(b2*z') + v3*(b3*z');
end

function [b1,b2,b3] = unitVectors(ST,Xo)
% initial condition of an electron drifting parallel to the local magnetic
% field.

F = analyticalFields(ST,Xo,true);
b = F.B/sqrt(F.B*F.B');

b1 = b;

b2 = cross(b,[0,0,1]);
b2 = b2/sqrt(b2*b2');

b3 = cross(b1,b2);
b3 = b3/sqrt(b3*b3');
end

function F = analyticalFields(ST,X,coords,q,mu,ppar)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.
% coordsys = true (Cartesian), coordsys = false (Cylindrical)

narginchk(2,6);

if coords
    R = sqrt(X(1)^2 + X(2)^2);
    phi = atan2(X(2),X(1));
    if (phi < 0); phi = phi + 2*pi; end
    Z = X(3);
else
    R = X(1);
    phi = X(2);
    Z = X(3);
end

if ST.F.init == true
    % Parameters of the analytical magnetic field
    % DIII-D
    jpb = -1;
    Ro = 1.5; % Major radius in meters.
    Zo = 0;
    Bo = 2.19;
    a = 0.5;% Minor radius in meters.
    qa = 5.0; % Safety factor at the separatrix (r=a)
    qo = 1.0; % Safety factor at the magnetic axis.
    lamb = a/sqrt(qa/qo - 1);
    Bpo = lamb*Bo/(qo*Ro);
    
    % Parameters of the analytical electric field
    Eo = 0.0;
    % Parameters of the analytical magnetic field
  
    F = ST.F;
    F.Bo = Bo;
    F.a = a;% 0.6;% Minor radius in meters.
    F.Ro = Ro; % Major radius in meters.
    F.Zo = Zo;
    F.qa = qa; % Safety factor at the separatrix (r=a)
    F.qo = qo;
    F.lamb = lamb;
    F.Bpo = Bpo;
    F.jpb = jpb;
    
    F.Eo = Eo;
    
    F.init = false;
    
    if (Bo > 0)
        disp('The toroidal magnetic field rotates clockwise')
    else
        disp('The toroidal magnetic field rotates counterclockwise')
    end
    if (jpb > 0)
        disp('The plasma current is parallel to the toroidal magnetic field')
    else
        disp('The plasma current is anti-parallel to the toroidal magnetic field')
    end
    disp(['q-factor at magnetic axis: ' num2str(qo)])
elseif ST.F.norm == true
    F = ST.F;
    
    F.Bo = F.Bo/ST.norm.B;
    F.a = F.a/ST.norm.l;
    F.Ro = F.Ro/ST.norm.l;
    F.Zo = F.Zo/ST.norm.l;
    F.lamb = F.lamb/ST.norm.l;
    F.Bpo = F.Bpo/ST.norm.B;
    
    F.Eo = F.Eo/ST.norm.E;
    
    F.norm = false;
else
    F = ST.F;
    
    % Magnetic field
    qprof = (ST.F.qo/ST.F.lamb^2)*(ST.F.lamb^2 + (R-ST.F.Ro)^2 + (Z-ST.F.Zo)^2);
    BR = -ST.F.jpb*ST.F.Bo*(Z-ST.F.Zo)/(qprof*R);
    Bphi = -ST.F.Ro*ST.F.Bo/R;
    BZ = ST.F.jpb*ST.F.Bo*(R-ST.F.Ro)/(qprof*R);
    
    Bmag = sqrt(BR^2 + Bphi^2 + BZ^2);
    % Magnetic field
    
    % Electric field
    ER = 0;
    Ephi = -ST.F.Ro*ST.F.Eo/R;
    EZ = 0;
    % Electric field
    
    
    if coords
        Bx = BR*cos(phi) - Bphi*sin(phi);
        By = BR*sin(phi) + Bphi*cos(phi);
        Bz = BZ;
        
        F.B = [Bx,By,Bz];
        
        Ex = ER*cos(phi) - Ephi*sin(phi);
        Ey = ER*sin(phi) + Ephi*cos(phi);
        Ez = EZ;
        
        F.E = [Ex,Ey,Ez];
    else
        F.B = [BR,Bphi,BZ];
        F.E = [ER,Ephi,EZ];
    end
    
    if (nargin > 3)
        G = sqrt( 1 + 2*mu*sqrt(F.B*F.B') + ppar^2 ); % This definition of G is only valid when ST.norm.v = c and ST.norm.m = ST.params.m !!
        
        dqprofdR = 2*ST.F.qo*(R-ST.F.Ro)/ST.F.lamb^2;
        dqprofdZ = 2*ST.F.qo*(Z-ST.F.Zo)/ST.F.lamb^2;
        
        dBRdR = ST.F.jpb*ST.F.Bo*(Z-ST.F.Zo)*(R*dqprofdR + qprof)/(qprof*R)^2;
        dBphidR = ST.F.Bo*ST.F.Ro/R^2;
        dBZdR = ST.F.jpb*ST.F.Bo*(qprof*R - (R-ST.F.Ro)*(R*dqprofdR + qprof))/(qprof*R)^2;
        
        dBdR = (BR*dBRdR + Bphi*dBphidR + BZ*dBZdR)/Bmag;

        dBRdZ = -ST.F.jpb*ST.F.Bo*(qprof - (Z-ST.F.Zo)*dqprofdZ)/(R*qprof^2);
        dBphidZ = 0;
        dBZdZ = -ST.F.jpb*ST.F.Bo*(R-ST.F.Ro)*dqprofdZ/(R*qprof^2);
        
        dBdZ = (BR*dBRdZ + Bphi*dBphidZ + BZ*dBZdZ)/Bmag;
        
        rotbR = Bphi*dBdZ/Bmag^2;
        rotbphi = (dBRdZ - dBZdR)/Bmag - (BR*dBdZ - BZ*dBdR)/Bmag^2;
        rotbZ = -Bphi*dBdR/Bmag^2;
        
        gradBR = dBdR;
        gradBphi = 0;
        gradBZ = dBdZ;
        
        HR = BR + ppar*rotbR/q;
        Hphi = Bphi + ppar*rotbphi/q;
        HZ = BZ + ppar*rotbZ/q;
        
        F.Bstar = [HR,Hphi,HZ];
        
        DR = ER - 2*mu*gradBR/(q*G);
        Dphi = Ephi - 2*mu*gradBphi/(q*G);
        DZ = EZ - 2*mu*gradBZ/(q*G);
        
        F.Estar = [DR,Dphi,DZ];
    end
end

end

% Monte-Calor collision operators

function cOp = initializeCollisionOperators(ST)
cOp = struct;

cOp.Te = (1.0E3)*ST.params.qe; % Background electron temperature in Joules
cOp.Ti = cOp.Te; % Background ion temperature in Joules
cOp.ne = 1.0E20; % Background electron density in 1/m^3
cOp.Zeff = 10.0; % Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
cOp.rD = ...
    sqrt( ST.params.ep*cOp.Te/(cOp.ne*ST.params.qe^2*(1 + cOp.Zeff*cOp.Te/cOp.Ti)) );
cOp.re = ST.params.qe^2/( 4*pi*ST.params.ep*ST.params.me*ST.params.c^2 );
cOp.Clog = 25.3 - 1.15*log10(1E-6*cOp.ne) + 2.3*log10(cOp.Te/ST.params.qe);
cOp.VTe = sqrt(2*cOp.Te/ST.params.me);
cOp.delta = cOp.VTe/ST.params.c;
cOp.Gamma = cOp.ne*ST.params.qe^4*cOp.Clog/(4*pi*ST.params.ep^2);
cOp.Tau = ST.params.me^2*ST.params.c^3/cOp.Gamma;
cOp.Tauc = ST.params.me^2*cOp.VTe^3/cOp.Gamma;
cOp.Ec = cOp.ne*ST.params.qe^3*cOp.Clog/(4*pi*ST.params.ep^2*ST.params.me*ST.params.c^2);

disp(['Relativistic collisional time: ' num2str(cOp.Tau) ' s'])
disp(['Thermal collisional time: ' num2str(cOp.Tauc) ' s'])

F = analyticalFields(ST,[ST.F.Ro,0,0],true);
Ef = F.E;
cOp.Vc = cOp.VTe*sqrt(0.5*cOp.Ec/sqrt(Ef*Ef'));

g = linspace(1,1.03,200);
energy = g*ST.params.m*ST.params.c^2;%*linspace(6E5,50E6,200)*ST.params.qe;
u = ST.params.c*sqrt(1 - (ST.params.m*ST.params.c^2./energy).^2);

% Normalization
cOp.Te = cOp.Te/(ST.norm.m*ST.norm.v^2);
cOp.Ti = cOp.Ti/(ST.norm.m*ST.norm.v^2);
cOp.ne = cOp.ne/ST.norm.l^3;
cOp.VTe = cOp.VTe/ST.norm.v;

cOp.Gamma = cOp.Gamma*(ST.norm.t/(ST.norm.m^2*ST.norm.v^3));
cOp.Tau = cOp.Tau/ST.norm.t;
cOp.Vc = cOp.Vc/ST.norm.v;

u = u/ST.norm.v;
% cOp.dt = cOp.dt/ST.norm.t;
% Normalization

cOp.x = @(v) v/cOp.VTe;

cOp.psi = @(v) 0.5*( erf(cOp.x(v)) - ...
    2*cOp.x(v).*exp(-cOp.x(v).^2)/sqrt(pi) )./cOp.x(v).^2;

cOp.CA = @(v) cOp.Gamma*cOp.psi(v)./v;

cOp.CF = @(v) cOp.Gamma*cOp.psi(v)/cOp.Te;

cOp.CB = @(v) (0.5*cOp.Gamma./v).*( cOp.Zeff + ...
    erf(cOp.x(v)) - cOp.psi(v) + 0.5*cOp.delta^4*cOp.x(v).^2 );

cOp.fun = @(v) 2*(1./cOp.x(v) + cOp.x(v)).*exp(-cOp.x(v).^2)/sqrt(pi) - ...
    erf(cOp.x(v))./cOp.x(v).^2 - cOp.psi(v);

cOp.timeStep = @(v) 2*cOp.CA(v)./cOp.CF(v).^2;

vo = ST.params.vo_params(1);

if cOp.timeStep(vo) > (cOp.Tau*ST.params.subcycling)
    cOp.dt = 10*cOp.timeStep(vo);
else
    cOp.dt = cOp.Tau*ST.params.subcycling;
end
cOp.subcyclingIter = floor(ST.norm.t*cOp.dt/ST.params.dt);

cOp.ratio = @(v) sqrt(cOp.dt)*cOp.CF(v)./sqrt(2*cOp.CA(v));


if ST.opt
    disp(['Collisional time ' num2str(cOp.Tau*ST.norm.t) ' (s)'])
    disp(['Time step of subcycling ' num2str(cOp.dt*ST.norm.t) ' (s)'])
    disp(['Calculate collisions each: ' num2str(cOp.subcyclingIter) ' leap-frog iterations'])
end

if ST.opt
    energy = 1E-6*energy/ST.params.qe;
    
%     figure;
%     subplot(3,2,1);
%     plot(energy,cOp.CA(u));
%     ylabel('$C_A$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
%     xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
%     subplot(3,2,2);
%     plot(energy,cOp.CB(u));
%     ylabel('$C_B$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
%     xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
%     subplot(3,2,3);
%     plot(energy,cOp.CF(u));
%     ylabel('$C_F$ ($e B_0/m_e^2 c $)','Interpreter','latex')
%     xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
%     subplot(3,2,4);
%     plot(energy,cOp.ratio(u));
%     ylabel('$C_Fdt/\sqrt{2C_A dt}$','Interpreter','latex')
%     xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
%     subplot(3,2,5);
%     plot(energy,cOp.timeStep(u)*ST.norm.t);
%     ylabel('$dt$ (s)','Interpreter','latex')
%     xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')   
%     subplot(3,2,6);
%     plot(energy,cOp.timeStep(u)/(ST.params.dt/ST.norm.t));
%     ylabel('Iterations','Interpreter','latex')
%     xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
    
    figure;
    subplot(3,1,1);
    plot(energy,cOp.CA(u));
    ylabel('$C_A$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
    xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
    subplot(3,1,2);
    plot(energy,cOp.CB(u));
    ylabel('$C_B$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
    xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
    subplot(3,1,3);
    plot(energy,cOp.CF(u));
    ylabel('$C_F$ ($e B_0/m_e^2 c $)','Interpreter','latex')
    xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
end
end

function [U,Wcoll] = collisionOperator(ST,X,V,dt)
% This function calculates the random kicks due to collisions
U = zeros(1,3);

v = sqrt(V*V');
gammap = 1/sqrt(1 - v^2);
p = gammap*v;

[b1,b2,b3] = unitVectors(ST,X);

U1 = gammap*V*b1';
U2 = gammap*V*b2';
U3 = gammap*V*b3';

xi = U1/(gammap*v);
pitch = acos(xi);

dWp = sqrt(dt)*sqrt(-2*log(1-rand))*cos(rand);
dWtheta = random('uniform',0,pi,1)*sqrt(dt);
dWxi = random('uniform',0,1,1)*sqrt(dt);
dWphi = random('uniform',0,2*pi,1)*sqrt(dt);

CA = ST.cOp.CA(v);
CB = ST.cOp.CB(v);
CF = ST.cOp.CF(v);

dp = ( -CF + 2*CA/p + ST.cOp.Gamma*ST.cOp.fun(v)/(p^2*sqrt(1 + p^2)) )*dt + ...
    sqrt(2*CA)*dWp;

dxi = -2*xi*CB*dt/p^2 - sqrt(2*CB*(1 - xi^2))*dWxi/p;

dphi = sqrt(2*CB)*dWphi/(p*sqrt(1 - xi^2));
if ~isfinite(dphi)
    dphi = 0;
end

dpitch = acos(xi + dxi);

x = [1,0,0];
y = [0,1,0];
z = [0,0,1];

dU1 = dp*cos(dpitch);
dU2 = dp*sin(dpitch)*cos(dphi);
dU3 = dp*sin(dpitch)*sin(dphi);


disp(['Collisions: ' num2str(U1) ' , ' num2str(U2) ' , ' num2str(U3)])

U(1) = (U1+dU1)*(b1*x') + (U2+dU2)*(b2*x') + (U3+dU3)*(b3*x');
U(2) = (U1+dU1)*(b1*y') + (U2+dU2)*(b2*y') + (U3+dU3)*(b3*y');
U(3) = (U1+dU1)*(b1*z') + (U2+dU2)*(b2*z') + (U3+dU3)*(b3*z');

Wcoll = (gammap - sqrt(1 + U*U'))/dt;
end

function [Uc,Wcoll] = CoulombCollision(ST,U,dt)
% This function calculates the random kicks due to collisions
Uc = zeros(1,3);
x = [1,0,0];
y = [0,1,0];
z = [0,0,1];

g = sqrt( 1 + U*U' );
V = U/g;
v = sqrt(V*V');

v1 = V/v;

v2 = cross(v1,y);
v2 = v2/sqrt(v2*v2');

v3 = cross(v1,v2);
v3 = v3/sqrt(v3*v3');

U1 = U*v1';
U2 = U*v2';
U3 = U*v3';

dW = random('norm',0,sqrt(dt),[3,1]);

CA = ST.cOp.CA(v);
CB = ST.cOp.CB(v);
CF = ST.cOp.CF(v);

dU = [-2.0*CF*dt + sqrt(2*CA)*dW(1);
    sqrt(2*CB)*dW(2);
    sqrt(2*CB)*dW(3)];

% disp(['Collisions: ' num2str(U1) ' , ' num2str(U2) ' , ' num2str(U3)])

Uc(1) = (U1+dU(1))*(v1*x') + (U2+dU(2))*(v2*x') + (U3+dU(3))*(v3*x');
Uc(2) = (U1+dU(1))*(v1*y') + (U2+dU(2))*(v2*y') + (U3+dU(3))*(v3*y');
Uc(3) = (U1+dU(1))*(v1*z') + (U2+dU(2))*(v2*z') + (U3+dU(3))*(v3*z');

Wcoll = (g - sqrt(1 + Uc*Uc'))/dt;
end

% LEAP-FROG PARTICLE PUSHER

function PP = particlePusherLeapfrog(ST,coll,rad)
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
kapp = zeros(1,ST.params.numSnapshots); % Approximation of the curvature
T = zeros(1,ST.params.numSnapshots); % Torsion
vpar = zeros(1,ST.params.numSnapshots); % parallel velocity
vperp = zeros(1,ST.params.numSnapshots); % perpendicular velocity
ppar = zeros(1,ST.params.numSnapshots); % parallel momentum
pperp = zeros(1,ST.params.numSnapshots); % perpendicular momentum
mu = zeros(1,ST.params.numSnapshots); % instantaneous magnetic moment
EK = zeros(1,ST.params.numSnapshots); % kinetic energy

fL = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
f1 = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
f2 = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
f3 = zeros(1,ST.params.numSnapshots); % Synchroton radiated power

WL = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
WR = zeros(3,ST.params.numSnapshots); % Synchroton radiated power
Psyn = zeros(1,ST.params.numSnapshots); % Synchroton radiated power
Wcoll = zeros(1,ST.params.numSnapshots); % Synchroton radiated power

% Normalization
X(1,:) = ST.params.Xo/ST.norm.l;
v(1,:) = ST.params.vo/ST.norm.v;
q = ST.params.q/ST.norm.q;
m = ST.params.m/ST.norm.m;

F = analyticalFields(ST,X(1,:),true);
B = F.B;
E = F.E;
% dt = ST.params.dt*ST.norm.wc;
dt = ST.params.dt/ST.norm.t;
ep = ST.params.ep*(ST.norm.m^2*ST.norm.v^3/(ST.norm.q^3*ST.Bo));
Kc = ST.params.Kc/(ST.norm.m*ST.norm.l*ST.norm.v^2/ST.norm.q^2);
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
ppar(1) = gamma*m*vpar(1);
pperp(1) = gamma*m*vperp(1);
mu(1) = m*gamma^2*vperp(1)^2/(2*B_mag);

% Curvature and torsion
vmag = sqrt( sum(v(1,:).^2) );
aux =  cross(v(1,:),E) + v(1,:)*sum(v(1,:).*B) - B*vmag^2;
k(1) = abs(q)*sqrt( sum(aux.^2) )/(gamma*m*vmag^3);

pitch = atan2(vperp(1),vpar(1));
if (pitch < 0)
    pitch = pitch + 2*pi;
end
% kapp(1) = abs(q)*B_mag*sin(pitch)/sqrt(ppar(1)^2 + pperp(1)^2);
kapp(1) = abs(q)*sin(pitch)/sqrt(ppar(1)^2 + pperp(1)^2);% Bo = 1 due to normalization
curv = k(1);
% Curvature and torsion

% % % Radiation damping force % % %

if rad
    Psyn(1) = -(2/3)*( Kc*q^2*gamma^4*vmag^4*curv^2 );
    
    F1 = 0;
    F2 = ( q^3/(6*pi*ep*m^2) )*( (E*v(1,:)')*E + cross(E,B) +...
        cross(B,cross(B,v(1,:))) );
    vec = E + cross(v(1,:),B);
    F3 = ( gamma^2*q^3/(6*pi*ep*m^2) )*( (E*v(1,:)')^2 - vec*vec' )*v(1,:);
    
    WL(1) = q*(E*v(1,:)');
    WR(1,1) = ( q^3/(6*pi*ep*m) )*(gamma*DEDt*v(1,:)');
    WR(2,1) = ( q^4/(6*pi*ep*m^2) )*( E*E' + cross(v(1,:),B)*E' );
    WR(3,1) = ( q^4/(6*pi*ep*m^2) )*( gamma^2*( (E*v(1,:)')^2 - vec*vec' ) );
    
    fL(1) = abs(q)*sqrt( (F1+F2+F3)*(F1+F2+F3)' );
    f1(1) = abs(q)*sqrt( F1*F1' );
    f2(1) = abs(q)*sqrt( F2*F2' );
    f3(1) = abs(q)*sqrt( F3*F3' );
end
% % % % % % % % % % % % % % % % % %

% initial half-step for position
DT = 0.5*dt;
V = v(1,:);
U = u(1,:);
XX = X(1,:) + DT*V; % Advance position
% initial half-step for position

a = q*dt/m;
dummyWcoll = 0;

for ii=2:ST.params.numSnapshots
    for jj=1:ST.params.cadence
        
        F = analyticalFields(ST,XX,true);
        B = F.B;
        E = F.E;
        
        U_R = U;
        U_L = U;
        
        % % % Leap-frog scheme for Lorentz force % % %
        U_hs = U_L + 0.5*a*(E + cross(V,B)); % Half step for velocity
        
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
        % % % Leap-frog scheme for Lorentz force % % %
        
        if rad
            % % % Leap-frog scheme for the radiation damping force % % %
            U_eff = 0.5*(U_L + U);
            gamma_eff = sqrt(1 + U_eff*U_eff');
            V_eff = U_eff/gamma_eff;
            
            F1 = 0;
            F2 = ( q^3/(6*pi*ep*m^2) )*( (E*V_eff')*E + cross(E,B) +...
                cross(B,cross(B,V_eff)) );
            vec = E + cross(V_eff,B);
            F3 = ( gamma_eff^2*q^3/(6*pi*ep*m^2) )*( (E*V_eff')^2 - vec*vec' )*V_eff;
            
            U_R = U_R + a*( F1 + F2 + F3 );
            U = U_L + U_R - U;
            % % % Leap-frog scheme for the radiation damping force % % %
        else
            U = U_L;
        end
        
        
        if coll && (mod((ii-1)*ST.params.cadence + jj,ST.cOp.subcyclingIter) == 0)
            % %  % Collisions %  % %
            [U,dummyWcoll] = CoulombCollision(ST,U,dt);
            % %  % Collisions %  % %
        end
        
        gamma = sqrt( 1 + U*U' ); % Comment or uncomment
        V = U/gamma;
        
        zeta_previous = atan2(XX(2),XX(1));
        if zeta_previous < 0
            zeta_previous = zeta_previous + 2*pi;
        end
        
        XX = XX + dt*V; % Advance position
        
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
    ppar(ii) = gamma*m*vpar(ii);
    pperp(ii) = gamma*m*vperp(ii);
    mu(ii) = m*gamma^2*vperp(ii)^2/(2*B_mag);
    
    k(ii) = curv;
    
    pitch = atan2(vperp(ii),vpar(ii));
    if (pitch < 0)
        pitch = pitch + 2*pi;
    end
    
    kapp(ii) = abs(q)*sin(pitch)/sqrt(ppar(ii)^2 + pperp(ii)^2);% Bo = 1 due to normalization
    
    if rad
        Psyn(ii) = -(2/3)*( Kc*q^2*gamma^4*vmag^4*curv^2 );
        
        fL(ii) = abs(q)*sqrt( (F1+F2+F3)*(F1+F2+F3)' );
        f1(ii) = abs(q)*sqrt( F1*F1' );
        f2(ii) = abs(q)*sqrt( F2*F2' );
        f3(ii) = abs(q)*sqrt( F3*F3' );
        
        WL(ii) = q*(E*V_eff');
        WR(1,ii) = ( q^3/(6*pi*ep*m) )*(gamma_eff*DEDt*V_eff');
        WR(2,ii) = ( q^4/(6*pi*ep*m^2) )*( E*E' + cross(V_eff,B)*E' );
        WR(3,ii) = ( q^4/(6*pi*ep*m^2) )*( gamma_eff^2*( (E*V_eff')^2 - vec*vec' ) );
    end
    
    if coll
        Wcoll(ii) = dummyWcoll;
    end
end


if ST.opt
    time = ST.time;%/(2*pi/ST.params.wc);
    
    % Relative error in energy conservation
    ERR_EK = 100*(EK - EK(1))./EK(1);
    ERR_mu = 100*(mu - mu(1))./mu(1);
    % Relative error in energy conservation
    
    x = linspace(min(time),max(time),10);
    y = zeros(1,10);
    
    h = figure;
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
    tmp = floor((ST.Eo/ST.params.qe)/1E6);
    saveas(h,['energy_Eo_' num2str(tmp) '_po_' num2str(ST.params.vo_params(2))],'fig')
       
    
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
    box on; grid on
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
    
    figure
    subplot(2,1,1)
    plot(time, ppar)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$p_\parallel$ [$mc$]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(2,1,2)
    plot(time, pperp)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$p_\perp$ [$mc$]','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
end

% Return position and velocity with SI units
PP.X = X*ST.norm.l;
PP.R = R*ST.norm.l;
PP.v = v*ST.norm.v;

PP.vpar = vpar;
PP.vperp = vperp;

PP.k = k/ST.norm.l; % curvature
PP.kapp = kapp/ST.norm.l; % curvature
PP.T = T/ST.norm.l; % curvature

PP.gamma = EK;
PP.mu = mu;

PP.Psyn = sum(WR,1)*ST.norm.m*ST.norm.v^3/ST.norm.l;
Psyn = Psyn*ST.norm.m*ST.norm.v^3/ST.norm.l;


POINCARE.R = POINCARE.R*ST.norm.l;
POINCARE.Z = POINCARE.Z*ST.norm.l;
POINCARE.k = POINCARE.k/ST.norm.l;
POINCARE.T = POINCARE.T/ST.norm.l;

PP.POINCARE = POINCARE;


if ST.opt
    figure
    subplot(2,1,1)
    plot(time,PP.k,'k',time,PP.kapp,'r')
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
    subplot(4,1,1)
    plot(time,fL,'k')
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|F_L|$','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(4,1,2)
    plot(time, f1,'r')
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|f_1|$','Interpreter','latex','FontSize',16)
    subplot(4,1,3)
    plot(time, f2,'b')
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|f_2|$','Interpreter','latex','FontSize',16)
    subplot(4,1,4)
    plot(time,f3,'m')
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|f_3|$','Interpreter','latex','FontSize',16)
    
    figure
    plot(time,WR(1,:),'r',time,WR(2,:),'b',time,WR(3,:),'k')
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
    
    figure
    subplot(4,1,1)
    plot(time, WL)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|W_L|$','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(4,1,2)
    plot(time,abs(PP.Psyn),'k',time,abs(Psyn),'r')
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
    subplot(4,1,3)
    plot(time, Wcoll)
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|W_{coll}|$','Interpreter','latex','FontSize',16)
    subplot(4,1,4)
    plot(time, abs(PP.Psyn./WL))
    box on
    grid on
    xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
    ylabel('$|W_R/W_L|$','Interpreter','latex','FontSize',16)
    disp(['Mean value: ' num2str(mean(sum(WR,1)./WL))])
end

    figure;
    subplot(2,2,[1,2])
    plot3(PP.X(:,1),PP.X(:,2),PP.X(:,3),'b',PP.X(1,1),PP.X(1,2),PP.X(1,3),'ko','MarkerFaceColor',[0,0,0])
%     hold on
%     plot3(PP.R(:,1),PP.R(:,2),PP.R(:,3),'r')
%     hold off
    axis equal
    box on
    legend({'Full orbit','$\vec{x}_0$','Guiding center'},'Interpreter','latex')
    xlabel('X','Interpreter','latex','FontSize',16)
    ylabel('Y','Interpreter','latex','FontSize',16)
    zlabel('Z','Interpreter','latex','FontSize',16)
    title(PP.method,'Interpreter','latex','FontSize',16)
    subplot(2,2,3)
    
    subplot(2,2,4)
    R = sqrt(PP.X(:,1).^2 + PP.X(:,2).^2);
    Z = PP.X(:,3);
    plot(R,Z,'b.','MarkerSize',2)
    axis equal;
    xlabel('R','Interpreter','latex','FontSize',16)
    ylabel('Z','Interpreter','latex','FontSize',16)

if ST.opt
    disp('Particle pusher: done!')
end
end

% Guiding-center

function GC = particlePusherGC(ST,coll,rad)
% Relativistic particle pusher for the guiding-center
GC = struct;

q = ST.params.q;
m = ST.params.m;

Xo = zeros(1,3); %(R,phi,Z)
Xo(1) = sqrt(ST.params.Xo(1)^2 + ST.params.Xo(2)^2);
Xo(2) = atan2(ST.params.Xo(2),ST.params.Xo(1));
if (Xo(2) < 0)
    Xo(2) = Xo(2) + 2*pi;
end
Xo(3) = ST.params.Xo(3);

Go = 1/sqrt(1-ST.params.vo*ST.params.vo'/ST.params.c^2);
po = ST.params.m*ST.params.c*sqrt(Go^2 - 1);
pparo = po*cos(deg2rad(ST.params.vo_params(2)));
ppero = sqrt(po^2 - pparo^2);
vo = sqrt(ST.params.vo*ST.params.vo');
vparo = vo*cos(deg2rad(ST.params.vo_params(2)));
vpero = vo*sin(deg2rad(ST.params.vo_params(2)));

% % % Normalization
Xo(1) = Xo(1)/ST.norm.l;
Xo(3) = Xo(3)/ST.norm.l;
po = po/ST.norm.p;
pparo = pparo/ST.norm.p;
ppero = ppero/ST.norm.p;
vparo = vparo/ST.norm.v;
vpero = vpero/ST.norm.v;
q = q/ST.norm.q;
m = m/ST.norm.m;
% % % Normalization

F = analyticalFields(ST,Xo,false);
B = F.B;
E = F.E;

mu = ppero^2/(2*m*sqrt((B*B')));

F = analyticalFields(ST,Xo,false,q,mu,pparo);
B = F.B;
E = F.E;
Bstar = F.Bstar;
Estar = F.Estar;
b = B/sqrt(B*B');
Bspar = dot(b,Bstar);

options = odeset('RelTol',1E-6,'AbsTol',1E-10,'Stats','on')%,'OutputFcn',@odeplot)
% options = odeset('RelTol',1E-3,'AbsTol',1E-6,'Stats','on')%,'OutputFcn',@odeplot)

DT = 2*pi/ST.params.wc;
% DT = 2*DT;
tfinal = ST.params.dt*ST.params.numIt;
numIt = floor(tfinal/DT);
tspan = linspace(0,tfinal,numIt)/ST.norm.t;
y0 = [Xo,pparo];

[~,y] = ode45(@(t,y) eqnsMotionGC(t,y,q,m,mu,ST),tspan,y0,options);  % NONSTIFF PROBLEM

R = y(:,1)*ST.norm.l;
phi = y(:,2);%mod(y(:,2),2*pi);
Z = y(:,3)*ST.norm.l;

X = zeros(numel(R),3);

X(:,1) = R.*cos(phi);
X(:,2) = R.*sin(phi);
X(:,3) = Z;

figure;
subplot(2,2,[1 2])
plot3(X(:,1),X(:,2),X(:,3),'Color',[0.47,0.67,0.19])
hold on;plot3(X(1,1),X(1,2),X(1,3),'ko','MarkerFaceColor',[0,0,0]);hold off
axis equal
box on
legend({'Guiding center','$\vec{x}_0$'},'Interpreter','latex')
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
subplot(2,2,3)

subplot(2,2,4)
plot(R,Z,'.','MarkerSize',2,'Color',[0.47,0.67,0.19])
axis equal;
xlabel('R','Interpreter','latex','FontSize',16)
ylabel('Z','Interpreter','latex','FontSize',16)
end

function dydt = eqnsMotionGC(t,y,q,m,mu,ST)
F = analyticalFields(ST,y(1:3),false,q,mu,y(4));
B = F.B;
b = B/sqrt(B*B');
Bs = F.Bstar;
Es = F.Estar;
% Bs = F.B;
% Es = F.E;
Bspar = dot(b,Bs);
G = sqrt( 1 + 2*mu*sqrt(B*B') + y(4)^2 ); % This definition of G is only valid when ST.norm.v = c and ST.norm.m = ST.params.m !!

dydt = zeros(4,1);

dydt(1) = y(4)*Bs(1)/(m*G*Bspar) + (Es(2)*b(3) - Es(3)*b(2))/Bspar;
dydt(2) = (y(4)*Bs(2)/(m*G*Bspar) + (Es(3)*b(1) - Es(1)*b(3))/Bspar)/y(1);
dydt(3) = y(4)*Bs(3)/(m*G*Bspar) + (Es(1)*b(2) - Es(2)*b(1))/Bspar;
dydt(4) = q*dot(Es,Bs)/Bspar;

end

% POST-PROCESSING

function PoincarePlots(ST)
if isfield(ST.PP,'X')
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
    try
        t = linspace(0,2*pi,100);
        x95 = ST.B.Ro + ST.B.a*cos(t);
        y95 = ST.B.a*sin(t);
        hold on
        plot(R,X(:,3),'k',x95,y95,'b',ST.B.Ro,0,'bx')
        if isfield(ST.PP,'R')
            plot(Rgc,Zgc,'g','LineWidth',2)
        end
        hold off
        axis equal
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        title('Poincare plot','Interpreter','latex','FontSize',16)
    catch
    end
end
end

function I = invariants(ST)

Ro = ST.F.Ro*ST.norm.l;
Bo = ST.F.Bo*ST.norm.B;
lambda = ST.F.lamb*ST.norm.l;
qo = ST.F.qo;

m = ST.params.m;
q = ST.params.q;

X = ST.PP.X;
gammap = ST.PP.gamma';
V = ST.PP.v;

% Toroidal coordinates
% r = radius, theta = poloidal angle, zeta = toroidal angle
r = sqrt( (sqrt(X(:,1).^2 + X(:,2).^2) - Ro).^2 + X(:,3).^2 );
eta = r/Ro;
theta = atan2(X(:,3),sqrt(X(:,1).^2 + X(:,2).^2) - Ro);
theta(theta < 0) = theta(theta < 0) + 2*pi;
zeta = atan2(X(:,1),X(:,2));
zeta(zeta < 0) = zeta(zeta < 0) + 2*pi;
% Toroidal coordinates


dzeta_dt = ...
    (X(:,2).*V(:,1) - X(:,1).*V(:,2))./( X(:,1).^2 + X(:,2).^2 );

wc = q*Bo./(m*gammap);

% T1 = dzeta_dt.*(1 + eta.*cos(theta)).^2;
% T2 = 0.5*wc.*eta.^2/qo;

T1 = dzeta_dt.*(1 + eta.*cos(theta)).^2;
T2 = lambda^2*wc.*log(1 + (r/lambda).^2)/(2*qo*Ro^2);


I = T1 - T2;
err = (I - I(1))/I(1);

figure;
subplot(2,1,1)
plot(ST.time,T1,'r',ST.time,T2,'b',ST.time,I,'k')
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Decomposition','Interpreter','latex','FontSize',16)
subplot(2,1,2)
plot(ST.time,err)
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('$\Delta p_\zeta$ (\%)','Interpreter','latex','FontSize',16)

end


