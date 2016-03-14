function ST = particleOrbits_new(pathToBField,ND,numIt,DT,cadence,xo,vo_params,opt)
% Here ro and vo are the initial position and velocity of the particle.
% The components of the initial velocity must be entered as fractions of
% the speed of light in vacuum.
% EXAMPLE:
% ST = particleOrbits_new('','2D',1E6,1E-2,10,[5,0,0],[0.05,45],false);

close all

ST = struct;
% Script parameters
ST.analytical = true; % true = analytical B field; false = xpander fields

ST.pathToBField = pathToBField; % Path to xpander fields
ST.ND = ND; % Dimensionality of xpander fields, 2D or 3D.
ST.opt = opt; % true = plot xpander fields, false = don't plot xpander fields

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
    B = analyticalB([1,1,1]);
    ST.Bo = B.Bo;
else
    ST.B = loadXpandField(ST);
    ST.Bo = ST.B.Bo;
end

% Initial position and velocity of tracer, in SI units
ST.params.Xo = xo; % Initial position
ST.params.vo_params = vo_params; % vo_params = [velocity magnitude, pitch angle]
[ST.params.vo, ST.params.vpar, ST.params.vperp] = initializeVelocity(ST);

% Particle's parameters
% ST.params.q = -ST.params.qe; % electrons
% ST.params.m = ST.params.me; % electrons
ST.params.q = 2*ST.params.qe; %alpha-particle
ST.params.m = 6.64424E-27; % alpha-particle
ST.params.wc = sqrt(1 - sum(ST.params.vo.^2)/ST.params.c^2)*abs(ST.params.q)*ST.Bo/ST.params.m;
ST.params.rL = ST.params.vperp/ST.params.wc; % Larmor radius of particle

disp(['Cyclotron frequency ' num2str(ST.params.wc) ' Hz'])
disp(['Larmor radius ' num2str(ST.params.rL) ' m'])

calculateDrifts(ST); % calculate characteristic velocities of drifts

% Normalisation parameters
ST.norm.q = abs(ST.params.q);
ST.norm.m = ST.params.m;
ST.norm.wc = ST.norm.q*ST.Bo/ST.norm.m;
ST.norm.l = ST.params.c/ST.norm.wc;
% Normalisation parameters

% Numerical parameters
ST.params.dt = DT*(2*pi/ST.params.wc);
ST.params.numIt = numIt;
ST.params.cadence = cadence;

ST.time = zeros(ST.params.numIt,1);
for ii=1:ST.params.numIt
    ST.time(ii) = (ii-1)*ST.params.dt;
end

% ST.PP = particlePusherLeapfrog(ST);
% PoincarePlots(ST.PP.X);

ST.PP = particlePusherMatlab(ST);
PoincarePlots(ST.PP.X);

DiegosInvariant(ST);

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

function D = calculateDrifts(ST)
% curvature drift velocity
Wpar = 0.5*ST.params.m*ST.params.vpar^2;
R = sqrt( sum(ST.params.Xo(1:2).^2) );
B = analyticalB(ST.params.Xo);
Bo = sqrt( sum(B.B.^2) );
Vcurv = 2*Wpar/( ST.params.q*Bo*R );
disp(['Curvature drift velocity: ' num2str(Vcurv) ' m/s'])
% curvature drift velocity
end



function [b,a] = unitVectors(ST)
% initial condition of an electron drifting parallel to the local magnetic
% field.

% tmp = weightedInterpolation(ST,ST.params.Xo);

if ST.analytical
    B = analyticalB(ST.params.Xo);
    b = B.B/sqrt(sum(B.B.^2));
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

function B = loadXpandField(ST)
% All quantities in SI units
B = struct;

if strcmp(ST.ND,'2D')
    data = load([ST.pathToBField 'ITER' ST.ND '.dat']);
    % X(R,Z,phi)
    B.NR = 150;
    B.NZ = 150;
    
    B.R = zeros(B.NR,B.NZ);
    B.Z = zeros(B.NR,B.NZ);
    
    B.B = zeros(B.NR,B.NZ); % magnitude
    B.BR = zeros(B.NR,B.NZ);
    B.BZ = zeros(B.NR,B.NZ);
    B.P = zeros(B.NR,B.NZ);
    
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
    B.B = sqrt(B.BR.^2 + B.BZ.^2 + B.Bphi.^2);
    
    B.Bo = max(max(max(B.B)));
    
    % geometry params
    B.grid.R = B.R(1,:);
    B.grid.Z = B.Z(:,1);
    
    B.grid.Rmin = min(min(B.R));
    B.grid.Rmax = max(max(B.R));
    
    B.grid.Zmin = min(min(B.Z));
    B.grid.Zmax = max(max(B.Z));
    
    B.grid.dR = mean(diff(B.grid.R));
    B.grid.dZ = mean(diff(B.grid.Z));
    % geometry params
    
    if ST.opt
        figure
        subplot(1,4,1)
        surfc(B.R,B.Z,B.BR,'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B_R$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(1,4,2)
        surfc(B.R,B.Z,B.BZ,'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B_Z$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(1,4,3)
        surfc(B.R,B.Z,B.Bphi,'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$B_\phi$ [T]','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        
        subplot(1,4,4)
        surfc(B.R,B.Z,B.P,'LineStyle','none')
        view([0,90])
        axis equal
        box on
        colorbar
        title('$P(R,Z)$','Interpreter','latex','FontSize',16)
        xlabel('R [m]','Interpreter','latex','FontSize',16)
        ylabel('Z [m]','Interpreter','latex','FontSize',16)
        colormap(jet)
        
    end
    
elseif strcmp(ST.ND,'3D')
    data = load([ST.pathToBField 'ITER' ST.ND '.dat']);
    % X(R,Z,phi)
    B.NR = 150;
    B.NZ = 150;
    B.Nphi = 100;
    
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

% testInterpolation(B);
B.SI = calculatescatteredInterpolant(ST,B);

end

function SI = calculatescatteredInterpolant(ST,B)
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
    
    SI.CELL = {SI.BR, SI.Bphi, SI.BZ};
    
elseif strcmp(ST.ND,'3D')
    R = reshape(B.R,[numel(B.R) 1]);
    Z = reshape(B.Z,[numel(B.Z) 1]);
    phi = reshape(B.phi,[numel(B.phi) 1]);
    
    
    DATA = reshape(B.BR,[numel(B.BR) 1]);
    SI.BR = scatteredInterpolant(R,Z,phi,DATA);
    clear DATA
    
    DATA = reshape(B.BZ,[numel(B.BZ) 1]);
    SI.BZ = scatteredInterpolant(R,Z,phi,DATA);
    clear DATA
    
    DATA = reshape(B.Bphi,[numel(B.Bphi) 1]);
    SI.Bphi = scatteredInterpolant(R,Z,phi,DATA);
    clear DATA
    
    SI.CELL = {SI.BR, SI.Bphi, SI.BZ};
    
else
    error('Use 2D or 3D');
end

disp('Scattered interpolant: done!')
end

function F = testInterpolation(FIELD)
% Interpolation test of the magnetic field components using Matlab
% scatteredInterpolant function.

R = reshape(FIELD.R(:,:,1),[numel(FIELD.R(:,:,1)) 1]);
Z = reshape(FIELD.Z(:,:,1),[numel(FIELD.Z(:,:,1)) 1]);
DATA = reshape(FIELD.BR(:,:,1),[numel(FIELD.BR(:,:,1)) 1]);

disp('Calculating scattered interpolant...')
SI = scatteredInterpolant(R,Z,DATA);
disp('Scattered interpolant: done!')

Rmin = min(R);
Rmax = max(R);
Zmin = min(Z);
Zmax = max(Z);

[RR,ZZ] = meshgrid(linspace(Rmin,Rmax,250),linspace(Zmin,Zmax,250));

F = SI(RR,ZZ);

figure
surfc(RR,ZZ,F,'LineStyle','none')
view([0,90])
axis equal
box on
colormap(jet)
colorbar
title('$B_R$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)

end

function BF = interpMagField(ST,X)

if strcmp(ST.ND,'2D')
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
    
elseif strcmp(ST.ND,'3D')
    
else
    error('Use 2D or 3D');
end

end

function F = weightedInterpolation(ST,X)
% The field 'ST.B' is interpolated following Garnett R. W. et al. (1992)

if strcmp(ST.ND,'2D')
    
    R = sqrt(X(1)^2 + X(2)^2) - ST.B.grid.Rmin;
    Z = X(3) + abs(ST.B.grid.Zmin);
    phi = atan2(X(2),X(1));
    phi(phi<0) = phi(phi<0) + 2*pi;
    
    iR = floor((R + 0.5*ST.B.grid.dR)/ST.B.grid.dR) + 1; % index in the radial direction
    iZ = floor((Z + 0.5*ST.B.grid.dZ)/ST.B.grid.dZ) + 1; % index in the Z direction
    
    R = R + ST.B.grid.Rmin;
    Z = Z - abs(ST.B.grid.Zmin);
    
    % if logic == true the particle is on the right (above) of the node at
    % iR (iZ), otherwishe, the particle is on the left (below).
    logR = ( R - ST.B.grid.R(iR) ) > 0;
    logZ = ( Z - ST.B.grid.Z(iZ) ) > 0;
    
    % Fit of the radial component
    vec = ST.B.grid.R;
    FIELD = ST.B.BR;
    
    a1 = invrtMat(vec,FIELD,iR,iZ,'R');
    if logR
        a2 = invrtMat(vec,FIELD,iR+1,iZ,'R');
    else
        a2 = invrtMat(vec,FIELD,iR-1,iZ,'R');
    end
    
    vec = ST.B.grid.Z;
    b1 = invrtMat(vec,FIELD,iR,iZ,'Z');
    if logZ
        b2 = invrtMat(vec,FIELD,iR,iZ+1,'Z');
    else
        b2 = invrtMat(vec,FIELD,iR,iZ-1,'Z');
    end
    
    dFdR = a1(2) + 2*a1(3)*ST.B.grid.R(iR);
    dFdZ = b1(2) + 2*b1(3)*ST.B.grid.Z(iZ);
    
    ddFddR = 2*a1(3);
    ddFddZ = 2*b1(3);
    
    FR1 = FIELD(iR,iZ) + dFdR*(R - ST.B.grid.R(iR)) + ...
        dFdZ*(Z - ST.B.grid.Z(iZ)) + 0.5*ddFddR*(R - ST.B.grid.R(iR))^2 + ...
        0.5*ddFddZ*(Z - ST.B.grid.Z(iZ))^2;
    
    dFdR = a2(2) + 2*a2(3)*ST.B.grid.R(iR);
    dFdZ = b2(2) + 2*b2(3)*ST.B.grid.Z(iZ);
    
    ddFddR = 2*a2(3);
    ddFddZ = 2*b2(3);
    
    FR2 = FIELD(iR,iZ) + dFdR*(R - ST.B.grid.R(iR)) + ...
        dFdZ*(Z - ST.B.grid.Z(iZ)) + 0.5*ddFddR*(R - ST.B.grid.R(iR))^2 + ...
        0.5*ddFddZ*(Z - ST.B.grid.Z(iZ))^2;
    
    FR = FR1 % Check this!!
    % Fit of the radial component
    
    
    disp('done')
    
    
    
elseif strcmp(ST.ND,'3D')
    % do something for the full 3D case.
else
    error('Use 2D or 3D');
end

end

function a = invrtMat(vec,F,iR,iZ,coord)
% parameters of interpolations
M = 5; % number of grid points used for interpolation
% parameters of interpolations

if strcmp(coord,'R')
    A = zeros(4,4);
    
    A11 = M;
    A12 = sum(vec(iR-2:iR+2));
    A14 = -0.5;
    A22 = sum(vec(iR-2:iR+2).^2);
    A23 = sum(vec(iR-2:iR+2).^3);
    A24 = -0.5*vec(iR);
    A33 = sum(vec(iR-2:iR+2).^4);
    A34 = -0.5*vec(iR)^2;
    
    A(1,:) = [A11, A12, A22, A14];
    A(2,:) = [A12, A22, A23, A24];
    A(3,:) = [A22, A23, A33, A34];
    A(4,:) = [-2*A14, -2*A24, -2*A34, 0];
    
    br = [sum( F(iR-2:iR+2,iZ) ); ...
        sum( F(iR-2:iR+2,iZ)'.*vec(iR-2:iR+2) ); ...
        sum( F(iR-2:iR+2,iZ)'.*(vec(iR-2:iR+2).^2) ); ...
        F(iR,iZ)];
    
    a = A\br; % coefficients of quadratic expansion along the R-axis
elseif strcmp(coord,'Z')
    A = zeros(4,4);
    
    A11 = M;
    A12 = sum(vec(iZ-2:iZ+2));
    A14 = -0.5;
    A22 = sum(vec(iZ-2:iZ+2).^2);
    A23 = sum(vec(iZ-2:iZ+2).^3);
    A24 = -0.5*vec(iZ);
    A33 = sum(vec(iZ-2:iZ+2).^4);
    A34 = -0.5*vec(iZ)^2;
    
    A(1,:) = [A11, A12, A22, A14];
    A(2,:) = [A12, A22, A23, A24];
    A(3,:) = [A22, A23, A33, A34];
    A(4,:) = [-2*A14, -2*A24, -2*A34, 0];
    
    bz = [sum( F(iR,iZ-2:iZ+2) ); ...
        sum( F(iR,iZ-2:iZ+2)'.*vec(iZ-2:iZ+2) ); ...
        sum( F(iR,iZ-2:iZ+2)'.*(vec(iZ-2:iZ+2).^2) ); ...
        F(iR,iZ)];
    
    a = A\bz; % coefficients of quadratic expansion along the Z-axis    
end

end

function B = analyticalB(X)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

B = struct;

% Parameters of the analytical magnetic field
Bo = 5;
a = 2;% 0.6;% Minor radius in meters.
Ro = 6.0; % Major radius in meters.
qa = 10; % Safety factor at the separatrix (r=a)
co = 0.5; % Extra parameter
lamb = a/co;
Bpo = (a/Ro)*(Bo/qa)*(1-co^2)/co;
% Parameters of the analytical magnetic field

% Toroidal coordinates
% r = radius, theta = poloidal angle, phi = toroidal angle
r = sqrt( (sqrt(X(1)^2 + X(2)^2) - Ro)^2 + X(3)^2 );
theta = atan2(X(3),sqrt(X(1)^2 + X(2)^2) - Ro);
if theta < 0
    theta = theta + 2*pi;
end
phi = atan2(X(1),X(2));
if phi < 0
    phi = phi + 2*pi;
end
% Toroidal coordinates

% Poloidal magnetic field
Bp = Bpo*(r/lamb)/( 1 + (r/lamb)^2 );

eta = r/Ro;
Br = 1/( 1 + eta*cos(theta) );

Bx = Br*( Bo*cos(phi) - Bp*sin(theta)*sin(phi) );
By = -Br*( Bo*sin(phi) + Bp*sin(theta)*cos(phi) );
Bz = Br*Bp*cos(theta);


B.Bo = Bo;
B.B = [Bx,By,Bz];
end

function PP = particlePusherLeapfrog(ST)
% Relativistic particle pusher.
% Vay, J.-L. 2008, PoP 15, 056701.

PP = struct;

X = zeros(ST.params.numIt,3); % (it,ii), ii=X,y,z
v = zeros(ST.params.numIt,3);
u = zeros(ST.params.numIt,3);

% Normalization
X(1,:) = ST.params.Xo/ST.norm.l;
v(1,:) = ST.params.vo/ST.params.c;
q = ST.params.q/ST.norm.q;
m = ST.params.m/ST.norm.m;
if ST.analytical
    tmp = analyticalB(X(1,:)*ST.norm.l);
    B = tmp.B/ST.Bo;
else
    B = interpMagField(ST,X(1,:)*ST.norm.l)/ST.Bo;
end
dt = ST.params.dt*ST.norm.wc;
% E = zeros(1,3);%ST.E/(ST.Bo*ST.params.c);
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
        tmp = analyticalB(X(ii,:)*ST.norm.l);
        B = tmp.B/ST.Bo;
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
% PP.EK = sqrt(1+sum(u.^2,2));
PP.ERR = 100*(PP.EK(1) - PP.EK)./PP.EK(1);

figure
plot3(PP.X(1:ST.params.cadence:end,1),PP.X(1:ST.params.cadence:end,2),...
    PP.X(1:ST.params.cadence:end,3))
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)

time = ST.time/(2*pi/ST.params.wc);

figure
plot(time, PP.ERR)
box on
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Energy conservation [\%]','Interpreter','latex','FontSize',16)

PP.v = PP.v*ST.params.c;

end

function PP = particlePusherRungeKutta(ST)
% Relativistic particle pusher.
% 4th-order Runge-Kutta method

PP = struct;

X = zeros(ST.params.numIt,3); % (it,ii), ii=X,y,z
v = zeros(ST.params.numIt,3);
u = zeros(ST.params.numIt,3);

% Normalization
X(1,:) = ST.params.Xo/ST.norm.l;
v(1,:) = ST.params.vo/ST.params.c;
q = ST.params.q/ST.norm.q;
m = ST.params.m/ST.norm.m;
if ST.analytical
    tmp = analyticalB(X(1,:)*ST.norm.l);
    B = tmp.B/ST.Bo;
else
    B = interpMagField(ST,X(1,:)*ST.norm.l)/ST.Bo;
end
dt = ST.params.dt*ST.norm.wc;
E = zeros(1,3);%ST.E/(ST.Bo*ST.params.c);

u(1,:) = v(1,:)/sqrt(1 - sum(v(1,:).^2));


end

function PP = particlePusherMatlab(ST)
% Function to advance in time the particle's trajectory using Matlab
% solvers for ODEs.

PP = struct;

X = zeros(ST.params.numIt,3); % (it,ii), ii=X,y,z
v = zeros(ST.params.numIt,3);
u = zeros(ST.params.numIt,3);

% Normalization
X(1,:) = ST.params.Xo/ST.norm.l;
v(1,:) = ST.params.vo/ST.params.c;
q = ST.params.q/ST.norm.q;
m = ST.params.m/ST.norm.m;
if ST.analytical
    tmp = analyticalB(X(1,:)*ST.norm.l);
    B = tmp.B/ST.Bo;
else
    B = interpMagField(ST,X(1,:)*ST.norm.l)/ST.Bo;
end
E = ST.E/(ST.Bo*ST.params.c);

u(1,:) = v(1,:)/sqrt(1 - sum(v(1,:).^2));

options = odeset('RelTol',1E-5,'AbsTol',1E-10);%,'Stats','on','OutputFcn',@odeplot)
% options = odeset('RelTol',1E-3,'AbsTol',1E-6);%,'Stats','on','OutputFcn',@odeplot)

tspan = ST.time*ST.norm.wc;
y0 = [X(1,:),u(1,:)];

if ST.analytical
    [t,y] = ode45(@(t,y) LorentzForce(t,y,q,m,ST.norm.l,ST.Bo,E),tspan,y0,options);  % NONSTIFF PROBLEM
else
    [t,y] = ode45(@(t,y) LorentzForceITER(t,y,q,m,ST.norm.l,ST.Bo,E,ST),tspan,y0,options);  % NONSTIFF PROBLEM
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
plot3(PP.X(1:ST.params.cadence:end,1),PP.X(1:ST.params.cadence:end,2),...
    PP.X(1:ST.params.cadence:end,3))
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)

time = ST.time/(2*pi/ST.params.wc);

figure
plot(time, PP.ERR)
box on
xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
ylabel('Energy conservation [\%]','Interpreter','latex','FontSize',16)

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
tmp_struct = analyticalB(y(1:3)*l);
B = tmp_struct.B/Bo;

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

function DI = DiegosInvariant(ST)
% All in Gaussian units
% Parameters of the analytical magnetic field
Bo = 1E4*5;
a = 1E2*2;% Minor radius in meters.
Ro = 1E2*6.0; % Major radius in meters.
qa = 10; % Safety factor at the separatrix (r=a)
co = 0.5; % Extra parameter
lamb = a/co;
Bpo = (a/Ro)*(Bo/qa)*(1-co^2)/co;
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


figure
subplot(2,1,1)
plot(ST.time,DI)
xlabel('Time $t$ [$\tau_c$]','Interpreter','latex','FontSize',16)
ylabel('Invariant','Interpreter','latex','FontSize',16)
subplot(2,1,2)
plot(ST.time,ERR)
xlabel('Time $t$ [$\tau_c$]','Interpreter','latex','FontSize',16)
ylabel('Relative error [\%]','Interpreter','latex','FontSize',16)


R = sqrt(X(:,1).^2 + X(:,2).^2);
figure
plot(R,X(:,3))
axis equal
xlabel('R [cm]','Interpreter','latex','FontSize',16)
ylabel('Z [cm]','Interpreter','latex','FontSize',16)
title('Particle orbit','Interpreter','latex','FontSize',16)
end

