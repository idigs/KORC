function integrateMagneticField_RK(pathToBField,fileType,ND,res,timeStepParams,numInitCond)
% name='fields/bfield_tracing/d3d_bfield_tracing-3E-3.dat';
% integrateMagneticField_RK(name,'SIESTA','3D',[100,79,149],[1E4,1E-2],10)
% ST = pOrbs('','','2D',[],[4E5,1E-2,10],[-1,1],[2.1,0,0],[0.99,0],true);

narginchk(6,6);

close all

% Plasma parameters and physical constants, all in SI units
ST.params = struct;
ST.params.kB = 1.38E-23; % Boltzmann constant
ST.params.mu0 = (4E-7)*pi; % Magnetic permeability
ST.params.ep = 8.854E-12;% Electric permitivity
ST.params.c = 2.9979E8; % Speed of light
ST.params.qe = 1.602176E-19; % Electron charge
ST.params.me = 9.109382E-31; % Electron mass

numIt = timeStepParams(1);
DT = timeStepParams(2);

if strcmp(pathToBField,'')
    ST.analytical = true;
else
    ST.analytical = false; % true = analytical B field; false = xpander fields
    ST.pathToBField = pathToBField; % Path to xpander fields
    ST.fileType = fileType;
end

if ST.analytical
    B = analyticalB([1,1,1],'initialize');
else
    B = loadMagneticField(pathToBField,fileType,ND,res);
end

P = figure;

Ros = linspace(1.61,2.1,numInitCond);
phio = 0;
Zos = 0; % linspace(-1,1,numInitCond);

for ii=1:numel(Ros)
    
    disp(num2str(ii))
    
    R = zeros(1,numIt);
    phi = zeros(1,numIt);
    Z = zeros(1,numIt);
    
    R(1) = Ros(ii);
    phi(1) = phio;
    Z(1) = Zos; % Zos(ii);
    
    options = odeset('RelTol',1E-6,'AbsTol',1E-10);
    
    tspan = (0:1:numIt)*DT;
    y0 = [R(1), phi(1), Z(1)];    
    
    if ST.analytical
        [s,y] = ode45(@(s,y) eqnsOfMotion1(s,y),tspan,y0,options);  % RK4
    else
        [~,y] = ode45(@(t,y) eqnsOfMotion2(t,y,B),tspan,y0,options);  % RK4
    end
    
    R = y(:,1);
    phi = mod(y(:,2),2*pi);
    Z = y(:,3);
    
    locs = find(abs(diff(phi)) > 6);
    figure(P)
    hold on
%     plot(R,Z,'r.','MarkerSize',6)
    plot(R(locs),Z(locs),'b.','MarkerSize',2)
    hold off
    
    if ST.analytical
        disp('Toroidal coordinate system')
        % Cylindrical to Cartesian
        x = R.*cos(phi);
        y = R.*sin(phi);
        z = Z;
        
        h = figure;
        subplot(3,1,1)
        plot3(x,y,z,'b')
        axis equal
        box on; grid on
        xlabel('$X$ [m]','Interpreter','latex','FontSize',16)
        ylabel('$Y$ [m]','Interpreter','latex','FontSize',16)
        zlabel('$Z$ [m]','Interpreter','latex','FontSize',16)
        
        tmp = sqrt(x.^2 + y.^2) - B.Ro;
        r = sqrt( tmp.^2 + z.^2 );
        eta = r/B.Ro;
        zeta = pi/2 - phi;
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        theta = atan2(Z,R - B.Ro);
        theta(theta<0) = theta(theta<0) + 2*pi;
        
        Bp = B.Bpo*(r/B.lamb)./( 1 + (r/B.lamb).^2 );
        Dtheta = [0;mod(diff(theta),2*pi)];
        zeta_theory = zeta(1) + cumsum( (B.Bo./Bp).*(eta./(1 + eta.*cos(theta))).*Dtheta );
        zeta_theory = mod(zeta_theory,2*pi);
        
        figure(h)
        subplot(3,1,2)
        plot(theta,zeta,'r.',theta,zeta_theory,'b.')
        axis([0 2*pi 0 2*pi])
        xlabel('$\theta$ [rad]','Interpreter','latex','FontSize',16)
        ylabel('$\zeta$ [rad]','Interpreter','latex','FontSize',16)
        
        
        % Theoretical curvature of the analytical magnetic field
        k = abs( B.Bo*Bp.*eta.*sin(theta)./(sqrt( B.Bo^2 + Bp.^2 ).*(1 + eta.*cos(theta)).^2 ));
        
        figure(h)
        subplot(3,1,3)
        plot(theta,k,'k.')
        xlim([0 2*pi])
        xlabel('$\theta$ [rad]','Interpreter','latex','FontSize',16)
        ylabel('$\kappa(\theta)$','Interpreter','latex','FontSize',16)
        % Theoretical curvature of the analytical magnetic field
        
        % Numerical curvature of the analytical magnetic field
        
        
        % Numerical curvature of the analytical magnetic field
        
    end
    
end

figure(P)
axis equal
xlabel('$R$ [m]','Interpreter','latex','FontSize',16)
ylabel('$Z$ [m]','Interpreter','latex','FontSize',16)
title('Poincare plot','Interpreter','latex','FontSize',16)

end

function B = analyticalB(Y,opt)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

narginchk(1,2);

% Parameters of the analytical magnetic field
Bo = 3;
a = 0.5;% Minor radius in meters.
Ro = 1.6; % Major radius in meters.
qa = 5; % Safety factor at the separatrix (r=a)
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
%     y0 = [R(1), phi(1), Z(1)]; 
    X = [Y(1)*cos(Y(2)),Y(1)*sin(Y(2)),Y(3)];    
    
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
    
    % Now we change to cylindrical coordinates
    % Cylindrical coordinates
    % R = radius, phi = azimuthal angle, Z = Z coordinate
    R = sqrt(X(1)^2 + X(2)^2);
    Z = X(3);
    phi = atan2(X(2),X(1));
    if phi < 0
        phi = phi + 2*pi;
    end
    
    % B = [B^R,B^phi,B^Z]
    B = [Bx*cos(phi) + By*sin(phi),...
        By*cos(phi) - Bx*sin(phi),...
        Bz];
end

end

function B = loadMagneticField(pathToBField,fileType,ND,res);
% All quantities in SI units
B = struct;
B.fileType = fileType;
B.ND = ND;

data = load(pathToBField);

switch fileType
    case 'RAW' % This option changes quite frequently
        B.R = data(:,1);
        B.Z = data(:,3);
        
        B.BR = data(:,4);
        B.Bphi = - data(:,5); % minus sign
        B.BZ = data(:,6);        
    case 'SIESTA'
        if strcmp(ND,'2D')
            NR = res(1);
            Nphi = res(2);
            NZ = res(3);
            
            R = data(:,1);
            Z = data(:,2);
            
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
        elseif strcmp(ND,'3D')
            B.NR = res(1);
            B.Nphi = res(2);
            B.NZ = res(3);
            
            R = data(:,1);
            Z = data(:,2);
            phi = ((0:1:B.Nphi-1) + 0.5)*(2*pi/B.Nphi);            
            R = reshape(R,B.NR,B.Nphi,B.NZ);
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
                B.phi(:,:,iphi) = phi(iphi);
                B.Z(:,:,iphi) = squeeze(Z(:,iphi,:));
                
                B.BR(:,:,iphi) = squeeze(BR(:,iphi,:));
                B.Bphi(:,:,iphi) = squeeze(Bphi(:,iphi,:));
                B.BZ(:,:,iphi) = squeeze(BZ(:,iphi,:));
            end
        else
            error('Please use 2D or 3D fields');
        end
        
    case 'VMEC'
        
        if strcmp(ND,'2D')
            B.NR = res(1);
            B.NZ = res(2);
            
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
            
        elseif strcmp(ND,'3D')
            error('Not ready for using 3D fields of VMEC!')
        else
            error('Please use 2D or 3D fields');
        end
        
    case 'XPANDER'
        
        if strcmp(ND,'2D')
            B.NR = res(1);
            B.NZ = res(2);
            
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
        elseif strcmp(ND,'3D')
            B.NR = res(1);
            B.Nphi = res(2);
            B.NZ = res(3);
            
            B.R = zeros(B.NR,B.NZ,B.Nphi);
            B.Z = zeros(B.NR,B.NZ,B.Nphi);
            B.phi = zeros(B.NR,B.NZ,B.Nphi);
            
            phi = ((0:1:B.Nphi-1) + 0.5)*(2*pi/B.Nphi);
            
            B.B = zeros(B.NR,B.NZ,B.Nphi); % magnitude
            B.BR = zeros(B.NR,B.NZ,B.Nphi);
            B.Bphi = zeros(B.NR,B.NZ,B.Nphi);
            B.BZ = zeros(B.NR,B.NZ,B.Nphi);
            
            for iphi = 1:B.Nphi;
                for iz=1:B.NZ
                    indi = (iphi-1)*B.NR*B.NZ + (iz-1)*B.NR + 1;
                    indf = (iphi-1)*B.NR*B.NZ + iz*B.NR;
                    B.R(:,iz,iphi) = data(indi:indf,1);
                    B.Z(:,iz,iphi) = data(indi:indf,3);
                    
                    B.BR(:,iz,iphi) = data(indi:indf,4);
                    B.Bphi(:,iz,iphi) = data(indi:indf,5);
                    B.BZ(:,iz,iphi) = data(indi:indf,6);
                end
                B.phi(:,:,iphi) = phi(iphi);
            end
        else
            error('Please use 2D or 3D fields');
        end
        
    otherwise
        error('Unknown file format!')
end

B.SI = calculateScatteredInterpolant(B);

B.R = [];
B.Z = [];
B.phi = [];

B.B = [];
B.BR = [];
B.Bphi = [];
B.BZ = [];
B.P = [];

end

function SI = calculateScatteredInterpolant(B)
% Cylindrical coordinates
% R = radius, phi = azimuthal angle, Z = z coordinate
% calculate interpolant of the field B
disp('*** Uniform grid for the magnetic field detected! ***')
disp('Switching to scattered interpolant...')
SI = struct;


if strcmp(B.ND,'2D')
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
    
elseif strcmp(B.ND,'3D')

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

disp('Scattered interpolant: done!')
end

function dyds = eqnsOfMotion1(s,y)
% The equation of motion is: dy/ds = b, where 'y' is the position vector
% 's' is the length of the orbit, and b is the unitary magnetic field
% vector (B = |B| b).
% y(1) = R;
% y(2) = phi;
% y(3) = Z;

B = analyticalB(y);
B = B/sqrt(sum(B.^2));

dyds = zeros(3,1);

dyds(1) = B(1);
dyds(2) = B(2)/y(1);
dyds(3) = B(3);
end

function dydt = eqnsOfMotion2(t,y,B_ST)
% dydt = B;
% y(1) = R;
% y(2) = phi;
% y(3) = Z;

B = interpMagField(B_ST,y);

dydt = zeros(3,1);

dydt(1) = B(1);
dydt(2) = B(2)/y(1);
dydt(3) = B(3);
end

function BF = interpMagField(B_ST,X)
% Cylindrical coordinates
% X(1) = R, X(2) = phi, X(3) = Z

BR = B_ST.SI.BR(X(1),X(3),mod(X(2),2*pi));
Bphi = B_ST.SI.Bphi(X(1),X(3),mod(X(2),2*pi));
BZ = B_ST.SI.BZ(X(1),X(3),mod(X(2),2*pi));

BF = [BR,Bphi,BZ];
end
