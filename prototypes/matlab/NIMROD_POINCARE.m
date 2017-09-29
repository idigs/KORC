function NIMROD_POINCARE(timeStepParams,numInitCond)

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
DS = timeStepParams(2);

ST.analytical = false;

load('NIMROD_LIMITED_1150.mat','S')

B = struct;
B.ND = '3D';
B.NR = numel(S.R);
B.NZ = numel(S.Z);
B.Nphi = numel(S.PHI);
B.R = S.R;
B.Z = S.Z;
B.phi = S.PHI;
B.BR = S.BR;
B.BZ = S.BZ;
B.Bphi = S.BPHI;

B.Ro = S.Ro;
B.Zo = S.Zo;
B.Bo = S.Bo;

R = zeros(size(S.BR));
PHI = zeros(size(S.BR));
Z = zeros(size(S.BR));

for rr=1:B.NR
    for pp=1:B.Nphi
        for zz=1:B.NZ
            R(rr,pp,zz) = S.R(rr);
            PHI(rr,pp,zz) = S.PHI(pp);
            Z(rr,pp,zz) = S.Z(zz);
        end
    end
end


B.SI.BR = scatteredInterpolant(reshape(R,[numel(R) 1]),reshape(PHI,[numel(PHI) 1]),reshape(Z,[numel(Z) 1]),reshape(S.BR,[numel(S.BR) 1]));
B.SI.Bphi = scatteredInterpolant(reshape(R,[numel(R) 1]),reshape(PHI,[numel(PHI) 1]),reshape(Z,[numel(Z) 1]),reshape(S.BPHI,[numel(S.BPHI) 1]));
B.SI.BZ = scatteredInterpolant(reshape(R,[numel(R) 1]),reshape(PHI,[numel(PHI) 1]),reshape(Z,[numel(Z) 1]),reshape(S.BZ,[numel(S.BZ) 1]));

% load('B.mat')

P = figure;

Ros = linspace(1.16,2.2,numInitCond);%1.6*ones(1,numInitCond);
phio = zeros(1,numInitCond);
Zos = -0.08*ones(1,numInitCond);%linspace(-0.8,0.8,numInitCond);

for ii=1:numInitCond
    
    disp(num2str(ii))
    
    R = zeros(1,numIt);
    phi = zeros(1,numIt);
    Z = zeros(1,numIt);
    
    R(1) = Ros(ii);
    phi(1) = phio(ii);
    Z(1) = Zos(ii);
    
    options = odeset('RelTol',1E-6,'AbsTol',1E-10);
    
    tspan = (0:1:numIt)*DS;
    yo = [R(1), phi(1), Z(1)];
    
    [~,Y] = ode45(@(t,y) eqnsOfMotion2(t,y,B),tspan,yo,options);  % RK4
    
    R = Y(:,1);
    phi = mod(Y(:,2),2*pi);
    Z = Y(:,3);
    
    locs = find(abs(diff(phi)) > 6);
    figure(P)
    hold on
    plot(R(locs),Z(locs),'.','MarkerSize',3)
    hold off
    
end

figure(P)
axis equal
xlabel('$R$ [m]','Interpreter','latex','FontSize',16)
ylabel('$Z$ [m]','Interpreter','latex','FontSize',16)
title('Poincare plot','Interpreter','latex','FontSize',16)
saveas(P,'poincare_limited_1150_2','fig')

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

BR = B_ST.SI.BR(X(1),mod(X(2),2*pi),X(3));
Bphi = B_ST.SI.Bphi(X(1),mod(X(2),2*pi),X(3));
BZ = B_ST.SI.BZ(X(1),mod(X(2),2*pi),X(3));

BF = [BR,Bphi,BZ];
end
