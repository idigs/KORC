function integrateMagneticField(timeStepParams,numInitCond,opt)
% integrateMagneticField([1E3,1E-2,8],30)

narginchk(2,3);

% close all

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
subSteps = timeStepParams(3);

if (nargin == 2)
    B = loadMagneticField('jfit_165365_1400.mat');
elseif  (nargin == 3) && strcmp(opt,'refine')
    B = loadAndRefineMagneticField('jfit_165365_1400.mat');
end
B.SI = calculatescatteredInterpolant(B);

P = figure;

Ros = [1.358];%linspace(1,1.823,numInitCond);
phio = pi;
Zos = 0;%linspace(-1,1,numInitCond);

for ii=1:numel(Ros)
    
    disp(num2str(ii))
    
    R = zeros(1,numIt);
    phi = zeros(1,numIt);
    Z = zeros(1,numIt);
    
    R(1) = Ros(ii);
    phi(1) = phio;
    Z(1) = Zos;%Zos(ii);
    
    for it=2:numIt
        
        tmp_R = zeros(1,subSteps);
        tmp_phi = zeros(1,subSteps);
        tmp_Z = zeros(1,subSteps);
        tmp_R(1) = R(it-1);
        tmp_phi(1) = phi(it-1);
        tmp_Z(1) = Z(it-1);
        
%         tmp_R = R(it-1);
%         tmp_phi = phi(it-1);
%         tmp_Z = Z(it-1);
        
        for nn=2:subSteps
            tmp_R(nn) = R(it-1) + ...
                DT*B.SI.BR(0.5*(R(it-1) + tmp_R(nn-1)),0.5*(Z(it-1) + tmp_Z(nn-1)));
            tmp_phi(nn) = phi(it-1) + DT*B.Ro*B.Bo/(0.5*(R(it-1) + tmp_R(nn-1)))^2;
            tmp_Z(nn) = Z(it-1) + ...
                DT*B.SI.BZ(0.5*(R(it-1) + tmp_R(nn-1)),0.5*(Z(it-1) + tmp_Z(nn-1)));
            
%             tmp_R = R(it-1) + ...
%                 DT*B.SI.BR(0.5*(R(it-1) + tmp_R),0.5*(Z(it-1) + tmp_Z));
%             tmp_phi = phi(it-1) + DT*B.Ro*B.Bo/(0.5*(R(it-1) + tmp_R));
%             tmp_Z = Z(it-1) + ...
%                 DT*B.SI.BZ(0.5*(R(it-1) + tmp_R),0.5*(Z(it-1) + tmp_Z));
        end
        
        R(it) = tmp_R(end);
        phi(it) = mod(tmp_phi(end),2*pi);
        Z(it) = tmp_Z(end);
    end
    
    locs = find(abs(diff(phi)) > 6);
    figure(P)
    hold on
    plot(R,Z,'k.','MarkerSize',6)
    plot(R(locs),Z(locs),'g.','MarkerSize',6)
    hold off
end

figure(P)
axis equal
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
title('Poincare plot','Interpreter','latex','FontSize',16)

end

function B = loadMagneticField(name)
% Here we calculate the magnetic field components using centered finite
% differences.
B = struct;

DATA = load(name);

B.NZ = numel(DATA.zg);
B.NR = numel(DATA.rg);

DZ = mean(diff(DATA.zg));
DR = mean(diff(DATA.rg));

B.BR = zeros(B.NZ,B.NR);
B.BZ = zeros(B.NZ,B.NR);

[B.R,B.Z] = meshgrid(DATA.rg,DATA.zg);

B.Bo = 2.19;
B.Ro = 1.695;
B.Bphi = - B.Bo*B.Ro./B.R;

B.BR(3:B.NZ-2,:) = ( - DATA.psig(5:B.NZ,:)/12 + 2*DATA.psig(4:B.NZ-1,:)/3 ...
    - 2*DATA.psig(2:B.NZ-3,:)/3 + DATA.psig(1:B.NZ-4,:)/12  )/DZ;
B.BR(1:2,:) = ( - 25*DATA.psig(1:2,:)/12 + 4*DATA.psig(2:3,:) - 3*DATA.psig(3:4,:) ...
    + 4*DATA.psig(4:5,:)/3 - DATA.psig(5:6,:)/4 )/DZ;
B.BR(B.NZ-1:B.NZ,:) = ( 25*DATA.psig(B.NZ-1:B.NZ,:)/12 - 4*DATA.psig(B.NZ-2:B.NZ-1,:) + 3*DATA.psig(B.NZ-3:B.NZ-2,:) ...
    - 4*DATA.psig(B.NZ-4:B.NZ-3,:)/3 + DATA.psig(B.NZ-5:B.NZ-4,:)/4 )/DZ;
B.BR = B.BR./B.R;

B.BZ(:,3:B.NR-2) = ( - DATA.psig(:,5:B.NR)/12 + 2*DATA.psig(:,4:B.NR-1)/3 ...
    - 2*DATA.psig(:,2:B.NR-3)/3 + DATA.psig(:,1:B.NR-4)/12  )/DR;
B.BZ(:,1:2) = ( - 25*DATA.psig(:,1:2)/12 + 4*DATA.psig(:,2:3) - 3*DATA.psig(:,3:4) ...
    + 4*DATA.psig(:,4:5)/3 - DATA.psig(:,5:6)/4 )/DR;
B.BZ(:,B.NR-1:B.NR) = ( 25*DATA.psig(:,B.NR-1:B.NR)/12 - 4*DATA.psig(:,B.NR-2:B.NR-1) + 3*DATA.psig(:,B.NR-3:B.NR-2) ...
    - 4*DATA.psig(:,B.NR-4:B.NR-3)/3 + DATA.psig(:,B.NR-5:B.NR-4)/4 )/DR;
B.BZ = - B.BZ./B.R;


% Here we plot the magnetic field components using centered finite
% differences.
figure
subplot(2,2,1)
surfc(B.R,B.Z,B.BR,'LineStyle','none')
view([0,90])
axis equal
box on
colorbar
title('$B^R$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
subplot(2,2,2)
surfc(B.R,B.Z,B.BZ,'LineStyle','none')
view([0,90])
axis equal
box on
colorbar
title('$B^Z$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
subplot(2,2,3)
surfc(B.R,B.Z,B.Bphi,'LineStyle','none')
view([0,90])
axis equal
box on
colorbar
title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
subplot(2,2,4)
contour(B.R,B.Z,DATA.psig,20)
view([0,90])
axis equal
box on
colorbar
title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
colormap(jet)
end

function B = loadAndRefineMagneticField(name)
% Here we calculate the magnetic field components using centered finite
% differences.
B = struct;

DATA = load(name);

[rg,zg] = meshgrid(DATA.rg,DATA.zg);
rg = reshape(rg,[numel(rg) 1]);
zg = reshape(zg,[numel(zg) 1]);
psig = reshape(DATA.psig,[numel(DATA.psig) 1]);
F = scatteredInterpolant(rg,zg,psig);


Rmin = min(DATA.rg);
Rmax = max(DATA.rg);
R = linspace(Rmin,Rmax,5*numel(DATA.rg));

Zmin = min(DATA.zg);
Zmax = max(DATA.zg);
Z = linspace(Zmin,Zmax,5*numel(DATA.zg));

B.NZ = numel(Z);
B.NR = numel(R);

DZ = mean(diff(Z));
DR = mean(diff(R));

[B.R,B.Z] = meshgrid(R,Z);
PSIp = F(B.R,B.Z);

B.BR = zeros(B.NZ,B.NR);
B.BZ = zeros(B.NZ,B.NR);

B.Bo = 2.19;
B.Ro = 1.695;
B.Bphi = - B.Bo*B.Ro./B.R;

B.BR(3:B.NZ-2,:) = ( - PSIp(5:B.NZ,:)/12 + 2*PSIp(4:B.NZ-1,:)/3 ...
    - 2*PSIp(2:B.NZ-3,:)/3 + PSIp(1:B.NZ-4,:)/12  )/DZ;
B.BR(1:2,:) = ( - 25*PSIp(1:2,:)/12 + 4*PSIp(2:3,:) - 3*PSIp(3:4,:) ...
    + 4*PSIp(4:5,:)/3 - PSIp(5:6,:)/4 )/DZ;
B.BR(B.NZ-1:B.NZ,:) = ( 25*PSIp(B.NZ-1:B.NZ,:)/12 - 4*PSIp(B.NZ-2:B.NZ-1,:) + 3*PSIp(B.NZ-3:B.NZ-2,:) ...
    - 4*PSIp(B.NZ-4:B.NZ-3,:)/3 + PSIp(B.NZ-5:B.NZ-4,:)/4 )/DZ;
B.BR = B.BR./B.R;

B.BZ(:,3:B.NR-2) = ( - PSIp(:,5:B.NR)/12 + 2*PSIp(:,4:B.NR-1)/3 ...
    - 2*PSIp(:,2:B.NR-3)/3 + PSIp(:,1:B.NR-4)/12  )/DR;
B.BZ(:,1:2) = ( - 25*PSIp(:,1:2)/12 + 4*PSIp(:,2:3) - 3*PSIp(:,3:4) ...
    + 4*PSIp(:,4:5)/3 - PSIp(:,5:6)/4 )/DR;
B.BZ(:,B.NR-1:B.NR) = ( 25*PSIp(:,B.NR-1:B.NR)/12 - 4*PSIp(:,B.NR-2:B.NR-1) + 3*PSIp(:,B.NR-3:B.NR-2) ...
    - 4*PSIp(:,B.NR-4:B.NR-3)/3 + PSIp(:,B.NR-5:B.NR-4)/4 )/DR;
B.BZ = - B.BZ./B.R;


% Here we plot the magnetic field components using centered finite
% differences.
figure
subplot(2,2,1)
surfc(B.R,B.Z,B.BR,'LineStyle','none')
view([0,90])
axis equal
box on
colorbar
title('$B^R$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
subplot(2,2,2)
surfc(B.R,B.Z,B.BZ,'LineStyle','none')
view([0,90])
axis equal
box on
colorbar
title('$B^Z$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
subplot(2,2,3)
surfc(B.R,B.Z,B.Bphi,'LineStyle','none')
view([0,90])
axis equal
box on
colorbar
title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
subplot(2,2,4)
contour(B.R,B.Z,PSIp,20)
view([0,90])
axis equal
box on
colorbar
title('$B^\phi$ [T]','Interpreter','latex','FontSize',16)
xlabel('R [m]','Interpreter','latex','FontSize',16)
ylabel('Z [m]','Interpreter','latex','FontSize',16)
colormap(jet)

data = zeros(1,6);
fileID = fopen([name '.txt'],'w');
for ir=1:B.NR
    for iz=1:B.NZ
        data(1) = B.R(iz,ir);
        data(2) = 0;
        data(3) = B.Z(iz,ir);
        
        data(4) = B.BR(iz,ir);
        data(5) = B.Bphi(iz,ir);
        data(6) = B.BZ(iz,ir);
        
        fprintf(fileID,'%3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\t %3.7f\n',data);
    end
end
fclose(fileID);
disp(['File: ``' name '.txt`` generated!'])

end

function SI = calculatescatteredInterpolant(B)
% Cylindrical coordinates
% R = radius, phi = azimuthal angle, Z = z coordinate
disp('Scattered interpolant...')
SI = struct;

R = reshape(B.R,[numel(B.R) 1]);
Z = reshape(B.Z,[numel(B.Z) 1]);

DATA = reshape(B.BR,[numel(B.BR) 1]);
SI.BR = scatteredInterpolant(R,Z,DATA);
clear DATA

DATA = reshape(B.Bphi,[numel(B.Bphi) 1]);
SI.Bphi = scatteredInterpolant(R,Z,DATA);
clear DATA

DATA = reshape(B.BZ,[numel(B.BZ) 1]);
SI.BZ = scatteredInterpolant(R,Z,DATA);
clear DATA

disp('Scattered interpolant: done!')
end
