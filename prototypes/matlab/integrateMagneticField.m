function integrateMagneticField(timeStepParams,xo)
% integrateMagneticField([1E3,1E-1,8],[1.5,0,0])

narginchk(2,2);

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
subSteps = timeStepParams(3);

B = loadMagneticField('jfit_165365_1400.mat');
B.SI = calculatescatteredInterpolant(B);

P = figure;

Ros = linspace(B.Ro,2.25,25);
Zos = linspace(-0.3,0.3,30);

for ii=1:numel(Zos)
    
    disp(num2str(ii))
    
    R = zeros(1,numIt);
    phi = zeros(1,numIt);
    Z = zeros(1,numIt);
    
    R(1) = 1.8;%Ros(ii);
    phi(1) = 0;
    Z(1) = Zos(ii);
    
    for it=2:numIt
        
        %     tmp_R = zeros(1,subSteps);
        %     tmp_phi = zeros(1,subSteps);
        %     tmp_Z = zeros(1,subSteps);
        %     tmp_R(1) = R(it-1);
        %     tmp_phi(1) = phi(it-1);
        %     tmp_Z(1) = Z(it-1);
        
        tmp_R = R(it-1);
        tmp_phi = phi(it-1);
        tmp_Z = Z(it-1);
        
        for nn=2:subSteps
            %         tmp_R(nn) = R(it-1) + ...
            %             DT*B.SI.BR(0.5*(R(it-1) + tmp_R(nn-1)),0.5*(Z(it-1) + tmp_Z(nn-1)));
            %         tmp_phi(nn) = phi(it-1) + ...
            %             DT*B.SI.Bphi(0.5*(R(it-1) + tmp_R(nn-1)),0.5*(Z(it-1) + tmp_Z(nn-1)));
            %         tmp_Z(nn) = Z(it-1) + ...
            %             DT*B.SI.BZ(0.5*(R(it-1) + tmp_R(nn-1)),0.5*(Z(it-1) + tmp_Z(nn-1)));
            
            tmp_R = R(it-1) + ...
                DT*B.SI.BR(0.5*(R(it-1) + tmp_R),0.5*(Z(it-1) + tmp_Z));
            tmp_phi = phi(it-1) + DT*B.Ro*B.Bo/(0.5*(R(it-1) + tmp_R));
            tmp_Z = Z(it-1) + ...
                DT*B.SI.BZ(0.5*(R(it-1) + tmp_R),0.5*(Z(it-1) + tmp_Z));
        end
        
        R(it) = tmp_R(end);
        phi(it) = mod(tmp_phi(end),2*pi);
        Z(it) = tmp_Z(end);
    end
    
    locs = find(abs(diff(phi)) > 6);
    figure(P)
    hold on
    plot(R(locs),Z(locs),'k.','MarkerSize',7)
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
