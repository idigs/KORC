clear all
% close all
clc


kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass

Bo = 1.131881701794893; % In teslas
Eo = 0.05606412092419791; % In V/m
ne = 2.898E19; % In particles/m3
Clog = 16.640712846994347;

Gamma = ne*qe^4*Clog/(4*pi*ep^2);

NE = 1E3;
NP = 1E3;

E_max = 20; % In MeV
E_min = 1;

pitch_max = 80; % In degrees
pitch_min = 0;  % In degrees


E = 1E6*linspace(E_min,E_max,NE)*qe;
g = E/(me*c^2);
v = c*sqrt(1-1./g.^2);
pitch = deg2rad(linspace(pitch_min,pitch_max,NP));

PTot = zeros(NE,NP);

for ee=1:NE
    for pp=1:NP
        PR = -(qe^4*Bo^2/(6*pi*ep*me^2*c^3))*(v(ee)*g(ee)).^2.*sin(pitch(pp)).^2;
        PE = qe*v(ee)*Eo;
        Pcoll = -(1 + g(ee))*Gamma/(g(ee)*v(ee)*me);
        
        PTot(ee,pp) = PR + PE + Pcoll;
    end
end

P_min = min(min(PTot));
P_max = max(max(PTot));

levels = linspace(P_min,P_max,10);
[~,I] = min(abs(levels));
if (levels(I) > 0)
    levels = [levels(1:I-1) 0 levels(I:end)];
else
    levels = [levels(1:I) 0 levels(I+1:end)];
end

yAxis = E/(1E6*qe);
xAxis = rad2deg(pitch);

figure
contourf(xAxis,yAxis,PTot,levels,'ShowText','On')
xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex')
colormap(jet);cb = colorbar;
ylabel(cb,'Total power (W)','Interpreter','latex')
