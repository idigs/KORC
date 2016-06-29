function cp = collisions_parameters(Te,ne,nz)
close all
% Collision parameters for RE in fusion plasmas
kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
e0 = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass

Te = Te*qe;
rD = sqrt( e0*Te/(ne*qe^2) );
re = qe^2/(4*pi*e0*me*c^2);

Eo = linspace(1E6,100E6,1E3);
gamma = qe*Eo/(me*c^2);
Eo = Eo/1E6;

% Iz = qe*[15.7596,27.62965,40.74,59.81,75.02]; % Argon
% Iz = qe*[21.5646,40.96296,63.45,97.12,126.21]; % Neon
% Zj = 1:1:numel(Iz);

% % % % Argon
% Iz = qe*(27.62965);
% Zj = 2.0;
% Zo = 18.0;

% % % Neon
Iz = qe*(21.5646);
Zj = 1.0;
Zo = 10.0;

nH = ne;
nef = ne + Zj*nz;
neb = ne + (Zo - Zj)*nz;

ae = zeros(size(Eo));
ai = zeros(size(Eo));

aux = (gamma - 1).*sqrt(gamma + 1);
lratio = rD/re;

Clog_ef = log( 0.5*lratio*aux./gamma );
Clog_eb = log( aux*(me*c^2/Iz) );

ae = nef*Clog_ef + neb*Clog_eb;

figure
semilogx(Eo,Clog_ef,Eo,Clog_eb,'--')
box on
grid on

Clog_eH = log( lratio*(gamma.^2 - 1)./gamma );
Clog_eZ = log( (lratio/Zj)*(Iz/(me*c^2)) );
Clog_eZo = log( (me*c^2/Iz)*(gamma.^2 - 1)./gamma );

ai = nH*Clog_eH + nz*(Clog_eZ*Zj^2 + Clog_eZo*Zo^2);

figure
semilogx(Eo,Clog_eH,Eo,Clog_eZ*ones(size(Eo)),'--',Eo,Clog_eZo,':')
box on
grid on

ae = ae/1E20;
ai = ai/1E20;

Zeff = (nH + nz*Zj^2)/(nH + nz*Zj);

figure
subplot(2,1,1)
semilogx(Eo,ae,'r',Eo,ai,'b')
box on
grid on
subplot(2,1,2)
semilogx(Eo,ai./ae,Eo,Zeff*ones(size(Eo)))
box on
grid on


end