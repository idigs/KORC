% Script for studying Coulomb logarithms of fast electrons collisioning
% with bound and free electrons as well as with ions.

function CoulombLogarithms(Te,E,iSpp)
close all

% Te: background electron temperature in eV.
% E: test electron's energy in eV.
% iSpp: Array containing ion species information
%   A: atomic number
%   nn: density of neutral impurities.
%   nZ: density of impurities with mean ionization level Z.
%   Iz: ionization energy of different impurities in eV.
%   In: ionization energy of neutral iSpp
%   Z: mean ionization level of impurities.

kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep0 = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
re = qe^2/(4*pi*ep0*me*c^2);

numImpurities = numel(fieldnames(iSpp));

% Unit conversion
TeJ = Te*qe; % Converted to Joules

EJ = E*qe; % Converted to Joules
g = EJ/(me*c^2);

lD = @(ne,Te) sqrt(ep0*(Te*qe)/(qe^2*ne));

CLog = @(ne,Te) 32.2 - 0.5*log(ne) + log(Te); % KORC's model

CLogee = @(fb) 20*(1 - fb/2);

CLogee_f = @(g,ne,Te) log( (g-1)*sqrt(g+1)*lD(ne,Te)/(2*g*re) ); % Mosher's model

CLogee_b = @(g,Iz) log( (g-1)*sqrt(g+1)*me*c^2./(Iz*qe) );

CLogeZ = @(g,Iz,ne,Te,Z) log( lD(ne,Te)*(Iz*qe)./(Z*re*me*c^2) );

CLogeZ0 = @(g,Iz) log( (g^2 - 1)*me*c^2./(g*(Iz*qe)) );

% Models for E critical

ECH = @(ne,Te) qe^3*ne*CLog(ne,Te)/(4*pi*ep0^2*me*c^2);

Ec4c = @(ne,fb) qe^3*ne*CLogee(fb)/(4*pi*ep0^2*me*c^2);

Ec4d = @(g,nef,Te,neb,Iz) qe^3*(nef*CLogee_f(g,nef,Te) + sum(neb.*CLogee_b(g,Iz)))/(4*pi*ep0^2*me*c^2);



% Models for Zeff

% Zeff_simple = (iZ.D.nz + sum(iSpp.nz.*iSpp.Z.^2))/ne;

nef = 0;
neb = 0;
neb_a = [];
Iz_a = [];
for ii=1:numImpurities
    nef = nef + sum(iSpp.(['spp' num2str(ii)]).nz.*iSpp.(['spp' num2str(ii)]).Z);
    neb = neb + iSpp.(['spp' num2str(ii)]).nn*iSpp.(['spp' num2str(ii)]).A + ...
        sum(iSpp.(['spp' num2str(ii)]).nz.*(iSpp.(['spp' num2str(ii)]).A-iSpp.(['spp' num2str(ii)]).Z));
    
    if ~isempty(iSpp.(['spp' num2str(ii)]).Iz)
        neb_a = [neb_a,iSpp.(['spp' num2str(ii)]).nn, iSpp.(['spp' num2str(ii)]).nz];
        Iz_a = [Iz_a,iSpp.(['spp' num2str(ii)]).In, iSpp.(['spp' num2str(ii)]).Iz];
    else
        neb_a = [neb_a,iSpp.(['spp' num2str(ii)]).nn];
        Iz_a = [Iz_a,iSpp.(['spp' num2str(ii)]).In];
    end
end
fb = neb/(nef + neb);

Ec1 = zeros(1,numel(g));
Ec2 = zeros(1,numel(g));
Ec3 = zeros(1,numel(g));
for gg=1:numel(g)
    Ec1(gg) = ECH(nef,Te);
    Ec2(gg) = Ec4c(nef+neb,fb);
    Ec3(gg) = Ec4d(g(gg),nef,Te,neb_a,Iz_a);
end


EAxis = E/1E6;

figure
subplot(2,1,1)
plot(EAxis,Ec1,'k',EAxis,Ec2,'r',EAxis,Ec3,'b')
legend({'$E_{CH}$','$E_{crit} (Parks)$','$E_{crit} (Mosher)$'},'Interpreter','latex')
xlabel('Energy (MeV)','Interpreter','latex')
ylabel('$E_{crit}$ (V/m)','Interpreter','latex')


end


% iSpp.spp1.A = 18; % Argon
% iSpp.spp1.nn = 0.002611*1E13*1E6; % neutral density in m^-3
% iSpp.spp1.In = 15.7596117; % Ionization energy of neutral Argon.
% iSpp.spp1.nz = [0.284,1.17]*1E13*1E6; % density of partially ionized impurities.
% iSpp.spp1.Z = [1,2]; % Average charge state.
% iSpp.spp1.Iz = [27.62967,40.735]; % Ionization energy.
% 
% iSpp.spp2.A = 1; % Deuterium
% iSpp.spp2.nn = 0.04481*1E13*1E6; % neutral density in m^-3
% iSpp.spp2.In = 15.46658; % Ionization energy of neutral Argon.
% iSpp.spp2.nz = 3.764*1E13*1E6; % Ionized deuterium.
% iSpp.spp2.Z = 1;
% iSpp.spp2.Iz = [];
