% Script for studying Coulomb logarithms of fast electrons collisioning
% with bound and free electrons as well as with ions.

function CoulombLogarithms(Te,E,iSpp)
% CoulombLogarithms(2.0,linspace(1E6,40E6,100),iSpp)
close all

% Te: background electron temperature in eV.
% E: test electron's kinetic energy in eV.
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

if isfield(iSpp,'H')
    numImpurities = numel(fieldnames(iSpp)) - 1;
else
    numImpurities = numel(fieldnames(iSpp));
end

% Unit conversion
TeJ = Te*qe; % Converted to Joules

EJ = E*qe + me*c^2; % Converted to Joules
g = EJ/(me*c^2);

lD = @(ne,Te) sqrt(ep0*(Te*qe)/(qe^2*ne));

CLog = @(ne,Te) 32.2 - 0.5*log(ne) + log(Te); % KORC's model

CLogee = @(fb) 20*(1 - fb/2);

CLogee_f = @(g,ne,Te) log( (g-1)*sqrt(g+1)*lD(ne,Te)/(2*g*re) ); % Mosher's model

CLogee_b = @(g,Iz) log( (g-1)*sqrt(g+1)*me*c^2./(Iz*qe) );

CLogeZ = @(g,Iz,ne,Te,Z) log( lD(ne,Te)*(Iz*qe)./(Z*re*me*c^2) );

CLogeZ0 = @(g,Iz) log( (g^2 - 1)*me*c^2./(g*(Iz*qe)) );

CLogeH = @(g,ne,Te) log((g.^2-1)*lD(ne,Te)/(g*re));

% Models for E critical

E_CH = @(ne,Te) qe^3*ne*CLog(ne,Te)/(4*pi*ep0^2*me*c^2); % CH

Ec_Parks = @(ne,fb) qe^3*ne*CLogee(fb)/(4*pi*ep0^2*me*c^2); % Parks' model

Ec_Mosher = @(g,nef,Te,neb,Iz) qe^3*(nef*CLogee_f(g,nef,Te) + sum(neb.*CLogee_b(g,Iz)))/(4*pi*ep0^2*me*c^2); % Mosher model



% Models for Zeff

Zeff_CH = @(ni,Zi,ne) sum(ni.*Zi.^2)/ne; % CH

Zeff_Parks = @(ne,fb,Te,ni,Ai) CLog(ne,Te)*sum(ni.*Ai.^2)/(ne*CLogee(fb));% Parks' model

Zeff_Mosher = @(g,nef,neb,Te,nH,nz,Iz,Z,Ai) ( nH*CLogeH(g,nef,Te) + ...
        sum(nz.*(Z.^2.*CLogeZ(g,Iz,nef,Te,Z) + Ai.^2.*CLogeZ0(g,Iz))) )/( nef*CLogee_f(g,nef,Te) + sum(neb.*CLogee_b(g,Iz)) );% Parks' model


% Actual calculation of E critical and Zeff

nef = 0;
neb = 0;
neb_a = [];
Iz_a = [];

nH = 0;
ni = [];
nz = [];
Z = [];
A = [];
Zi = [];
nAi = [];
Ai = [];
for ii=1:numImpurities
    if isfield(iSpp,'H')
        nH = iSpp.H.nz;
    end
    
    nef = nef + sum(iSpp.(['spp' num2str(ii)]).nz.*iSpp.(['spp' num2str(ii)]).Z);
    neb = neb + iSpp.(['spp' num2str(ii)]).nn*iSpp.(['spp' num2str(ii)]).A + ...
        sum(iSpp.(['spp' num2str(ii)]).nz.*(iSpp.(['spp' num2str(ii)]).A-iSpp.(['spp' num2str(ii)]).Z));

    if ~isempty(iSpp.(['spp' num2str(ii)]).Iz)
        neb_a = [neb_a, iSpp.(['spp' num2str(ii)]).nz];
        Iz_a = [Iz_a, iSpp.(['spp' num2str(ii)]).Iz];
        nz = [nz,iSpp.(['spp' num2str(ii)]).nz];
        Z = [Z,iSpp.(['spp' num2str(ii)]).Z];
        A = [A, iSpp.(['spp' num2str(ii)]).A];
    end
   
    ni = [ni,iSpp.(['spp' num2str(ii)]).nz];
    Zi = [Zi,iSpp.(['spp' num2str(ii)]).Z];
    nAi = [nAi, sum(iSpp.(['spp' num2str(ii)]).nz)];
    Ai = [Ai, iSpp.(['spp' num2str(ii)]).A];
end

if isfield(iSpp,'H')
    nef = nef + nH;
    ni = [ni,nH];
    Zi = [Zi,1];
    nAi = [nAi,nH];
    Ai = [Ai,1];
end

fb = neb/(nef + neb);

Ec1 = zeros(1,numel(g));
Ec2 = zeros(1,numel(g));
Ec3 = zeros(1,numel(g));
Ec4 = zeros(1,numel(g));

Zeff1 = zeros(1,numel(g));
Zeff2 = zeros(1,numel(g));
Zeff3 = zeros(1,numel(g));
Zeff4 = zeros(1,numel(g));

Coeff1 = zeros(1,numel(g));
Coeff2 = zeros(1,numel(g));
Coeff3 = zeros(1,numel(g));
Coeff4 = zeros(1,numel(g));

for gg=1:numel(g)
    Ec1(gg) = E_CH(nef,Te);
    Ec2(gg) = E_CH(nef+neb,Te);
    Ec3(gg) = Ec_Parks(nef+neb,fb);
    Ec4(gg) = Ec_Mosher(g(gg),nef,Te,neb_a,Iz_a);
    

    Zeff1(gg) = Zeff_CH(ni,Zi,nef);
    Zeff2(gg) = Zeff_CH(ni,Zi,nef+neb);
    Zeff3(gg) = Zeff_Parks(nef+neb,fb,Te,nAi,Ai);
    Zeff4(gg) = Zeff_Mosher(g(gg),nef,neb_a,Te,nH,nz,Iz_a,Z,A);
    
    Coeff1(gg) = 2/( Ec1(gg)*(1 + Zeff1(gg)) );
    Coeff2(gg) = 2/( Ec2(gg)*(1 + Zeff2(gg)) );
    Coeff3(gg) = 2/( Ec3(gg)*(1 + Zeff3(gg)) );
    Coeff4(gg) = 2/( Ec4(gg)*(1 + Zeff4(gg)) );
end




EAxis = E/1E6;

figure
subplot(3,2,1)
plot(EAxis,Ec1,'k',EAxis,Ec2,'c',EAxis,Ec3,'r',EAxis,Ec4,'b')
grid minor;xlim([min(EAxis) max(EAxis)]);
legend({'CH ($n_{ef}$)','CH $(n_{ef} + n_{eb})$','Parks','Mosher'},'Interpreter','latex')
xlabel('Kinetic energy (MeV)','Interpreter','latex')
ylabel('$E_{c}$ (V/m)','Interpreter','latex')

subplot(3,2,3)
plot(EAxis,Zeff1,'k',EAxis,Zeff2,'c',EAxis,Zeff3,'r',EAxis,Zeff4,'b')
grid minor;xlim([min(EAxis) max(EAxis)]);
legend({'CH ($n_{ef}$)','CH $(n_{ef} + n_{eb})$','Parks','Mosher'},'Interpreter','latex')
xlabel('Kinetic energy (MeV)','Interpreter','latex')
ylabel('$\langle Z \rangle$','Interpreter','latex')

subplot(3,2,5)
plot(EAxis,Coeff1,'k',EAxis,Coeff2,'c',EAxis,Coeff3,'r',EAxis,Coeff4,'b')
grid minor;xlim([min(EAxis) max(EAxis)]);
legend({'CH ($n_{ef}$)','CH $(n_{ef} + n_{eb})$','Parks','Mosher'},'Interpreter','latex')
xlabel('Kinetic energy (MeV)','Interpreter','latex')
ylabel('$2/E_{c}(1+\langle Z \rangle)$','Interpreter','latex')

subplot(3,2,2)
plot(EAxis,Ec1./Ec1,'k',EAxis,Ec2./Ec1,'c',EAxis,Ec3./Ec1,'r',EAxis,Ec4./Ec1,'b')
grid minor;xlim([min(EAxis) max(EAxis)]);
legend({'CH ($n_{ef}$)','CH $(n_{ef} + n_{eb})$','Parks','Mosher'},'Interpreter','latex')
xlabel('Kinetic energy (MeV)','Interpreter','latex')
ylabel('$E_{c}/E_{CH}$ (V/m)','Interpreter','latex')

subplot(3,2,4)
plot(EAxis,Zeff1./Zeff1,'k',EAxis,Zeff2./Zeff1,'c',EAxis,Zeff3./Zeff1,'r',EAxis,Zeff4./Zeff1,'b')
grid minor;xlim([min(EAxis) max(EAxis)]);
legend({'CH ($n_{ef}$)','CH $(n_{ef} + n_{eb})$','Parks','Mosher'},'Interpreter','latex')
xlabel('Kinetic energy (MeV)','Interpreter','latex')
ylabel('$\langle Z \rangle/Z_{eff}$','Interpreter','latex')

subplot(3,2,6)
plot(EAxis,Coeff1./Coeff1,'k',EAxis,Coeff2./Coeff1,'c',EAxis,Coeff3./Coeff1,'r',EAxis,Coeff4./Coeff1,'b')
grid minor;xlim([min(EAxis) max(EAxis)]);
legend({'CH ($n_{ef}$)','CH $(n_{ef} + n_{eb})$','Parks','Mosher'},'Interpreter','latex')
xlabel('Kinetic energy (MeV)','Interpreter','latex')
ylabel('$E_{CH}(1+Z_{eff})/E_{c}(1+\langle Z \rangle)$','Interpreter','latex')
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
