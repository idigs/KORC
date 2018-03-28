% Script for studying Coulomb logarithms of fast electrons collisioning
% with bound and free electrons as well as with ions.

function CLog = CoulombLogarithms(Te,E)
% Te: background electron temperature in eV.
% E: test electron's energy in eV.
% iZ: cell array with the inpurity variables.
% nZ: density of neutral impurities.
% nZ: density of impurities with mean ionization level Z.
% Iz: ionization energy of different impurities in eV.
% Z: mean ionization level of impurities.

% Deuterium:
% Z = [0]; % Corresponding ionization level of D.
% Iz = [15.46658]; % Ionization energy of unique electron in eV.
% Argon:
% Z = [0,1,2]; % Corresponding ionization level of Ar.
% Iz = [15.7596117,27.62967,40.735]; % Ionization energy of different Ar impurities.

iZ = loadStruct();

CLog = struct;

kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep0 = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
re = qe^2/(4*pi*ep0*me*c^2);

% Unit conversion
TeJ = Te*qe; % Converted to Joules

EJ = E*qe; % Converted to Joules
g = EJ/(me*c^2);

iZ.Ar.InJ = iZ.Ar.In*qe; % Converted to Joules
iZ.Ar.IzJ = iZ.Ar.Iz*qe; % Converted to Joules
iZ.D.InJ = iZ.D.In*qe; % Converted to Joules

ne = iZ.D.nz + sum(iZ.Ar.nz.*iZ.Ar.Z);

lD = @(ne,Te) sqrt(ep0*(Te*qe)/(qe^2*ne));

CLogee_f1 = @(ne,Te) 32.2 - 0.5*log(ne) + log(Te); % KORC's model

CLogee_f2 = @(g,ne,Te) log( (g-1)*sqrt(g+1)*lD(ne,Te)/(2*g*re) ); % Mosher's model

CLogee_b = @(g,Iz) log( (g-1)*sqrt(g+1)*me*c^2./(Iz*qe) );

CLogeZ = @(g,Iz,ne,Te,Z) log( lD(ne,Te)*(Iz*qe)./(Z*re*me*c^2) );

CLogeZ0 = @(g,Iz) log( (g^2 - 1)*me*c^2./(g*(Iz*qe)) );

Zeff_simple = (iZ.D.nz + sum(iZ.Ar.nz.*iZ.Ar.Z.^2))/ne;

% A = iZ.D.nz*iZ.D.Z^2*CLogee_f1(ne,Te) + iZ.Ar.A^2*sum(iZ.Ar.nz.*CLogeZ0(g,iZ.Ar.Iz)) +...
%     sum(iZ.Ar.nz.*iZ.Ar.Z.^2.*CLogeZ(g,iZ.Ar.Iz,ne,Te,iZ.Ar.Z));
A = iZ.D.nz*iZ.D.Z^2*CLogee_f1(ne,Te) + sum(iZ.Ar.nz.*iZ.Ar.Z.^2)*CLogee_f1(ne,Te);

B = ne*CLogee_f1(ne,Te) + iZ.D.nn*CLogee_b(g,iZ.D.In) +...
    iZ.Ar.nn*CLogee_b(g,iZ.Ar.In) + sum(iZ.Ar.nz.*CLogee_b(g,iZ.Ar.Iz));

Zeff = A/B;

% Final output

CLog.CLogee_f1 = CLogee_f1;
CLog.CLogee_f2 = CLogee_f2;
CLog.CLogee_b = CLogee_b;
CLog.CLogeZ = CLogeZ;
CLog.CLogeZ0 = CLogeZ0;

CLog.description = struct;
CLog.description.CLogee_f1 = 'Collisions with free electrons: KORC model.';
CLog.description.CLogee_f12 = 'Collisions with free electrons: Mosher model.';
CLog.description.CLogee_b = 'Collisions with bound electrons: Mosher model.';
CLog.description.CLogeZ = 'Collisions with partially ionized impurities. The chosen Z must correspond to Iz.';
CLog.description.CLogeZ = 'Collisions with nuclei of partially ionized impurities. Iz is the same as for CLogeZ.';
end


function iZ = loadStruct()
iZ = struct;

iZ.Ar.A = 18; % Argon
iZ.Ar.nn = 0.002611*1E13*1E6; % neutral density in m^-3
iZ.Ar.In = 15.7596117; % Ionization energy of neutral Argon.
iZ.Ar.nz = [0.284,1.17]*1E13*1E6; % density of partially ionized impurities.
iZ.Ar.Z = [1,2]; % Average charge state.
iZ.Ar.Iz = [27.62967,40.735]; % Ionization energy.

iZ.D.A = 1; % Deuterium
iZ.D.nn = 0.04481*1E13*1E6; % neutral density in m^-3
iZ.D.In = 15.46658; % Ionization energy of neutral Argon.
iZ.D.nz = 3.764*1E13*1E6; % Ionized deuterium.
iZ.D.Z = 1;
end

