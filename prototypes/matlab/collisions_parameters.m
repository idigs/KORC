function cp = collisions_parameters(Te,ne,nz,eta,E,B)
close all
% Collision parameters for RE in fusion plasmas
kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
e0 = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
% E = 5.0; % Background electric field
% B = 2.19 % Background magnetic field
eta = eta*pi/180; % Pitch angle

Te = Te*qe;
rD = sqrt( e0*Te/(ne*qe^2) );
re = qe^2/(4*pi*e0*me*c^2);

% Eo = linspace(0,1E2,1E4) + me*c^2/qe;
x = linspace(70,80,1E4);
% x = linspace(40,80,1E4);
% fx = (exp(x) - min(exp(x)))/max(exp(x) - min(exp(x)));
fx = exp(x)/max(exp(x));
Eo = 500E6*fx + me*c^2/qe;
gamma = qe*Eo/(me*c^2);
find(gamma==1,1,'last')
Eo(gamma==1) = [];
gamma(gamma==1) = [];
v = c*sqrt(1 - 1./gamma.^2);
Eo = (Eo -  + me*c^2/qe)/1E6;
% Iz = qe*[15.7596,27.62965,40.74,59.81,75.02]; % Argon
% Iz = qe*[21.5646,40.96296,63.45,97.12,126.21]; % Neon
% Zj = 1:1:numel(Iz);

% % % % Argon
Iz = qe*27.62965; %qe*(15.7596);
Zj = 2.0;
Zo = 18.0;

% % % % Neon
% Iz = qe*(21.5646);
% Zj = 1.0;
% Zo = 10.0;

nH = ne;
nef = ne + Zj*nz;
neb = (Zo - Zj)*nz;

ae = zeros(size(Eo));
ai = zeros(size(Eo));

aux = (gamma - 1).*sqrt(gamma + 1);
lratio = rD/re;

Clog_ef = log( 0.5*lratio*aux./gamma );
Clog_eb = log( aux*(me*c^2/Iz) );

ae = nef*Clog_ef + neb*Clog_eb;

figure
semilogx(Eo,Clog_ef,Eo,Clog_eb,'--')
xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
ylabel('$\log{(\Lambda_e)}$','Interpreter','latex','FontSize',16)
box on
grid on

Clog_eH = log( lratio*(gamma.^2 - 1)./gamma );
Clog_eZ = log( (lratio/Zj)*(Iz/(me*c^2)) );
Clog_eZo = log( (me*c^2/Iz)*(gamma.^2 - 1)./gamma );

ai = nH*Clog_eH + nz*(Clog_eZ*Zj^2 + Clog_eZo*Zo^2);

figure
semilogx(Eo,Clog_eH,Eo,Clog_eZ*ones(size(Eo)),'--',Eo,Clog_eZo,':')
xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
ylabel('$\log{(\Lambda_i)}$','Interpreter','latex','FontSize',16)
box on
grid on


% Collision force magnitude
Fcolle = 4*pi*me*c^4*re^2*ae.*gamma.*(gamma + 1)./(v.*gamma).^2;
Fcollef = 4*pi*me*c^4*re^2*(nef*Clog_ef).*gamma.*(gamma + 1)./(v.*gamma).^2;
Fcolleb = 4*pi*me*c^4*re^2*(neb*Clog_eb).*gamma.*(gamma + 1)./(v.*gamma).^2;
Fcolli = 4*pi*me*c^4*re^2*ai./(v.^2.*gamma);
FcollH = 4*pi*me*c^4*re^2*(nH*Clog_eH)./(v.^2.*gamma);
FcollZj = 4*pi*me*c^4*re^2*(nz*Clog_eZ*Zj^2)./(v.^2.*gamma);
FcollZo = 4*pi*me*c^4*re^2*(nz*Clog_eZo*Zo^2)./(v.^2.*gamma);
FE = E*qe*ones(size(FcollZo));

A2 = (v*cos(eta)).^2*(E/c)^4;
B2 = v*B^2*sin(eta);
C2 = (gamma/c).^4.*( v.*((E*v*cos(eta)/c).^2 - E^2 - (v*B*sin(eta)).^2) ).^2;
AC = (gamma/c).^2.*(E*v*cos(eta)/c).^2.*((E*v*cos(eta)/c).^2 - E^2 - (v*B*sin(eta)).^2);

FR = (2*qe^2*re/(3*me*c))*sqrt( A2 + B2 + C2 + 2*AC );

figure
loglog(Eo,Fcolle,'r-',Eo,Fcollef,'r:',Eo,Fcolleb,'r--',...
    Eo,Fcolli,'k-',Eo,FcollH,'k:',Eo,FcollZj,'k--',Eo,FcollZo,'k-.',...
    Eo,FE,'b-',Eo,FR,'g-')
xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
ylabel('$F_{coll}$ (N)','Interpreter','latex','FontSize',16)
box on
grid on
% Collision force magnitude


% Power balance
Pcolle = 4*pi*me*c^4*re^2*ae.*(gamma + 1).*gamma.^2./(v.*gamma.^3);
Pcollef = 4*pi*me*c^4*re^2*(nef*Clog_ef).*(gamma + 1).*gamma.^2./(v.*gamma.^3);
Pcolleb = 4*pi*me*c^4*re^2*(neb*Clog_eb).*(gamma + 1).*gamma.^2./(v.*gamma.^3);

Pcolli = 4*pi*me*c^4*re^2*ai./(v.*gamma);
PcollH = 4*pi*me*c^4*re^2*(nH*Clog_eH)./(v.*gamma);
PcollZj = 4*pi*me*c^4*re^2*(nz*Clog_eZ*Zj^2)./(v.*gamma);
PcollZo = 4*pi*me*c^4*re^2*(nz*Clog_eZo*Zo^2)./(v.*gamma);

PE = qe*E*v*cos(eta);

A2 = (v*cos(eta)).^2*(E/c)^4;
B2 = v*B^2*sin(eta);
C2 = (gamma/c).^4.*( v.*((E*v*cos(eta)/c).^2 - E^2 - (v*B*sin(eta)).^2) ).^2;
AC = (gamma/c).^2.*(E*v*cos(eta)/c).^2.*((E*v*cos(eta)/c).^2 - E^2 - (v*B*sin(eta)).^2);

PR = (2*qe^2*re/(3*me*c))*( E^2 + gamma.^2.*((E*v*cos(eta)/c).^2 - E^2 - (v*B*sin(eta)).^2) );
PR = abs(PR);

figure
loglog(Eo,Pcolle,'r-',Eo,Pcollef,'r:',Eo,Pcolleb,'r--',...
    Eo,Pcolli,'k-',Eo,PcollH,'k:',Eo,PcollZj,'k--',Eo,PcollZo,'k-.',...
    Eo,PE,'b-',Eo,PR,'g-')
xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
ylabel('$P$ (Watts)','Interpreter','latex','FontSize',16)
box on
grid on
% Power balance

ae = ae/1E20;
ai = ai/1E20;

Zeff = (nH + nz*Zj^2)/(nH + nz*Zj);

figure
subplot(2,1,1)
semilogx(Eo,ae,'r',Eo,ai,'b')
xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
ylabel('$\alpha_{e,i}$','Interpreter','latex','FontSize',16)
box on
grid on
subplot(2,1,2)
semilogx(Eo,ai./ae,Eo,Zeff*ones(size(Eo)))
xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
ylabel('$Z_{coll}, Z_{eff}$','Interpreter','latex','FontSize',16)
box on
grid on

end