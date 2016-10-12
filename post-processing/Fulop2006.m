% Rejection method for the momentum distribution function of Eq. (9) of
% Stahl et al. 2013. 

clear all
close all
clc

NE = 150;
NP = 200;

NPAR = 100;
NPERP = 50;

Npcls = 1E4;

pitch_min = 0; % in degrees
pitch_max = 80; % in degrees
Emax = 40E6; % Machimum energy in eV

% Plasma parameters and physical constants, all in SI units
kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
re = Kc*qe^2/(me*c^2);

% Parameters of fRE
ne = 3E20; % background electron density in m^-3
Zeff = 1.0; % Effective ion charge
Ec = 0.15; % Critical electric field in V/m
Epar = 10*Ec; % Parallel electric field in V/m
% Epar = 2.0; 
Ebar = Epar/Ec;
Tp = 10; % Background temperature in eV
Tp = Tp*qe; % in Joules (kB*T)
lambdaD = sqrt(ep*Tp/(ne*qe^2));
bmin = Zeff/(12*pi*ne*lambdaD^2);
Clog = log(lambdaD/bmin);
Tau = 1/(4*pi*re^2*ne*c*Clog);

pr = sqrt(Ebar - 1);
Er = c*sqrt((1 + pr^2)*(me*c)^2);
disp(['The minimum energy is: ' num2str(1E-6*Er/qe) ' MeV'])

Emax = Emax*qe; % In Joules
pmax = sqrt(Emax^2/(me*c^2)^2 - 1); % normalized to me*c

cz = sqrt( 3*(Zeff + 5)/pi )*Clog;
alpha = (Ebar - 1)/(1 + Zeff);

E = linspace(Er,Emax,NE);
p = sqrt( (E/c).^2 - (me*c)^2 );
p = p/(me*c);
pitch = (pi/180)*linspace(pitch_min,pitch_max,NP);
chi = cos(pitch);
chimin = min(chi);

% % % FRE(p,chi) % % %
fo = alpha/cz;
n = alpha*cz -1;
go = 97.8436;
eo = 15;
Vol = cz*alpha*(1 - exp(-eo))/eo;
C1 = 0.5*alpha;
C2 = 1/cz - C1;
C3 = 0.5*alpha/cz;

% % % Marginal distribution gRE % % %
% Pc =  @(x) fo*x./(C2*x.^2 + C1).^2; % comparison function
Co = @(x) fo./x;
Do = @(x) C2*x + C1./x;
P = @(x,a) Co(x).*(a./Do(x) + 1./Do(x).^2).*exp(-a.*Do(x)); % Marginal distribution function
Pc = @(y,b) fo*y.*exp(-y/cz)*(1 - exp(-b))/b;
chi_deviate = @(x) sqrt( fo*C1./(C2*(fo - 2*C2*C1*x)) - C1/C2 );

disp(['Normalization: ' num2str(trapz(fliplr(chi),P(fliplr(chi),0)))])

urnd = rand(1,Npcls);
deviate = chi_deviate(urnd);
IL = find(deviate < chimin);
while numel(IL) ~= 0
    urnd = rand(1,numel(IL));
    deviate(IL) = chi_deviate(urnd);
    IL = find(deviate < chimin);
end

h = histogram(deviate,'Normalization','pdf');
dchi = mean(diff(h.BinEdges));
xDeviate = 0.5*dchi + h.BinEdges(1:end-1);
yDeviate = h.Values;

figure
subplot(2,1,1)
% semilogy(chi,P(chi,0),'r--',chi,P(chi,pr),'b',xDeviate,yDeviate,'g.:')
semilogy(chi,P(chi,0),'r--',xDeviate,yDeviate,'b.:')
xlabel('$\chi$','Interpreter','latex')
ylabel('$\int f_{RE}(\chi,p) dp$','Interpreter','latex')
subplot(2,1,2)
plot(p,Pc(p,eo),'k')
xlabel('$p$ ($m_ec$)','Interpreter','latex')
ylabel('$\int f_c(\chi,p) d\chi$','Interpreter','latex')
% % % Marginal distribution gRE % % %
%% Bidimensional distribution

F = @(x,y) fo*y.*exp( -y.*(C2*x + C1./x) )./x;
% G = @(x,y) fo*y.*exp( -y/cz ).*exp(-go^2*(1-x).^2);
G = @(x,y) fo*y.*exp( -y/cz ).*exp(-eo*(1-x));

fRE = zeros(NE,NP);
fc = zeros(NE,NP);
for ii=1:NE
    for jj=1:NP
        fRE(ii,jj) = F(chi(jj),p(ii));
        fc(ii,jj) = G(chi(jj),p(ii));
    end
end
% % % FRE(p,chi) % % %


% Figures
figure
pitch = (180/pi)*pitch; % degrees
E = 1E-6*E/qe; % MeV

xAxis = pitch;

A = log10(fRE);
B = log10(fc);
levels = [0.2,0.1,0,-1,-2,-3,-4,-5,-6];
contour(xAxis,p,A,levels,'ShowText','on')
hold on;contour(xAxis,p,B,levels,'ShowText','on','LineColor',[0,0,0]);hold off
xlabel('$\chi$','Interpreter','latex')
ylabel('$p$ ($m_ec$)','Interpreter','latex')
box on;
colormap(jet)
hc = colorbar;
ylabel(hc,'$f_{RE}(\chi,p)$','Interpreter','latex','FontSize',16)
