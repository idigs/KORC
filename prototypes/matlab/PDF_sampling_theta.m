% Rejection method for the momentum distribution function of Eq. (9) of
% Stahl et al. 2013. 
clear all
close all
clc

NE = 300;
NP = 300;

Npcls = 1E6;

pitch_min = 0; % in degrees
pitch_max = 20; % in degrees
Emax = 60E6; % Machimum energy in eV

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
ne = 5E19; % background electron density in m^-3
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

pitch_max = pi*pitch_max/180;
pitch = linspace(-pitch_max,pitch_max,NP);

% % % FRE(p,chi) % % %
fo = alpha/cz;
n = alpha*cz -1;
eo = 15;
Vol = cz*alpha*(1 - exp(-eo))/eo;
C1 = 0.5*alpha;
C2 = 1/cz - C1;
C3 = 0.5*alpha/cz;

% Bidimensional PDFs
F = @(x,y) fo*y.*exp( -y.*(C2*x + C1./x) )./x;
% Bidimensional PDFs

fRE = zeros(NE,NP);
fc = zeros(NE,NP);
for ii=1:NE
    for jj=1:NP
        x = cos(pitch(jj));
        fRE(ii,jj) = F(x,p(ii));
    end
end
% % % FRE(p,chi) % % %
%% Sampling of distribution
Nsamples = 2*Npcls;

pitch_sampled = zeros(1,Nsamples);
p_sampled = zeros(1,Nsamples);

pitch_sampled(1) = 0.0;
[~,I] = max(fRE(:,1));
p_sampled(1) = p(I);

sigma_p = 1;
sigma_pitch = 0.1;

ii=2;
while (ii <= Nsamples)
    %
    pitch_test = pitch_sampled(ii-1) + random('norm',0,sigma_pitch);
    while (pitch_test > pitch_max) || (pitch_test < -pitch_max)
        pitch_test = pitch_sampled(ii-1) + random('norm',0,sigma_pitch);
    end
    
    chi_test = cos(pitch_test);
    p_test = p_sampled(ii-1) + random('norm',0,sigma_p);
    while (p_test < pr) || (p_test > pmax)
        p_test = p_sampled(ii-1) + random('norm',0,sigma_p);
    end
    
    %
%     pitch_test = random('unif',-pitch_max,pitch_max);
%     chi_test = cos(pitch_test);
%     p_test = random('unif',pr,pmax);
    
    %
    chi_sampled = cos(pitch_sampled(ii-1));
    ratio = F(chi_test,p_test)/F(chi_sampled,p_sampled(ii-1));
    
    if ( ratio >= 1.0 )
        pitch_sampled(ii) = pitch_test;
        p_sampled(ii) = p_test;
        
        ii  = ii + 1;
    elseif (random('uniform',0,1) < ratio)
        pitch_sampled(ii) = pitch_test;
        p_sampled(ii) = p_test;
        
        ii  = ii + 1;
    end
end
p_sampled(1:Npcls) = [];
pitch_sampled(1:Npcls) = [];

pitch_sampled(pitch_sampled<0)  = -pitch_sampled(pitch_sampled<0);

h = figure;
subplot(2,1,1)
histogram2((180/pi)*pitch_sampled,p_sampled,...
    'FaceColor','flat','Normalization','pdf','LineStyle','none',...
    'DisplayStyle','tile')
colormap(jet)
colorbar
view([0,90])
axis([0 20 pr pmax])
xlabel('$\theta$','Interpreter','latex')
ylabel('$p$ ($m_ec$)','Interpreter','latex')
%% Bidimensional distribution (plot)
% Figures
pitch = (180/pi)*pitch; % degrees
E = 1E-6*E/qe; % MeV

% A = log10(fRE);
A = fRE;
xAxis = pitch;

levels = linspace(0,max(max(A)),10);
subplot(2,1,2)
contourf(xAxis,p,A,levels,'ShowText','on')
xlabel('$\theta$','Interpreter','latex')
ylabel('$p$ ($m_ec$)','Interpreter','latex')
box on;
axis([0 20 pr pmax])
colormap(jet)
hc = colorbar;
ylabel(hc,'$f_{RE}(\chi,p)$','Interpreter','latex','FontSize',16)