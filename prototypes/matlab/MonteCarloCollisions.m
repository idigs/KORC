function CO = MonteCarloCollisions(DT,numIt,np,Te,ne,Zeff,ERE)
% CO = MonteCarloCollisions(1E-2,5,1E4,10E3,1E19,1.0);
% CO = MonteCarloCollisions(1E-3,10000,100,1E3,1E20,0.0);
close all

CO.params = struct;
CO.params.kB = 1.38E-23; % Boltzmann constant
CO.params.Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
CO.params.mu0 = (4E-7)*pi; % Magnetic permeability
CO.params.ep = 8.854E-12;% Electric permitivity
CO.params.c = 2.9979E8; % Speed of light
CO.params.qe = -1.602176E-19; % Electron charge
CO.params.me = 9.109382E-31; % Electron mass


CO.b1 = [1,0,0];
CO.b2 = [0,1,0];
CO.b3 = [0,0,1];

CO.Bo = 0.0; % 2.19 ; % in Teslas
CO.Eo = 0.0; % in V/m
CO.B = [CO.Bo,0,0];
CO.E = [CO.Eo,0,0];
CO.DT = DT;
CO.numIt = numIt;
CO.Te = Te;
CO.ne = ne;
CO.Zeff = Zeff;
CO.VTe = sqrt(2*CO.Te*abs(CO.params.qe)/CO.params.me);

CO.np = np; % Number of particles

CO = initializeCollisionOperators(CO);

partiallyIonisedCollisionOperator(CO);

CO = normalize(CO);

% V = ThermalDistribution(CO);

ERE = CO.params.me*CO.params.c^2 + ERE *abs(CO.params.qe);
u = sqrt(1 - (CO.params.me*CO.params.c^2./(ERE)).^2);

V = repmat([u;0;0],[1,CO.np]);

if (CO.Bo > 0)
    CO.params.wc = sqrt(1 - u^2)*abs(CO.params.qe)*CO.Bo/CO.params.me;
    disp(['Gyro-period: ' num2str(2*pi/CO.params.wc) ' s'])
    
    dt = DT*(2*pi/CO.params.wc);
    dt = dt/CO.norm.t;
else
    dt = CO.cop.dt;
end

x = linspace(-1,1,500);
fx = exp(-x.^2/CO.VTe^2)/(CO.VTe*sqrt(pi));

snapshot = floor(numIt/5);
sh = figure;

a = -dt;

tic 

for ii=1:CO.numIt
    
    for pp=1:CO.np
        if (CO.Bo > 0)
            U = V(:,pp)/sqrt(1 - V(:,pp)'*V(:,pp));
            
            Up = [U(1) + a*(CO.Eo - 0.5*V(3,pp)*CO.Bo);
                U(2) + 0.5*a*V(3,pp)*CO.Bo;
                U(3) - 0.5*a*V(2,pp)*CO.Bo];
            
            tau = 0.5*a*CO.Bo;
            gp = sqrt(1 + Up'*Up);
            sig = gp^2 - tau^2;
            us = 0.5*a*CO.Bo*Up(1);
            g = sqrt(0.5)*sqrt( sig + sqrt(sig^2 + 4*(tau^2 + us^2)) );
            t = tau/g;
            s = 1/(1 + t^2);
            
            U = s*[Up(1) + Up(1)*t^2 - Up(3)*t;
                Up(2) + Up(3)*t;
                Up(3) - Up(2)*t];
            
            V(:,pp) = U/sqrt(1 + U'*U);
        end
        
        V(:,pp) = collisionOperator(CO,V(:,pp)',dt);
    end
    
    
    if mod(ii,snapshot) == 0
        figure(sh);
        subplot(3,5,ii/snapshot)
        hold on;plot3(V(1,:),V(2,:),V(3,:),'r.');
        grid on;box on;axis equal;axis([-1,1,-1,1,-1,1]);hold off
        view([150,10])
        title(['$t=$' num2str(ii*dt*CO.norm.t) ' s'],'interpreter','latex')
        
        subplot(3,5,ii/snapshot+5)
        hold on
        histogram(V(1,:),25,'Normalization','pdf','LineStyle','none')
        histogram(V(2,:),25,'Normalization','pdf','LineStyle','none')
        histogram(V(3,:),25,'Normalization','pdf','LineStyle','none')
        plot(x,fx,'k','LineWidth',2)
        hold off
        legend({'$f(v_x)$','$f(v_y)$','$f(v_z)$','$f_M(v)$'},'Interpreter','latex')
        box on
        grid on
        
        subplot(3,5,ii/snapshot+10)
        g = 1./sqrt(1-sum(V.^2,1));
%         E = (g*CO.params.me*CO.params.c^2 - CO.params.me*CO.params.c^2)/(abs(CO.params.qe)*1E6);
        E = (g*CO.params.me*CO.params.c^2 - CO.params.me*CO.params.c^2);
        kBTe = CO.Te*abs(CO.params.qe);
        EAxis = linspace(0,10*kBTe,500);
        fE = 2*(1/kBTe)^(1.5)*sqrt(EAxis/pi).*exp(-EAxis/kBTe);
%         EAxis = EAxis/(abs(CO.params.qe)*1E6);
        hold on
        histogram(E,25,'Normalization','pdf','LineStyle','none')
        plot(EAxis,fE,'k','LineWidth',2)
        hold off
        xlabel('$\mathcal{E}$ (MeV)','interpreter','latex')
        box on
        grid on
    end
end

toc

CO.V = V;
end

function CO = initializeCollisionOperators(CO)
% NOTE: it is assumed that in the normalization vo and mo are the speed of
% light and the electron mass, respectively.
CO.cop = struct;

CO.cop.Te = CO.Te*abs(CO.params.qe); % Background electron temperature in Joules

CO.cop.Ti = CO.Te; % Background ion temperature in Joules
CO.cop.ne = CO.ne; % Background electron density in 1/m^3
CO.cop.Zeff = CO.Zeff; % Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
CO.cop.rD = ...
    sqrt( CO.params.ep*CO.cop.Te/(CO.cop.ne*abs(CO.params.qe)^2*(1 + CO.cop.Zeff*CO.cop.Te/CO.cop.Ti)) );
CO.cop.re = abs(CO.params.qe)^2/( 4*pi*CO.params.ep*CO.params.me*CO.params.c^2 );
CO.cop.Clog = 25.3 - 1.15*log10(1E-6*CO.cop.ne) + 2.3*log10(CO.cop.Te/abs(CO.params.qe));
CO.cop.VTe = CO.VTe;
CO.cop.delta = CO.cop.VTe/CO.params.c;
CO.cop.Gamma = CO.cop.ne*abs(CO.params.qe)^4*CO.cop.Clog/(4*pi*CO.params.ep^2);

CO.cop.Tauc = CO.params.me^2*CO.VTe^3/CO.cop.Gamma;

CO.cop.Tau = CO.params.me^2*CO.params.c^3/CO.cop.Gamma;
% CO.cop.Tau = CO.cop.Tauc;

CO.cop.Ec = CO.cop.ne*abs(CO.params.qe)^3*CO.cop.Clog/(4*pi*CO.params.ep^2*CO.params.me*CO.params.c^2);
CO.cop.ED = CO.cop.ne*abs(CO.params.qe)^3*CO.cop.Clog/(4*pi*CO.params.ep^2*CO.cop.Te);

Ef = CO.E;
CO.cop.Vc = CO.cop.VTe*sqrt(0.5*CO.cop.Ec/sqrt(Ef*Ef'));

CO.cop.c = CO.params.c;

% g = linspace(1,5,1000); % relativistic gamma factor
% Er = CO.params.me*CO.params.c^2; % Rest energy
% E = g*CO.params.me*CO.params.c^2;%*linspace(6E5,50E6,200)*ST.params.qe;

Er = CO.params.me*CO.params.c^2; % Rest energy
E = Er + linspace(1.0E-3,50.0E6,1E5)*abs(CO.params.qe);
g = E/Er;

u = CO.params.c*sqrt(1 - Er^2./E.^2);


disp(['Relativistic collisional time: ' num2str(CO.cop.Tau) ' s'])
disp(['Thermal collisional time: ' num2str(CO.cop.Tauc) ' s'])


% % % NORMALIZATIN PARAMETERS % % % 
CO.norm.B = CO.Bo;
CO.norm.q = abs(CO.params.qe);
CO.norm.m = CO.params.me;
if (CO.Bo > 0)
    CO.norm.wc = abs(CO.params.qe)*CO.Bo/CO.params.me;
    CO.norm.t = 1/CO.norm.wc;
else
    CO.norm.t = CO.cop.Tau;
end
CO.norm.l = CO.params.c*CO.norm.t;
CO.norm.v = CO.params.c;
CO.norm.E = CO.norm.m*CO.norm.v^2;
CO.norm.p = CO.norm.m*CO.norm.v;
CO.norm.n = 1/CO.norm.l^3;
% % % NORMALIZATIN PARAMETERS % % % 

% % % NORMALIZATION % % %
CO.cop.Te = CO.cop.Te/(CO.norm.m*CO.norm.v^2);
CO.cop.Ti = CO.cop.Ti/(CO.norm.m*CO.norm.v^2);
CO.cop.ne = CO.cop.ne/CO.norm.l^3;
CO.cop.VTe = CO.cop.VTe/CO.norm.v;

CO.cop.Gamma = CO.cop.Gamma*(CO.norm.t/(CO.norm.m^2*CO.norm.v^3));
CO.cop.Tau = CO.cop.Tau/CO.norm.t;
CO.cop.Vc = CO.cop.Vc/CO.norm.v;

CO.cop.c = CO.cop.c/CO.norm.v;

u = u/CO.norm.v;
% % % NORMALIZATION % % %

CO.cop.dt = CO.DT*CO.cop.Tau; % % % % TIME STEP % % % %

CO.cop.g = @(v) 1./sqrt(1 - v.^2);

CO.cop.x = @(v) v/CO.cop.VTe;

CO.cop.psi = @(v) 0.5*( erf(CO.cop.x(v)) - 2*CO.cop.x(v).*exp(-CO.cop.x(v).^2)/sqrt(pi) )./CO.cop.x(v).^2;

CO.cop.CA = @(v) CO.cop.Gamma*CO.cop.psi(v)./v;

CO.cop.CF = @(v) CO.cop.Gamma*CO.cop.psi(v)/CO.cop.Te;

CO.cop.CB = @(v) (0.5*CO.cop.Gamma./v).*( CO.cop.Zeff + ...
    erf(CO.cop.x(v)) - CO.cop.psi(v) + 0.5*CO.cop.delta^4*CO.cop.x(v).^2 );

% CO.cop.CB = @(v) (0.5*CO.cop.Gamma./v).*( erf(CO.cop.x(v)) - CO.cop.psi(v) );

CO.cop.vD = @(v) 2*CO.cop.CB(v)./(v.*CO.cop.g(v)).^2;

CO.cop.vpar = @(v) 2*CO.cop.CA(v)./(v.*CO.cop.g(v)).^2;

CO.cop.vS = @(v) 2*CO.cop.CF(v)./(v.*CO.cop.g(v));

E = 1E-6*(E-Er)/abs(CO.params.qe);
xAxis = E;

figure;
subplot(2,1,1);
plot(xAxis,CO.cop.CA(u),'k',xAxis,CO.cop.CB(u),'r',xAxis,CO.cop.CF(u),'b');
legend({'$C_A$ ($e B_0/m_e^3 c^2 $)','$C_B$ ($e B_0/m_e^3 c^2 $)','$C_F$ ($e B_0/m_e^2 c $)'},...
    'Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
box on

subplot(2,1,2)
semilogy(xAxis,CO.cop.vS(u)/CO.norm.t,'k',...
    xAxis,CO.cop.vD(u)/CO.norm.t,'r',...
    xAxis,CO.cop.vpar(u)/CO.norm.t,'b');
legend({'$\nu_S$','$\nu_D$','$\nu_\parallel$'},'Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
ylabel('$\nu$ (s$^{-1}$)','Interpreter','latex')
box on


% ylabel('$C_A$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
% xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
% subplot(3,1,2);
% plot(xAxis,CO.cop.CB(u));
% ylabel('$C_B$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
% xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
% subplot(3,1,3);
% plot(xAxis,CO.cop.CF(u));
% ylabel('$C_F$ ($e B_0/m_e^2 c $)','Interpreter','latex')
% xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
end

function PICO = partiallyIonisedCollisionOperator(CO)
c = CO.params.c;
ne = CO.ne/(1E20);
Te = CO.Te/1E3;
VTe = CO.VTe;

ZAr = 18; % Argon atomic number
IzAr = [15.7596,27.62965,40.74,59.81,75.02]; % Ionisation energy of Argon (eV)
aAr = [0.353,0.329,0.306,0.283,0.260,0.238];
nAr = CO.ne*repmat(1.0,1,5);

IzNe = [21.5646,40.96296,63.45,97.12,126.21]; % Ionisation energy of Neon (eV)

k = 5;
alpha = 1/137;

Clog0 = 14.9 - 0.5*log(ne) + log(Te);

% functions

Clog_ei = @(g) Clog0 + log(1 + (2*c*sqrt(g.^2 - 1)/VTe).^k)/k;

Clog_ee = @(g) Clog0 + log(1 + (2*c^2*(g - 1)/VTe^2).^(k/2))/k;

% yj = @(p,aj) 2*aj.*p/alpha;
yj = @(g,aj) 2*aj.*sqrt(g.^2 - 1)/alpha;

Gj = @(g,Zj,Zoj,aj) (2/3)*(Zj^2 - Zoj^2).*log(yj(g,aj).^1.5 + 1) -...
                    (2/3)*(Zj-Zoj).^2.*yj(g,aj).^1.5./(3*yj(g,aj).^1.5 + 1);

Tau_ei = @(g) 4*pi*CO.params.ep^2*CO.params.me^2*CO.params.c^3./(CO.ne*CO.params.qe^4*Clog_ei(g));

Tau_ee = @(g) 4*pi*CO.params.ep^2*CO.params.me^2*CO.params.c^3./(CO.ne*CO.params.qe^4*Clog_ee(g));

Zeff = @(nj,Zoj) nj.*Zoj.^2/CO.ne;

hj = @(g,Ij) sqrt((g.^2 - 1).*(g - 1))./Ij;

vDCS = @(g,nj,Zoj) (Zeff(nj,Zoj)./Tau_ei(g)).*g./(g.^2 - 1).^1.5;

vSCS = @(g) (1./Tau_ei(g)).*g.^2./(g.^2 - 1).^1.5;

vD = @(g,nj,Zj,Zoj,aj) vDCS(g,nj,Zoj).*(1 + ( nj.*Gj(g,Zj,Zoj,aj)./(CO.ne*Clog_ei(g)) )./Zeff(nj,Zoj));

vS = @(g,nj,Zj,Zoj,Ij) vSCS(g).*( 1 + nj.*(Zj-Zoj).*(log(1 + hj(g,Ij).^k)/k - 1 + 1./g.^2)./(CO.ne*Clog_ee(g)) );


Er = CO.params.me*CO.params.c^2; % Rest energy (Joules)
E = Er + linspace(1.0E-3,50.0E6,1E5)*abs(CO.params.qe);
g = E/Er;
p = sqrt((E/c).^2 - (CO.params.me*CO.params.c)^2);
p = p/(CO.params.me*CO.params.c);

EAxis = E/(abs(CO.params.qe*1E6));

xAxis = p;

for ii=1:5
    figure
    subplot(2,2,1)
    plot(xAxis,Clog_ei(g),'r',xAxis,Clog_ee(g),'k')
    legend({'$\log\Lambda_{ei}$','$\log\Lambda_{ee}$'},'Interpreter','latex')
    ylabel('$\log\Lambda$','Interpreter','latex')
    xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
    box on;xlim([min(xAxis) max(xAxis)])
    
    subplot(2,2,2)
    plot(xAxis,Tau_ei(g),'r',xAxis,Tau_ee(g),'k')
    legend({'$\tau_{ei}$','$\tau_{ee}$'},'Interpreter','latex')
    ylabel('$\tau_c$ (s)','Interpreter','latex')
    xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
    box on;xlim([min(xAxis) max(xAxis)])
    
    subplot(2,2,3)
    semilogy(xAxis,vDCS(g,nAr(ii),ii),'r',...
        xAxis,vD(g,nAr(ii),ZAr,ii,aAr(ii)),'r--',...
        xAxis,vSCS(g),'k',...
        xAxis,vS(g,nAr(ii),ZAr,ii,IzAr(ii)*abs(CO.params.qe)/Er),'k--')
    legend({'$\nu_{D,CS}$','$\nu_D$','$\nu_{S,CS}$','$\nu_S$'},'Interpreter','latex')
    ylabel('$\nu_{CS}$ (s$^-1$)','Interpreter','latex')
    xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
    box on;xlim([min(xAxis) max(xAxis)])
    
    subplot(2,2,4)
    semilogy(xAxis,vD(g,nAr(ii),ZAr,ii,aAr(ii))./vDCS(g,nAr(ii),ii),'r',...
        xAxis,vS(g,nAr(ii),ZAr,ii,IzAr(ii)*abs(CO.params.qe)/Er)./vSCS(g),'k')
    legend({'$\nu_D/\nu_{D,CS}$','$\nu_S/\nu_{S,CS}$'},'Interpreter','latex')
    ylabel('$\nu/\nu_{CS}$','Interpreter','latex')
    xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
    box on;xlim([min(xAxis) max(xAxis)])
end

end

function CO = normalize(CO)

CO.B = CO.B/CO.norm.B;
CO.Bo = CO.Bo/CO.norm.B;
CO.E = CO.E/CO.norm.E;
CO.Eo = CO.Eo/CO.norm.E;
CO.VTe = CO.VTe/CO.norm.v;

end

function V = SampleThermalDistribution(CO)
vmax = 1.0; % A fraction of the speed of light
sv = CO.VTe/10;

f = @(U) exp(-0.5*sum(U.^2)/CO.VTe^2); % VTe = sqrt(T/m)
% f = @(U) exp(-sum(U.^2)/CO.VTe^2); % VTe = sqrt(2T/m)

V = zeros(3,CO.np);
V(:,1) = CO.VTe;

ii=2;
while (ii <= CO.np)
    U = V(:,ii-1) + random('norm',0,sv,3,1);
    while any(U > vmax)
        U = V(:,ii-1) + random('norm',0,sv,3,1);
    end
    ratio = f(U)/f(V(:,ii-1));
    
    if ( ratio >= 1.0 )
        V(:,ii) = U;
        ii  = ii + 1;
    elseif (ratio > random('uniform',0,1))
        V(:,ii) = U;
        ii  = ii + 1;
    end
end

end

function V = ThermalDistribution(CO)
vmax = 1.0; % A fraction of the speed of light
sv = CO.VTe/5;

f = @(U) exp(-(U/CO.VTe).^2); % VTe = sqrt(2T/m)

V = zeros(3,CO.np);
V(:,1) = CO.VTe;

for jj=1:3
    ii=2;
    while (ii <= CO.np)
        U = V(jj,ii-1) + random('norm',0,sv);
        ratio = f(U)/f(V(jj,ii-1));
        
        if ( ratio >= 1.0 )
            V(jj,ii) = U;
            ii  = ii + 1;
        elseif (ratio > random('uniform',0,1))
            V(jj,ii) = U;
            ii  = ii + 1;
        end
    end
end

end

function U = collisionOperator(CO,V,dt)
% This function calculates the random kicks due to collisions
U = zeros(1,3);

v = sqrt(V*V');
g = CO.cop.g(v);

v1 = V/v;
v2 = cross(v1,CO.b2);
v2 = v2/sqrt(v2*v2');
v3 = cross(v1,v2);
v3 = v3/sqrt(v3*v3');

U1 = g*(V*v1');
U2 = g*(V*v2');
U3 = g*(V*v3');

% dW = random('norm',0,sqrt(dt),[3,1]);
dW = sqrt(3*dt)*random('unif',-1,1,[3,1]);

CA = CO.cop.CA(v);
CB = CO.cop.CB(v);
CF = CO.cop.CF(v);

dU = [-2.0*CF*dt + sqrt(2*CA)*dW(1);
    sqrt(2.0*CB)*dW(2);
    sqrt(2.0*CB)*dW(3)];

U(1) = (U1+dU(1))*(v1*CO.b1') + (U2+dU(2))*(v2*CO.b1') + (U3+dU(3))*(v3*CO.b1');
U(2) = (U1+dU(1))*(v1*CO.b2') + (U2+dU(2))*(v2*CO.b2') + (U3+dU(3))*(v3*CO.b2');
U(3) = (U1+dU(1))*(v1*CO.b3') + (U2+dU(2))*(v2*CO.b3') + (U3+dU(3))*(v3*CO.b3');

U = U/sqrt( 1 + U*U' );

end
