function CO = MonteCarloCollisions(DT,numIt,np,Te,ne,Zeff)
% CO = MonteCarloCollisions(1E-2,5,1E4,10E3,1E19,1.0);
close all

CO.params = struct;
CO.params.kB = 1.38E-23; % Boltzmann constant
CO.params.Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
CO.params.mu0 = (4E-7)*pi; % Magnetic permeability
CO.params.ep = 8.854E-12;% Electric permitivity
CO.params.c = 2.9979E8; % Speed of light
CO.params.qe = 1.602176E-19; % Electron charge
CO.params.me = 9.109382E-31; % Electron mass

CO.b1 = [1,0,0];
CO.b2 = [0,1,0];
CO.b3 = [0,0,1];

CO.Bo = 2.19;
CO.B = [CO.Bo,0,0];
CO.E = [0,0,0];
CO.DT = DT;
CO.numIt = numIt;
CO.Te = Te;
CO.ne = ne;
CO.Zeff = Zeff;
CO.VTe = sqrt(CO.Te*CO.params.qe/CO.params.me);

CO.np = np; % Number of particles

CO = initializeCollisionOperators(CO);

CO = normalize(CO);

parametricStudy(CO);

Vo = ThermalDistribution(CO);
V = Vo;
% V = zeros(3,CO.np);
% V(1,:) = 0.1344;
% V(2,:) = 0.1344;

h = figure;
subplot(2,1,1)
hold on
histogram(V(1,:),'Normalization','pdf')
histogram(V(2,:),'Normalization','pdf')
histogram(V(3,:),'Normalization','pdf')
x = linspace(-1,1,100);
fx = exp(-0.5*x.^2/CO.VTe^2)/(CO.VTe*sqrt(2*pi));
hold on;plot(x,fx);hold off
hold off

for ii=1:CO.numIt
    for pp=1:CO.np
        V(:,pp) = collisionOperator(CO,V(:,pp)',CO.cop.dt);
    end
end

U = sqrt(sum(V.^2,1));

figure(h);
subplot(2,1,2)
hold on
histogram(V(1,:),'Normalization','pdf')
histogram(V(2,:),'Normalization','pdf')
histogram(V(3,:),'Normalization','pdf')
histogram(U,'Normalization','pdf')
hold off

end

function CO = initializeCollisionOperators(CO)
% NOTE: it is assumed that in the normalization vo and mo are the speed of
% light and the electron mass, respectively.
CO.cop = struct;

CO.cop.Te = CO.Te*CO.params.qe; % Background electron temperature in Joules

CO.cop.Ti = CO.Te; % Background ion temperature in Joules
CO.cop.ne = CO.ne; % Background electron density in 1/m^3
CO.cop.Zeff = CO.Zeff; % Full nuclear charge of each impurity: Z=1 for D, Z=10 for Ne
CO.cop.rD = ...
    sqrt( CO.params.ep*CO.cop.Te/(CO.cop.ne*CO.params.qe^2*(1 + CO.cop.Zeff*CO.cop.Te/CO.cop.Ti)) );
CO.cop.re = CO.params.qe^2/( 4*pi*CO.params.ep*CO.params.me*CO.params.c^2 );
CO.cop.Clog = 25.3 - 1.15*log10(1E-6*CO.cop.ne) + 2.3*log10(CO.cop.Te/CO.params.qe);
CO.cop.VTe = sqrt(2*CO.cop.Te/CO.params.me);
CO.cop.delta = CO.cop.VTe/CO.params.c;
CO.cop.Gamma = CO.cop.ne*CO.params.qe^4*CO.cop.Clog/(4*pi*CO.params.ep^2);
CO.cop.Tau = CO.params.me^2*CO.params.c^3/CO.cop.Gamma;
CO.cop.Ec = CO.cop.ne*CO.params.qe^3*CO.cop.Clog/(4*pi*CO.params.ep^2*CO.params.me*CO.params.c^2);
CO.cop.ED = CO.cop.ne*CO.params.qe^3*CO.cop.Clog/(4*pi*CO.params.ep^2*CO.cop.Te);

Ef = CO.E;
CO.cop.Vc = CO.cop.VTe*sqrt(0.5*CO.cop.Ec/sqrt(Ef*Ef'));

CO.cop.c = CO.params.c;

Ek = linspace(1.0,10E6,1E6)*CO.params.qe; % Kinetic energy
Er = CO.params.me*CO.params.c^2; % Rest energy
ET = Ek + Er; % Total energy
u = CO.params.c*sqrt(1 - Er^2./ET.^2);

% % % NORMALIZATIN PARAMETERS % % % 
CO.norm.B = CO.Bo;
CO.norm.q = abs(CO.params.qe);
CO.norm.m = CO.params.me;
CO.norm.t = CO.cop.Tau;
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

CO.cop.psi = @(v) 0.5*( erf(CO.cop.x(v)) - ...
    2*CO.cop.x(v).*exp(-CO.cop.x(v).^2)/sqrt(pi) )./CO.cop.x(v).^2;

CO.cop.CA = @(v) CO.cop.Gamma*CO.cop.psi(v)./v;

CO.cop.CF = @(v) CO.cop.Gamma*CO.cop.psi(v)/CO.cop.Te;

CO.cop.CB = @(v) (0.5*CO.cop.Gamma./v).*( CO.cop.Zeff + ...
    erf(CO.cop.x(v)) - CO.cop.psi(v) + 0.5*CO.cop.delta^4*CO.cop.x(v).^2 );

CO.cop.dvdp = @(v) 1./CO.cop.g(v).^3;

CO.cop.fun = @(v) 2*(1./CO.cop.x(v) + CO.cop.x(v)).*exp(-CO.cop.x(v).^2)/sqrt(pi) - ...
    erf(CO.cop.x(v))./CO.cop.x(v).^2 - CO.cop.psi(v);

CO.cop.f = @(v) ( 2*(1./CO.cop.x(v) + CO.cop.x(v)).*exp(-CO.cop.x(v).^2)/sqrt(pi) -...
    erf(CO.cop.x(v))./CO.cop.x(v).^2 - CO.cop.psi(v) )./(v.^2.*CO.cop.g(v).^3);

CO.cop.h = @(v) -CO.cop.CF(v) + ...
    2*CO.cop.CA(v)./(v.*CO.cop.g(v)) + CO.cop.Gamma*CO.cop.f(v);

Ek = 1E-6*Ek/CO.params.qe;
xAxis = u;

figure;
subplot(3,2,1);
plot(xAxis,CO.cop.CA(u));
ylabel('$C_A$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,2,2);
plot(xAxis,CO.cop.CB(u));
ylabel('$C_B$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,2,3);
plot(xAxis,CO.cop.CF(u));
ylabel('$C_F$ ($e B_0/m_e^2 c $)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,2,4);
semilogx(Ek,u);
ylabel('$v$ ($c$)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,2,[5 6]);
plot(xAxis,CO.cop.h(u),'b',xAxis,CO.cop.h(u)+CO.cop.CF(u),'r');
ylabel('Drag (blue)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
end

function CO = normalize(CO)

CO.B = CO.B/CO.norm.B;
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

f = @(U) exp(-0.5*(U/CO.VTe).^2); % VTe = sqrt(T/m)

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

function parametricStudy(CO)
Ek = linspace(1.0,10E6,1E6)*CO.params.qe; % Kinetic energy
Er = CO.params.me*CO.params.c^2; % Rest energy
ET = Ek + Er; % Total energy
u = CO.params.c*sqrt(1 - Er^2./ET.^2);

u = u/CO.norm.v; % Normalization
xi = 1.0;
dt = linspace(1E-4,0.1,1E3);

p = @(v) CO.cop.g(v).*v;

upperb1 = @(xi,v) -2*xi.*CO.cop.CB(v)./p(v).^2;
upperb2 = @(xi,v) -sqrt(2*CO.cop.CB(v).*(1 - xi.^2))./p(v); 


Ek = Ek/CO.params.qe;
Ek = Ek/1E6;

figure
subplot(2,1,1)
semilogy(u,upperb1(1.0,u),'b',u,upperb1(0.9,u),'r',u,upperb1(0.8,u),'g',...
    u,upperb1(0.7,u),'c')
subplot(2,1,2)
semilogy(u,upperb2(1.0,u),'b',u,upperb2(0.9,u),'r',u,upperb2(0.8,u),'g',...
    u,upperb2(0.7,u),'c')



end

function U = collisionOperator(CO,V,dt)
% This function calculates the random kicks due to collisions
U = zeros(1,3);

v = sqrt(V*V');
if v > 0
    g = CO.cop.g(v);
    p = g*v;
    
    U1 = g*(V*CO.b1');
    U2 = g*(V*CO.b2');
    U3 = g*(V*CO.b3');
    
    xi = U1/p;
    A = acos(xi);
    
    dWp = random('normal',0,dt);
    dWxi = random('uniform',0,1,1)*sqrt(dt);
    dWphi = random('uniform',0,2*pi,1)*sqrt(dt);
    
    CA = CO.cop.CA(v);
    CB = CO.cop.CB(v);
    
    dp = CO.cop.h(v)*dt;% + sqrt(2*CA)*dWp;
    
    dxi = -2*xi*CB*dt/p^2;% - sqrt(2*CB*(1 - xi^2))*dWxi/p;
    
    dA = (CB*cot(A)/p^2)*dt + sqrt(2*CB)*dWxi/p;
    dB = sqrt(2*CB)*dWphi/(p*sin(A));
    
    dphi = sqrt(2*CB)*dWphi/(p*sqrt(1 - xi^2));
    if ~isfinite(dphi)
        dphi = 0;
    end
    
%     f1 = @(dt) CO.cop.h(v)*dt + sqrt(2*CA)*sqrt(dt);
%     f2 = @(dt) CO.cop.h(v)*dt;
%     f3 = @(dt) -2*xi*CB*dt/p^2 - sqrt(2*CB*(1 - xi^2))*sqrt(dt)/p;
%     f4 = @(dt) -2*xi*CB*dt/p^2;
%     f5 = @(dt) sqrt(2*CB)*sqrt(dt)/(p*sqrt(1 - xi^2));
%     
%     t = linspace(1E-4,0.1,1E3);
%     figure
%     subplot(3,1,1)
%     plot(t,f1(t),'k',t,f2(t),'r')
%     subplot(3,1,2)
%     plot(t,f3(t),'k',t,f4(t),'r')
%     subplot(3,1,3)
%     plot(t,f5(t))
    
    pitch = acos(xi + dxi);
    
    dU1 = dp*cos(pitch);
    dU2 = dp*sin(pitch)*cos(dphi);
    dU3 = dp*sin(pitch)*sin(dphi);
    
    U(1) = U1 + dU1;
    U(2) = U2 + dU2;
    U(3) = U3 + dU3;
    
    U = U/sqrt( 1 + U*U' );
end
end
