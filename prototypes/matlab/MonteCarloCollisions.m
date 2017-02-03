function CO = MonteCarloCollisions(DT,numIt,np,Te,ne,Zeff)
% CO = MonteCarloCollisions(5E-2,100,1000.0,1E19,1.0);
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

V = SampleThermalDistribution(CO);

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

figure(h);
subplot(2,1,2)
hold on
histogram(V(1,:),'Normalization','pdf')
histogram(V(2,:),'Normalization','pdf')
histogram(V(3,:),'Normalization','pdf')
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

Ek = linspace(1.0E3,1E6,1000)*CO.params.qe; % Kinetic energy
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
Ek = u;

figure;
subplot(3,2,1);
plot(Ek,CO.cop.CA(u));
ylabel('$C_A$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,2,2);
plot(Ek,CO.cop.CB(u));
ylabel('$C_B$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,2,3);
plot(Ek,CO.cop.CF(u));
ylabel('$C_F$ ($e B_0/m_e^2 c $)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,2,4);
plot(Ek,CO.cop.g(u));
ylabel('$f(v)$','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,2,[5 6]);
plot(Ek,CO.cop.h(u),'b',Ek,CO.cop.h(u)+CO.cop.CF(u),'r');
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
    
    R = random('normal',0,dt,1,3);   
    dWp = R(1);
    dWxi = R(2);
    dWphi = R(3);
%     dWxi = random('uniform',0,1,1)*sqrt(dt);
%     dWphi = random('uniform',0,2*pi,1)*sqrt(dt);
    
    CA = CO.cop.CA(v);
    CB = CO.cop.CB(v);
    
    dp = CO.cop.h(v)*dt + sqrt(2*CA)*dWp;
    
    dxi = -2*xi*CB*dt/p^2 - sqrt(2*CB*(1 - xi^2))*dWxi/p;
    
    dA = (CB*cot(A)/p^2)*dt + sqrt(2*CB)*dWxi/p;
    dB = sqrt(2*CB)*dWphi/(p*sin(A));
    
    dphi = sqrt(2*CB)*dWphi/(p*sqrt(1 - xi^2));
    if ~isfinite(dphi)
        dphi = 0;
    end
    
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