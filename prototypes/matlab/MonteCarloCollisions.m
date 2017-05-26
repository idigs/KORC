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
CO.VTe = sqrt(2*CO.Te*CO.params.qe/CO.params.me);

CO.np = np; % Number of particles

CO = initializeCollisionOperators(CO);

CO = normalize(CO);

% V = ThermalDistribution(CO);

u = sqrt(1 - (CO.params.me*CO.params.c^2./(ERE*CO.params.qe)).^2);

V = repmat([u;0;0],[1,CO.np]);

x = linspace(-1,1,500);
fx = exp(-x.^2/CO.VTe^2)/(CO.VTe*sqrt(pi));

snapshot = floor(numIt/5);
sh = figure;

for ii=1:CO.numIt
    for pp=1:CO.np
        V(:,pp) = collisionOperator(CO,V(:,pp)',CO.cop.dt);
    end
    if mod(ii,snapshot) == 0
        figure(sh);
        subplot(2,5,ii/snapshot)
        hold on;plot3(V(1,:),V(2,:),V(3,:),'r.');
        grid on;box on;axis equal;axis([-1,1,-1,1,-1,1]);hold off
        view([150,10])
        title(['$t=$' num2str(ii*CO.cop.dt*CO.norm.t) ' s'],'interpreter','latex')
        
        subplot(2,5,ii/snapshot+5)
        hold on
        histogram(V(1,:),25,'Normalization','pdf','LineStyle','none')
        histogram(V(2,:),25,'Normalization','pdf','LineStyle','none')
        histogram(V(3,:),25,'Normalization','pdf','LineStyle','none')
        plot(x,fx,'k','LineWidth',2)
        hold off
        legend({'$f(v_x)$','$f(v_y)$','$f(v_z)$','$f_M(v)$'},'Interpreter','latex')
        box on
        grid on
    end
end

CO.V = V;
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
CO.cop.VTe = CO.VTe;
CO.cop.delta = CO.cop.VTe/CO.params.c;
CO.cop.Gamma = CO.cop.ne*CO.params.qe^4*CO.cop.Clog/(4*pi*CO.params.ep^2);
CO.cop.Tauc = CO.params.me^2*CO.VTe^3/CO.cop.Gamma;

CO.cop.Tau = CO.params.me^2*CO.params.c^3/CO.cop.Gamma;
% CO.cop.Tau = CO.cop.Tauc;

CO.cop.Ec = CO.cop.ne*CO.params.qe^3*CO.cop.Clog/(4*pi*CO.params.ep^2*CO.params.me*CO.params.c^2);
CO.cop.ED = CO.cop.ne*CO.params.qe^3*CO.cop.Clog/(4*pi*CO.params.ep^2*CO.cop.Te);

Ef = CO.E;
CO.cop.Vc = CO.cop.VTe*sqrt(0.5*CO.cop.Ec/sqrt(Ef*Ef'));

CO.cop.c = CO.params.c;

g = linspace(1,5,1000); % relativistic gamma factor
E = g*CO.params.me*CO.params.c^2;%*linspace(6E5,50E6,200)*ST.params.qe;
Er = CO.params.me*CO.params.c^2; % Rest energy
u = CO.params.c*sqrt(1 - Er^2./E.^2);

disp(['Relativistic collisional time: ' num2str(CO.cop.Tau) ' s'])
disp(['Thermal collisional time: ' num2str(CO.cop.Tauc) ' s'])


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

CO.cop.psi = @(v) 0.5*( erf(CO.cop.x(v)) - 2*CO.cop.x(v).*exp(-CO.cop.x(v).^2)/sqrt(pi) )./CO.cop.x(v).^2;

CO.cop.CA = @(v) CO.cop.Gamma*CO.cop.psi(v)./v;

CO.cop.CF = @(v) CO.cop.Gamma*CO.cop.psi(v)/CO.cop.Te;

CO.cop.CB = @(v) (0.5*CO.cop.Gamma./v).*( CO.cop.Zeff + ...
    erf(CO.cop.x(v)) - CO.cop.psi(v) + 0.5*CO.cop.delta^4*CO.cop.x(v).^2 );

% CO.cop.CB = @(v) (0.5*CO.cop.Gamma./v).*( erf(CO.cop.x(v)) - CO.cop.psi(v) );

E = 1E-6*E/CO.params.qe;
xAxis = E;

figure;
subplot(3,1,1);
plot(xAxis,CO.cop.CA(u));
ylabel('$C_A$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,1,2);
plot(xAxis,CO.cop.CB(u));
ylabel('$C_B$ ($e B_0/m_e^3 c^2 $)','Interpreter','latex')
xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
subplot(3,1,3);
plot(xAxis,CO.cop.CF(u));
ylabel('$C_F$ ($e B_0/m_e^2 c $)','Interpreter','latex')
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

dW = random('norm',0,sqrt(dt),[3,1]);

CA = CO.cop.CA(v);
CB = CO.cop.CB(v);
CF = CO.cop.CF(v);

dU = [-2.0*CF*dt + sqrt(2*CA)*dW(1);
    sqrt(2*CB)*dW(2);
    sqrt(2*CB)*dW(3)];

U(1) = (U1+dU(1))*(v1*CO.b1') + (U2+dU(2))*(v2*CO.b1') + (U3+dU(3))*(v3*CO.b1');
U(2) = (U1+dU(1))*(v1*CO.b2') + (U2+dU(2))*(v2*CO.b2') + (U3+dU(3))*(v3*CO.b2');
U(3) = (U1+dU(1))*(v1*CO.b3') + (U2+dU(2))*(v2*CO.b3') + (U3+dU(3))*(v3*CO.b3');

U = U/sqrt( 1 + U*U' );

end
