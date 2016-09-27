clear all
close all

% parameters
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
Bo = 3.285; % Magnetic field in T

lambda_min = 450E-9; % in meters
lambda_max = 950E-9; % in meters
Nlambda = 200;
lambda = linspace(lambda_min,lambda_max,Nlambda);

NE = 200;
% Npitch = 9;

upper_integration_limit = 100;

Emin = 1E6;
Emax = 80E6;
E = linspace(Emin,Emax,NE); % Energy in eV
E = E*qe; % Energy in J
Eo = me*c^2; % rest energy in J

p = sqrt(E.^2 - (Eo)^2)/c;

% pitch = (pi/180)*linspace(1,90,Npitch); % Pitch angle in radians
pitch = (pi/180)*[5,10,20,30,40,50,60,70,80];
Npitch = numel(pitch);

k = zeros(NE,Npitch);
for ii=1:NE
    k(ii,:) = qe*Bo*sin(pitch)./p(ii);
end
R = 1./k;

% Below this point all variables are in c.g.s. units
Bo = 1E4*Bo;
c = 1E2*c;
qe = 3E9*abs(qe);
m = 1E3*me;

lambda = 1E2*lambda; % in cm

E = 1E7*E; % Energy in erg
Eo = 1E7*Eo; % Energy in erg

R = 1E2*R;

lambdac = zeros(NE,Npitch);
Psyn = zeros(Nlambda,NE,Npitch);

for ii=1:NE
    lambdac(ii,:) = (4*pi/3)*R(ii,:).*(Eo./E(ii)).^3;
end

C0 = 4*pi*c*qe^2/sqrt(3);
fun = @(x) besselk(5/3,x);
for ii=1:NE
    for pp=1:Npitch
        for ll=1:Nlambda
            if (lambda(ll) > lambdac(ii,pp))
                Psyn(ll,ii,pp) = 0;
            else
                lower_integration_limit = lambdac(ii,pp)/lambda(ll);
                
                if (lower_integration_limit > upper_integration_limit)
                    Psyn(ll,ii,pp) = 0;
                else
                    Q = integral(fun,lower_integration_limit,upper_integration_limit);
                    if (Q < 0)
                        error('Negative value at integral K(5/3,x)')
                    end
                    C1 = (Eo/E(ii))^2/lambda(ll)^3;
                    Psyn(ll,ii,pp) =  C0*C1*Q;
                end
            end
        end
    end
    
end

qe = qe/3E9;
E = E/1E7;
E = E/qe;

lambdac = lambdac/1E2;
lambda = lambda/1E2;
Psyn = 1E-7*Psyn;
% Beyond this point all variables are in SI units

E = E/1E6; % in MeV
lambdac = lambdac/1E-9; % in nm
pitch = 180*pitch/pi;

A = log10(lambdac);

figure
surf(pitch,E,A,'LineStyle','none')
axis([min(pitch) max(pitch) min(E) max(E) min(min(A)) max(max(A))])
box on
colormap(jet(512))
xlabel('Pitch angle $\theta$ ($^\circ$)','Interpreter','latex')
ylabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')

lambda = lambda/1E-9;

figure
for ii=1:Npitch
    subplot(3,3,ii)
    A = squeeze(Psyn(:,:,ii));
    surf(E,lambda,A,'LineStyle','none')
    axis([min(E) max(E) min(lambda) max(lambda) min(min(A)) max(max(A))])
    box on; view([0,90])
    cmp = colormap(jet(1024));
    new_colormap = cmp;
    new_colormap(1,:) = [1,1,1];
    colormap(new_colormap);
    colorbar
    shading interp
    title(['$\theta=$' num2str(pitch(ii)) '$^\circ$'],'Interpreter','latex')
    ylabel('$\lambda$ (nm)','Interpreter','latex')
    xlabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
end


