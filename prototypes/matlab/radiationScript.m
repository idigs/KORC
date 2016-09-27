clear all
close all

% parameters
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
Bo = 3.285; % Magnetic field in T

lambda_min = 450E-9; % in meters
lambda_max = 950E-9; % in meters
% lambda_min = 1E-9; % in meters
% lambda_max = 1000E-9; % in meters
Nlambda = 100;
lambda = linspace(lambda_min,lambda_max,Nlambda);

NE = 100;
% Npitch = 9;

upper_integration_limit = 100;

Emin = 1E6;
Emax = 80E6;
E = linspace(Emin,Emax,NE); % Energy in eV
E = E*qe; % Energy in J
Eo = me*c^2; % rest energy in J
gammap = E/Eo;

p = sqrt(E.^2 - (Eo)^2)/c;

% pitch = (pi/180)*linspace(1,90,Npitch); % Pitch angle in radians
pitch = (pi/180)*[5,10,20,30,40,50,60,70,80];
Npitch = numel(pitch);

% Emission angle
Npsi = 10;
psi = (pi/180)*linspace(0,5,Npsi);


% Curvature of the electron orbit
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
k = 1E-2*k;

lambdac = zeros(NE,Npitch);
wc = zeros(NE,Npitch);
Psyn = zeros(Nlambda,NE,Npitch);
Psyn_tot = zeros(NE,Npitch);

for ii=1:NE
    lambdac(ii,:) = (4*pi/3)*R(ii,:).*(Eo./E(ii)).^3;
    wc(ii,:) = 2*pi*c./lambdac(ii,:);
end

% Spectral distribution of Psyn
C0 = 4*pi*c*qe^2/sqrt(3);
fun = @(x) besselk(5/3,x);
for ii=1:NE
    for pp=1:Npitch
        for ll=1:Nlambda
            lower_integration_limit = lambdac(ii,pp)/lambda(ll);
            if (lambda(ll) < lambdac(ii,pp)) && (lower_integration_limit < upper_integration_limit)
                Q = integral(fun,lower_integration_limit,upper_integration_limit);
                if (Q < 0)
                    error('Negative value at integral K(5/3,x)')
                end
                C1 = (Eo/E(ii))^2/lambda(ll)^3;
                Psyn(ll,ii,pp) =  C0*C1*Q;  
            end
            
            Psyn_tot(ii,pp) = trapz(lambda,squeeze(Psyn(:,ii,pp)));
        end
    end    
end

% Angular distribution of Psyn
Psyn_psi = zeros(NE,Npitch,Npsi);
Psyn_psi_tot = zeros(NE,Npitch);
for ii=1:NE
    for pp=1:Npitch
        x = (gammap(ii)*psi).^2;
        A0 = c*qe^2*k(ii,pp)^2*gammap(ii)^5;
                
        Psyn_psi(ii,pp,:) = A0*(1 + x).^(-5/2).*(7/16 + ...
            (5/16)*x./(1+x));
        
        Psyn_psi_tot(ii,pp) = trapz(psi,squeeze(Psyn_psi(ii,pp,:)));
    end    
end

Psyn_psi_lambda = zeros(Nlambda,NE,Npitch,Npsi);
Psyn_psi_lambda_tot = zeros(NE,Npitch);
% Angular and spectral distribution of Psyn
for ii=1:NE
    for pp=1:Npitch
        tmp = zeros(1,Nlambda);
        for ll=1:Nlambda
            if (lambda(ll) < lambdac(ii,pp))
                x = (gammap(ii)*psi).^2;
                lambda_ratio = lambdac(ii,pp)/lambda(ll);
                C0 = 3*c*qe^2*k(ii,pp)/(2*pi*lambda(ll)^2);
                zeta = 0.5*lambda_ratio*(1 + x).^(3/2);
                
                Psyn_psi_lambda(ll,ii,pp,:) = ...
                    C0*lambda_ratio^2*gammap(ii)^2*(1 + x).^2.*( besselk(2/3,zeta).^2 + ...
                    (x./(1 + x)).*besselk(1/3,zeta).^2 );
            end
            
            tmp(ll) = trapz(psi,squeeze(Psyn_psi_lambda(ll,ii,pp,:)));
        end
        Psyn_psi_lambda_tot(ii,pp) = trapz(lambda,tmp);
    end
end


qe = qe/3E9;
E = E/1E7;
E = E/qe;

lambdac = lambdac/1E2;
lambda = lambda/1E2;
Psyn = 1E-7*Psyn;
Psyn_psi = 1E-7*Psyn_psi;
Psyn_psi_lambda = 1E-7*Psyn_psi_lambda;

Psyn_psi_lambda_tot = 1E-7*Psyn_psi_lambda_tot;
Psyn_psi_tot = 1E-7*Psyn_psi_tot;
Psyn_tot = 1E-7*Psyn_tot;
% Beyond this point all variables are in SI units

E = E/1E6; % in MeV
lambdac = lambdac/1E-9; % in nm
pitch = 180*pitch/pi;

% A = log10(lambdac);
% 
% figure
% surf(pitch,E,A,'LineStyle','none')
% axis([min(pitch) max(pitch) min(E) max(E) min(min(A)) max(max(A))])
% box on
% colormap(jet(512))
% xlabel('Pitch angle $\theta$ ($^\circ$)','Interpreter','latex')
% ylabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')

lambda = lambda/1E-9;
psi = 180*psi/pi;

figure
subplot(1,3,1)
surf(pitch,E,Psyn_tot,'LineStyle','none')
axis([min(pitch) max(pitch) min(E) max(E) min(min(Psyn_tot)) max(max(Psyn_tot))])
box on; view([0,90])
cmp = colormap(jet(1024));
colorbar
xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
ylabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
title('$\int P_{syn}(\mathcal{E},\theta,\lambda) d\lambda$','Interpreter','latex')
subplot(1,3,2)
surf(pitch,E,Psyn_psi_tot,'LineStyle','none')
axis([min(pitch) max(pitch) min(E) max(E) min(min(Psyn_psi_tot)) max(max(Psyn_psi_tot))])
box on; view([0,90])
cmp = colormap(jet(1024));
colorbar
xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
ylabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
title('$\int P_{syn}(\mathcal{E},\theta,\psi) d\psi$','Interpreter','latex')
subplot(1,3,3)
surf(pitch,E,Psyn_psi_lambda_tot,'LineStyle','none')
axis([min(pitch) max(pitch) min(E) max(E) min(min(Psyn_psi_lambda_tot)) max(max(Psyn_psi_lambda_tot))])
box on; view([0,90])
cmp = colormap(jet(1024));
colorbar
xlabel('$\theta$ ($^\circ$)','Interpreter','latex')
ylabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
title('$\int \int P_{syn}(\mathcal{E},\theta,\lambda,\psi) d\psi d\lambda$','Interpreter','latex')
%%
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

% figure % Psyn_psi = zeros(NE,Npitch,Npsi);
% for ii=1:Npitch
%     subplot(3,3,ii)
%     A = log10(squeeze(Psyn_psi(:,ii,:)));
%     surf(psi,E,A,'LineStyle','none')
%     axis([min(psi) max(psi) min(E) max(E) min(min(A)) max(max(A))])
%     box on; view([0,90])
%     cmp = colormap(jet(1024));
% %     new_colormap = cmp;
% %     new_colormap(1,:) = [1,1,1];
% %     colormap(new_colormap);
%     colorbar
%     shading interp
%     title(['$\theta=$' num2str(pitch(ii)) '$^\circ$'],'Interpreter','latex')
%     xlabel('$\psi$ ($^\circ$)','Interpreter','latex')
%     ylabel('Energy $\mathcal{E}$ (MeV)','Interpreter','latex')
% end
%%

for jj=Npsi:-1:1
    figure('name',['psi = ' num2str(psi(jj)) ' (degrees)'])
    % figure
    for ii=1:Npitch
        subplot(3,3,ii)
        A = squeeze(Psyn_psi_lambda(:,:,ii,jj));
%         A = zeros(Nlambda,NE);
%         for aa=1:NE
%             for bb=1:Nlambda
%                 A(bb,aa) = trapz(psi,squeeze(Psyn_psi_lambda(bb,aa,ii,:)));
%             end
%         end
        surf(E,lambda,A,'LineStyle','none')
        try
            axis([min(E) max(E) min(lambda) max(lambda) min(min(A)) max(max(A))])
        catch
            axis([min(E) max(E) min(lambda) max(lambda)])
        end
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
end