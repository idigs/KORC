function radiatedPower(numTracers)

close all

% Energy of electron, in eV.
Eo = 3E6;
c = 2.9979E8; % Speed of light, in m/s
qe = 1.602176E-19; % Electron charge, in Coulombs
me = 9.109382E-31; % Electron mass, in kg.
vo = sqrt( 1 - (me*c^2/(Eo*qe))^2 );

a = 0.5;
Ro = 1.6;

[xo,yo,zo] = torusmap(a,Ro,rand(3,numTracers));

% radius = rand(1,numTracers);
% angle = 2*pi*rand(1,numTracers);
% 
% min_angle = 1; % in degrees
% max_angle = 3; % in degrees
% % angleo = min_angle + (max_angle - min_angle)*rand(1,numTracers); % initial azimuthal angle
% angleo = 10; % initial azimuthal angle
% 
% xo = (Ro + a*sqrt(radius).*cos(angle)).*cos(angleo*pi/180);
% yo = (Ro + a*sqrt(radius).*cos(angle)).*sin(angleo*pi/180);
% zo = a*sqrt(radius).*sin(angle);


h = figure;
plot3(xo,yo,zo,'k.','MarkerSize',6)
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
savefig(h,'initial_distribution_tracers.fig')

tracers_RZ = cell(1,numTracers);
tracers_k = cell(1,numTracers);
% ko = zeros(1,numTracers);

parfor ii=1:numTracers
    try
        
    ST = particleOrbits('','','2D',[],[1E5,1E-2,10],[-1,1],[xo(ii),yo(ii),zo(ii)],[vo,30],false);
    
    close all
    
    zeta = atan2(ST.PP.X(:,2),ST.PP.X(:,1));
    zeta(zeta<0) = zeta(zeta<0) + 2*pi;
    locs = find(abs(diff(zeta)) > 6);
    
    R = sqrt(ST.PP.X(locs,1).^2 + ST.PP.X(locs,2).^2);
    Z = ST.PP.X(locs,3);
    
    tracers_RZ{ii} = [R, Z];
    tracers_k{ii} = ST.PP.k(locs);
%     ko(ii) = ST.PP.k(2);
    catch
        disp('Exception!')
    end

    disp(['Tracer No. ' num2str(ii)])
end

% Ro = sqrt(xo.^2 + yo.^2); % initial radial coordinate
% 
% S = 6*ones(numTracers,1);
% C = k.^2/max(abs(k.^2));
% 
% h = figure;
% subplot(2,1,1)
% scatter3(RZ(1,:),RZ(2,:),k.^2,S,C)
% box on
% axis on
% xlabel('R','Interpreter','latex','FontSize',16)
% ylabel('Z','Interpreter','latex','FontSize',16)
% zlabel('$\kappa^2(R,Z)$','Interpreter','latex','FontSize',16)
% 
% subplot(2,1,2)
% scatter3(Ro,zo,k-ko,S,C)
% box on
% axis on
% xlabel('R','Interpreter','latex','FontSize',16)
% ylabel('Z','Interpreter','latex','FontSize',16)
% zlabel('$\kappa - \kappa_o$','Interpreter','latex','FontSize',16)
% colormap(jet)
% colorbar

end