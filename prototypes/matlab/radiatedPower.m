function radiatedPower(numTracers)

% close all

a = 0.5;
Ro = 1.6;

radius = rand(1,numTracers);
angle = 2*pi*rand(1,numTracers);

min_angle = 1; % in degrees
max_angle = 3; % in degrees
% angleo = min_angle + (max_angle - min_angle)*rand(1,numTracers); % initial azimuthal angle
angleo = 10; % initial azimuthal angle

xo = (Ro + a*sqrt(radius).*cos(angle)).*cos(angleo*pi/180);
yo = (Ro + a*sqrt(radius).*cos(angle)).*sin(angleo*pi/180);
zo = a*sqrt(radius).*sin(angle);


h = figure;
plot3(xo,yo,zo,'k.','MarkerSize',6)
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
savefig(h,'initial_distribution_tracers.fig')

RZ = zeros(2,numTracers);
k = zeros(1,numTracers);
ko = zeros(1,numTracers);

for ii=1:numTracers
    try
        
    ST = particleOrbits('','','2D',[],[1E4,1E-2,1],[-1,1],[xo(ii),yo(ii),zo(ii)],[0.985387037755342,5],true);
    
    close all
    
    zeta = atan2(ST.PP.X(:,2),ST.PP.X(:,1));
    zeta(zeta<0) = zeta(zeta<0) + 2*pi;
    locs = find(abs(diff(zeta)) > 6);
    
    R = sqrt(ST.PP.X(locs,1).^2 + ST.PP.X(locs,2).^2);
    Z = ST.PP.X(locs,3);
    
    RZ(:,ii) = [R(1); Z(1)];
    k(ii) = ST.PP.k(locs(1));
    ko(ii) = ST.PP.k(2);
    catch
        disp('Exception!')
    end

    disp(['Tracer No. ' num2str(ii)])
end

Ro = sqrt(xo.^2 + yo.^2); % initial radial coordinate

S = 6*ones(numTracers,1);
C = k.^2/max(abs(k.^2));

h = figure;
subplot(2,1,1)
scatter3(RZ(1,:),RZ(2,:),k.^2,S,C)
box on
axis on
xlabel('R','Interpreter','latex','FontSize',16)
ylabel('Z','Interpreter','latex','FontSize',16)
zlabel('$\kappa^2(R,Z)$','Interpreter','latex','FontSize',16)

subplot(2,1,2)
scatter3(Ro,zo,k-ko,S,C)
box on
axis on
xlabel('R','Interpreter','latex','FontSize',16)
ylabel('Z','Interpreter','latex','FontSize',16)
zlabel('$\kappa - \kappa_o$','Interpreter','latex','FontSize',16)
colormap(jet)
colorbar

end