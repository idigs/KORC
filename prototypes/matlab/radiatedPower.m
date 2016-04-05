function radiatedPower(numTracers)

close all

a = 0.5;
Ro = 1.6;

radius = a*rand(1,numTracers);
angle = 2*pi*rand(1,numTracers);

min_angle = 1; % in degrees
max_angle = 3; % in degrees
angleo = min_angle + (max_angle - min_angle)*rand(1,numTracers); % initial azimuthal angle

xo = (Ro + sqrt(radius).*cos(angle)).*cos(angleo*pi/180);
yo = (Ro + sqrt(radius).*cos(angle)).*sin(angleo*pi/180);
zo = sqrt(radius).*sin(angle);


h = figure;
plot3(xo,yo,zo,'k.','MarkerSize',6)
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
savefig(h,'initial_distribution_tracers.fig')

RZ = zeros(2,numTracers);
k = zeros(1,numTracers);

for ii=1:numTracers
    ST = particleOrbits('','','2D',[],[1E3,1E-2,1],[-1,1],[xo(ii),yo(ii),zo(ii)],[0.99,1],false);
    
    zeta = atan2(ST.PP.X(:,2),ST.PP.X(:,1));
    zeta(zeta<0) = zeta(zeta<0) + 2*pi;
    locs = find(abs(diff(zeta)) > 6);
    
    try
    R = sqrt(ST.PP.X(locs,1).^2 + ST.PP.X(locs,2).^2);
    Z = ST.PP.X(locs,3);
    
    RZ(:,ii) = [R(1); Z(1)];
    k(ii) = ST.PP.k(locs(1));
    catch
        disp('Exception!')
    end
%     R
    disp(['Tracer No. ' num2str(ii)])
end

S = 6*ones(numTracers,1);
C = k.^2/max(abs(k.^2));

scatter3(RZ(1,:),RZ(2,:),k.^2,S,C)
xlabel('R','Interpreter','latex','FontSize',16)
ylabel('Z','Interpreter','latex','FontSize',16)
zlabel('$\kappa^2$','Interpreter','latex','FontSize',16)
colormap(jet)
colorbar

end