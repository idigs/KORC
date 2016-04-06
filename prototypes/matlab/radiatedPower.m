function radiatedPower(numTracers,poolsize)

close all

% Energy of electron, in eV.
Eo = 3E6;
c = 2.9979E8; % Speed of light, in m/s
qe = 1.602176E-19; % Electron charge, in Coulombs
me = 9.109382E-31; % Electron mass, in kg.
vo = sqrt( 1 - (me*c^2/(Eo*qe))^2 );

% Minor and major of analytical toroidal magnetic field
a = 0.5;
Ro = 1.6;

% Uniform distribution of points in the torus defined by the flux surfaces
% of the analytical magnetic field.
[xo,yo,zo] = torusmap(a,Ro,rand(3,numTracers));

h = figure;
plot3(xo,yo,zo,'k.','MarkerSize',6)
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
savefig(h,'initial_distribution_tracers.fig')

RZ = cell(1,numTracers);
k = cell(1,numTracers);
ko = zeros(1,numTracers);
pitch = cell(1,numTracers);

poolobj = gcp('nocreate');
if isempty(poolobj)
    poolobj = parpool(poolsize);
end

parfor ii=1:numTracers
    try
        
        ST = particleOrbits_ProductionRuns('','','2D',[],[1E6,1E-2,10],[-1,1],[xo(ii),yo(ii),zo(ii)],[vo,50],false);
        
        close all
        
        zeta = atan2(ST.PP.X(:,2),ST.PP.X(:,1));
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        locs = find(abs(diff(zeta)) > 6);
        
        R = sqrt(ST.PP.X(locs,1).^2 + ST.PP.X(locs,2).^2);
        Z = ST.PP.X(locs,3);
        
        pitch{ii} = atan2(ST.PP.vperp(locs),ST.PP.vpar(locs));
        
        RZ{ii} = [R, Z];
        k{ii} = ST.PP.k(locs);
        ko(ii) = ST.PP.k(1);
        
    catch
        disp('Exception!')
    end
    
    disp(['Tracer No. ' num2str(ii)])
    
    
end


Ro = sqrt(xo.^2 + yo.^2); % initial radial coordinate


for ii=2:numTracers
    k_max_previous = max(k{ii-1});
    k_max_current = max(k{ii});
    if k_max_previous < k_max_current
        k_max = k_max_current;
    else
        k_max = k_max_previous;
    end
end

h = figure;
subplot(2,1,1)
hold on
for ii=1:numTracers
    C = k{ii}.^2/k_max^2;
    scatter3(RZ{ii}(:,1),RZ{ii}(:,2),k{ii}.^2,6,C)
end
hold off
box on
axis on
xlabel('R','Interpreter','latex','FontSize',16)
ylabel('Z','Interpreter','latex','FontSize',16)
zlabel('$\kappa^2(R,Z)$','Interpreter','latex','FontSize',16)

subplot(2,1,2)
C = ko.^2/max(ko)^2;
scatter3(Ro,zo,ko,6,C)
box on
axis on
xlabel('R','Interpreter','latex','FontSize',16)
ylabel('Z','Interpreter','latex','FontSize',16)
zlabel('$\kappa_o$','Interpreter','latex','FontSize',16)
colormap(jet)
% colorbar

for ii=2:numTracers
    pitch_max_previous = max(pitch{ii-1});
    pitch_max_current = max(pitch{ii});
    if pitch_max_previous < pitch_max_current
        pitch_max = pitch_max_current;
    else
        pitch_max = pitch_max_previous;
    end
end

g = figure;
hold on
for ii=1:numTracers
    C = pitch{ii}/pitch_max;
    scatter3(RZ{ii}(:,1),RZ{ii}(:,2),pitch{ii},6,C)
end
hold off
box on
axis on
xlabel('R','Interpreter','latex','FontSize',16)
ylabel('Z','Interpreter','latex','FontSize',16)
zlabel('$\theta_{v_\perp/v_\parallel}$','Interpreter','latex','FontSize',16)
colormap(jet)

delete(poolobj);

end