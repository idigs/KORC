function radiatedPower(numTracers,poolsize)

close all

% Number of time iterations for calculating electrons' orbits
numTimeIt = 1E6;
cadence = 2E4;

% Initial pitch angle
pitcho = [10,20,50];

% Energy of electron, in eV.
Eo = [3E6];
c = 2.9979E8; % Speed of light, in m/s
qe = 1.602176E-19; % Electron charge, in Coulombs
me = 9.109382E-31; % Electron mass, in kg.

% Minor and major of analytical toroidal magnetic field
a = 0.5;
Ro = 1.6;

[xo,yo,zo] = torusmap(a,Ro,rand(3,numTracers));

h = figure;
plot3(xo,yo,zo,'k.','MarkerSize',6)
axis equal
xlabel('X','Interpreter','latex','FontSize',16)
ylabel('Y','Interpreter','latex','FontSize',16)
zlabel('Z','Interpreter','latex','FontSize',16)
savefig(h,'poloidal_plane_figures/initial_distribution_tracers.fig')
close(h)

for ee=1:numel(Eo)
    % Uniform distribution of points in the torus defined by the flux surfaces
    % of the analytical magnetic field.
    vo = sqrt( 1 - (me*c^2/(Eo(ee)*qe))^2 );
    
    for pp=1:numel(pitcho)
        
        RZ = cell(1,numTracers);
        k = cell(1,numTracers);
        T = cell(1,numTracers);
        ko = zeros(1,numTracers);
        pitch = cell(1,numTracers);
        
        EK = cell(1,numTracers);
        angularMomentum = cell(1,numTracers);
        mu = cell(1,numTracers);
        
        poolobj = gcp('nocreate');
        if isempty(poolobj)
            poolobj = parpool(poolsize);
        end
        
        local_pitcho = pitcho(pp);
        
        parfor pii=1:numTracers
            try
                ST = ...
                    particleOrbits_ProductionRuns('','','2D',[],[numTimeIt,1E-2,cadence],[-1,1],[xo(pii),yo(pii),zo(pii)],[vo,local_pitcho],false);
               
                zeta = atan2(ST.PP.X(:,2),ST.PP.X(:,1));
                zeta(zeta<0) = zeta(zeta<0) + 2*pi;
                locs = find(abs(diff(zeta)) > 6);
                
                R = sqrt(ST.PP.X(locs,1).^2 + ST.PP.X(locs,2).^2);
                Z = ST.PP.X(locs,3);
                
                pitch{pii} = atan2(ST.PP.vperp(locs),ST.PP.vpar(locs));
                
                RZ{pii} = [ST.PP.POINCARE.R, ST.PP.POINCARE.Z];
                k{pii} = ST.PP.POINCARE.k;
                T{pii} = ST.PP.POINCARE.T;
                ko(pii) = ST.PP.k(1);
                
                EK{pii} = ST.PP.EK;
                mu{pii} = ST.PP.mu);
                angularMomentum{pii} = ST.PP.angularMomentum;
                
            catch
                disp('Exception!')
            end
            disp(['Energy No.' num2str(ee) ' Pitch No.' num2str(pp) ' Tracer No. ' num2str(pii)])
        end
        
        filename = ['poloidal_plane_figures/var_Eo_' num2str(Eo(ee)) ...
            '_po_' num2str(pitcho(pp)) '.mat'];
        save(filename,'RZ','k','T','ko','EK','mu','angularMomentum','pitch')
        
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
            C = k{ii}/k_max;
            scatter3(RZ{ii}(:,1),RZ{ii}(:,2),k{ii},6,C)
        end
        hold off
        box on; axis on; grid on
        xlabel('R','Interpreter','latex','FontSize',16)
        ylabel('Z','Interpreter','latex','FontSize',16)
        zlabel('$\kappa(R,Z)$','Interpreter','latex','FontSize',16)
        
        subplot(2,1,2)
        C = ko/max(ko);
        scatter3(Ro,zo,ko,6,C)
        box on; axis on; grid on
        xlabel('R','Interpreter','latex','FontSize',16)
        ylabel('Z','Interpreter','latex','FontSize',16)
        zlabel('$\kappa_o$','Interpreter','latex','FontSize',16)
        colormap(jet)
        savefig(h,['poloidal_plane_figures/curvature_Eo_' num2str(Eo(ee)) '_po_' num2str(pitcho(pp)) '.fig'])
        close(h)
        
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
        box on; axis on; grid on
        xlabel('R','Interpreter','latex','FontSize',16)
        ylabel('Z','Interpreter','latex','FontSize',16)
        zlabel('$\theta_{v_\perp/v_\parallel}$ [rad]','Interpreter','latex','FontSize',16)
        colormap(jet)
        savefig(g,['poloidal_plane_figures/pitch_Eo_' num2str(Eo(ee)) '_po_' num2str(pitcho(pp)) '.fig'])
        close(g)
    end
end

delete(poolobj);

end