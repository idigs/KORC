function pitchAngleDependency(numAngles)

% Initial pitch angles used in the diagnostic, in degrees.
pitchAngle = linspace(0,80,numAngles);

% Energy of electron, in eV.
Eo = 3E6;
c = 2.9979E8; % Speed of light, in m/s
qe = 1.602176E-19; % Electron charge, in Coulombs
me = 9.109382E-31; % Electron mass, in kg.
vo = sqrt( 1 - (me*c^2/(Eo*qe))^2 );

cadence = 20;

h = figure;
hh = figure;
tvsR = figure;
strLegend = cell(1,numAngles);

xo = [1.35,1.85];
for jj=1:numel(xo)
    for ii=1:numAngles
        ST = particleOrbits('','','2D',[],[1E6,1E-2,cadence],[-1,1],[xo(jj),0,0],[vo,pitchAngle(ii)],false);
        
        zeta = atan2(ST.PP.X(:,2),ST.PP.X(:,1));
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        locs = find(abs(diff(zeta)) > 6);
        
        R = sqrt(ST.PP.X(:,1).^2 + ST.PP.X(:,2).^2);
        Z = ST.PP.X(:,3);
        
        time = ST.time/(2*pi/ST.params.wc);
        
        theta = atan2(ST.PP.vperp,ST.PP.vpar);
        
        figure(h)
        hold on
        plot(R(ST.params.inds),ST.PP.k(ST.params.inds),'.','MarkerSize',2)
        hold off
        
        figure(hh)
        hold on
        plot(theta(ST.params.inds),ST.PP.k(ST.params.inds),'.','MarkerSize',2)
        hold off
        
        figure(tvsR)
        hold on
        plot(R(ST.params.inds),theta(ST.params.inds),'.','MarkerSize',2)
        hold off
        
        
        g = figure;
        subplot(3,1,1)
        plot(time(ST.params.inds), ST.PP.k(ST.params.inds))
        box on
        grid on
        xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
        ylabel('Curvature $\kappa(t)$','Interpreter','latex','FontSize',16)
        
        subplot(3,1,2)
        plot(theta(ST.params.inds), ST.PP.k(ST.params.inds),'.','MarkerSize',2)
        box on
        grid on
        xlabel('$\theta_{v_\perp/v_\parallel}$ [rad]','Interpreter','latex','FontSize',16)
        ylabel('Curvature $\kappa(t)$','Interpreter','latex','FontSize',16)
        
        subplot(3,1,3)
        plot(time(ST.params.inds), ST.PP.T(ST.params.inds),time(ST.params.inds),0*time(ST.params.inds),'k--')
        box on
        grid on
        xlabel('Time $t$ [$\tau_e$]','Interpreter','latex','FontSize',16)
        ylabel('Torsion $\tau(t)$','Interpreter','latex','FontSize',16)
        savefig(g,['pitchAngleDependency/xo_' num2str(xo(jj)) '_' num2str(ii) '.fig'])
        close(g)
        
        strLegend{ii} = [num2str( pitchAngle(ii) ) '$^\circ$'];
        
        disp(['Tracer No. ' num2str(ii)])
    end
end

figure(h)
box on; grid on

ylabel('$\kappa(R)$','Interpreter','latex','FontSize',16)
set(legend(strLegend),'Interpreter','latex','Location','southwest','FontSize',14)
savefig(h,['pitchAngleDependency/xo_' num2str(xo(jj)) '_k_vs_R.fig'])

figure(hh)
box on; grid on
xlabel('$\theta_{v_\perp/v_\parallel}$ [rad]','Interpreter','latex','FontSize',16)
ylabel('$\kappa(\theta_{v_\perp/v_\parallel})$','Interpreter','latex','FontSize',16)
set(legend(strLegend),'Interpreter','latex','Location','southwest','FontSize',14)
savefig(hh,['pitchAngleDependency/xo_' num2str(xo(jj)) '_k_vs_theta.fig'])

figure(tvsR)
box on; grid on
xlabel('$R$ [m]','Interpreter','latex','FontSize',16)
ylabel('$\theta_{v_\perp/v_\parallel}$ [rad]','Interpreter','latex','FontSize',16)
set(legend(strLegend),'Interpreter','latex','Location','southwest','FontSize',14)
savefig(tvsR,['pitchAngleDependency/xo_' num2str(xo(jj)) '_theta_vs_R.fig'])

end