function ST = syntheticCameraFortran(path,lambdas,visible,range)
% ST = syntheticCameraFortran('/media/l8c/FantomHD/SimulationOutputs/Avalanche/Z1/',[400E-9,900E-9],'on',[99,100])
close all

ST = struct;
ST.path = path;
ST.visible = visible;
ST.lambdas = lambdas;

ST.range = range;
ST.num_snapshots = ST.range(2) - ST.range(1) + 1;

ST.params = loadSimulationParameters(ST);

ST.time = ...
    ST.params.simulation.dt*double(ST.params.simulation.output_cadence)*double(ST.range(1):1:ST.range(2));

ST.data = loadData(ST);

% plotSyntheticCameraAnalysis(ST)

generateFigures(ST)

end

function params = loadSimulationParameters(ST)
params = struct;

info = h5info([ST.path 'simulation_parameters.h5']);

for ii=1:length(info.Groups)
    for jj=1:length(info.Groups(ii).Datasets)
        name = info.Groups(ii).Name(2:end);
        subname = info.Groups(ii).Datasets(jj).Name;
        params.(name).(subname) = ...
            h5read(info.Filename,['/' name '/' subname]);
    end
end

info = h5info([ST.path 'synthetic_camera.h5']);

for ii=1:length(info.Groups)
    for jj=1:length(info.Groups(ii).Datasets)
        name = info.Groups(ii).Name(2:end);
        subname = info.Groups(ii).Datasets(jj).Name;
        params.(name).(subname) = ...
            h5read(info.Filename,['/' name '/' subname]);
    end
end

try
    info = h5info([ST.path 'avalanche_parameters.h5']);
    
    for ii=1:length(info.Groups)
        for jj=1:length(info.Groups(ii).Datasets)
            name = info.Groups(ii).Name(2:end);
            subname = info.Groups(ii).Datasets(jj).Name;
            params.(name).(subname) = ...
                h5read(info.Filename,['/' name '/' subname]);
        end
    end
catch
end
end

function data = loadData(ST)
data = struct;

it = ST.range(1):1:ST.range(2);

NX = ST.params.synthetic_camera_params.num_pixels(1);
NY = ST.params.synthetic_camera_params.num_pixels(2);
Nl = ST.params.synthetic_camera_params.Nlambda;

filename = [ST.path 'synthetic_camera_snapshots.h5'];

if (ST.params.synthetic_camera_params.integrated_opt == 0)
    list = {'Psyn_angular_pixel','np_angular_pixel','Psyn_lambda_pixel','np_lambda_pixel','Psyn_pplane','np_pplane'};
    
    for ll=1:length(list)
        disp(['Loading ' list{ll}])
        for ss=1:ST.params.simulation.num_species
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,Nl,ST.num_snapshots);
            for ii=1:numel(it)
                dataset = ...
                    ['/' num2str(it(ii)*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                data.(['sp' num2str(ss)]).(list{ll})(:,:,:,ii) = h5read(filename, dataset);
            end
        end
    end
end


if (ST.params.synthetic_camera_params.integrated_opt == 1)
    list={'Psyn_angular_pixel','np_angular_pixel','Psyn_lambda_pixel','np_lambda_pixel','Psyn_pplane','np_pplane',...
        'PTot_pplane','Psyn_pplane','np_pplane','P_lambda'};
else
    list={'PTot_pplane'};
end

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    for ss=1:ST.params.simulation.num_species
        if strcmp(list{ll},'P_lambda')
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(Nl,ST.num_snapshots);
            for ii=1:numel(it) % Here
                dataset = ...
                    ['/' num2str(it(ii)*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                data.(['sp' num2str(ss)]).(list{ll})(:,ii) = h5read(filename, dataset);
            end
        else
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,ST.num_snapshots);
            for ii=1:numel(it) % Here
                dataset = ...
                    ['/' num2str(it(ii)*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                data.(['sp' num2str(ss)]).(list{ll})(:,:,ii) = h5read(filename, dataset);
            end
        end
    end
end


end

function plotSyntheticCameraAnalysis(ST)
disp('Plotting snapshots...')

lambda = ST.params.synthetic_camera_params.lambda;
xAxis = ST.params.synthetic_camera_params.pixels_nodes_x;
yAxis = ST.params.synthetic_camera_params.pixels_nodes_y;

NX = ST.params.synthetic_camera_params.num_pixels(1);
NY = ST.params.synthetic_camera_params.num_pixels(2);

RAxis = ST.params.poloidal_plane_params.nodes_R;
ZAxis = ST.params.poloidal_plane_params.nodes_Z;
NR = ST.params.poloidal_plane_params.grid_dims(1);
NZ = ST.params.poloidal_plane_params.grid_dims(1);

[~,i1] = min(abs(lambda - ST.lambdas(1)));
[~,i2] = min(abs(lambda - ST.lambdas(2)));
Nl = i2 - i1 + 1;

for ss=1:ST.params.simulation.num_species
    disp(['Species: ' num2str(ss)])
    Psyn_lambda_pixel = zeros(NX,NY,Nl);
    np_lambda_pixel = zeros(NX,NY,Nl);
    Psyn_angular_pixel = zeros(NX,NY);
    np_angular_pixel = zeros(NX,NY);
    
    Psyn_lambda_pplane = zeros(NX,NY,Nl);
    np_lambda_pplane = zeros(NX,NY,Nl);
    Psyn_pplane = zeros(NX,NY);
    np_pplane = zeros(NX,NY);
    
    for ii=1:NX
        for jj=1:NY
            Psyn_lambda_pixel(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel(ii,jj,i1:i2,:),4);
            np_lambda_pixel(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel(ii,jj,i1:i2,:),4);
            
            Psyn_angular_pixel(ii,jj) = trapz(lambda(i1:i2),Psyn_lambda_pixel(ii,jj,:));
            np_angular_pixel(ii,jj) = sum(np_lambda_pixel(ii,jj,:),3);
        end
    end
    
    for ii=1:NR
        for jj=1:NZ
            Psyn_lambda_pplane(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_pplane(ii,jj,i1:i2,:),4);
            np_lambda_pplane(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_pplane(ii,jj,i1:i2,:),4);
            
            Psyn_pplane(ii,jj) = trapz(lambda(i1:i2),Psyn_lambda_pplane(ii,jj,:));
            np_pplane(ii,jj) = sum(np_lambda_pplane(ii,jj,:),3);
        end
    end
    
    % Convert from m to nm
    axis_lambda = 1E9*lambda(i1:i2);
    Psyn_lambda_pixel = 1E-9*Psyn_lambda_pixel;
    
    h = figure;
    subplot(4,2,[1 3])
    surfc(xAxis,yAxis,Psyn_angular_pixel','LineStyle','none')
    colormap(jet); hc = colorbar('Location','eastoutside');
    xlabel(hc,'$P_{syn}$ (Photon/s)','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
    
    figure(h);
    subplot(4,2,5)
%     hold on
%     for ii=1:NX
%         for jj=1:NY
%             plot(axis_lambda,squeeze(Psyn_lambda_pixel(ii,jj,:)))
%         end
%     end
%     hold off
    box on;
    ylabel('$P_{syn}$ (Photon/s/nm)','FontSize',12,'Interpreter','latex')
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
      
    figure(h);
    subplot(4,2,[2 4])
    surfc(xAxis,yAxis,np_angular_pixel','LineStyle','none')
    colormap(jet);  hc = colorbar('Location','eastoutside');
    xlabel(hc,'Number of RE','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$y$-axis','FontSize',14,'Interpreter','latex')
    xlabel('$x$-axis','FontSize',14,'Interpreter','latex')
    
    figure(h);
    subplot(4,2,6)
%     hold on
%     for ii=1:NX
%         for jj=1:NY
%             plot(axis_lambda,squeeze(np_lambda_pixel(ii,jj,:)))
%         end
%     end
%     hold off
    box on;
    ylabel('Number of RE','FontSize',12,'Interpreter','latex')
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')    
    
    figure(h);
    subplot(4,2,[7 8])
    yyaxis left 
    set(gca,'YColor',[0,0,1])
    fy = squeeze(sum(sum(Psyn_lambda_pixel,1),2));
%     fy = fy/max(fy);
    plot(axis_lambda,fy,'b','LineWidth',2)
    ylabel('$P_{syn}$ (Photon/s/nm)','FontSize',12,'Interpreter','latex')
    ylim([0 max(fy)])
    yyaxis right
    set(gca,'YColor',[1,0,0])
    fy = squeeze(sum(sum(np_lambda_pixel,1),2));
    plot(axis_lambda,fy,'r','LineWidth',2)
    ylim([0 max(fy)])
    ylabel('Number of RE','FontSize',12,'Interpreter','latex')
    box on;
    xlim([min(axis_lambda) max(axis_lambda)])
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
    
%     saveas(h,[ST.path 'SyntheticCameraFortran_pixel_ss_' num2str(ss)],'fig')
    
    
    g = figure;
    subplot(4,2,[1 3])
    surfc(RAxis,ZAxis,Psyn_pplane','LineStyle','none')
    colormap(jet); hc = colorbar('Location','eastoutside');
    xlabel(hc,'$P_{syn}$ (Photon/s)','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$Z$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$R$-axis','FontSize',12,'Interpreter','latex')
    
    figure(g);
    subplot(4,2,5)
    hold on
    for ii=1:NR
        for jj=1:NZ
            plot(axis_lambda,squeeze(Psyn_lambda_pplane(ii,jj,:)))
        end
    end
    hold off
    box on;
    ylabel('$P_{syn}$ (Photon/s/nm)','FontSize',12,'Interpreter','latex')
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
      
    figure(g);
    subplot(4,2,[2 4])
    surfc(RAxis,ZAxis,np_pplane','LineStyle','none')
    colormap(jet);  hc = colorbar('Location','eastoutside');
    xlabel(hc,'Number of RE','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$Z$-axis','FontSize',14,'Interpreter','latex')
    xlabel('$R$-axis','FontSize',14,'Interpreter','latex')
    
    figure(g);
    subplot(4,2,6)
    hold on
    for ii=1:NR
        for jj=1:NZ
            plot(axis_lambda,squeeze(np_lambda_pplane(ii,jj,:)))
        end
    end
    hold off
    box on;
    ylabel('Number of RE','FontSize',12,'Interpreter','latex')
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')    
    
    figure(g);
    subplot(4,2,[7 8])
    yyaxis left 
    set(gca,'YColor',[0,0,1])
    fy = squeeze(sum(sum(Psyn_lambda_pplane,1),2));
    fy = fy/max(fy);
    plot(axis_lambda,fy,'b','LineWidth',2)
    ylabel('$P_{syn}$ (Photon/s/nm)','FontSize',12,'Interpreter','latex')
    ylim([0 max(fy)])
    yyaxis right
    set(gca,'YColor',[1,0,0])
    fy = squeeze(sum(sum(np_lambda_pplane,1),2));
    plot(axis_lambda,fy,'r','LineWidth',2)
    ylim([0 max(fy)])
    ylabel('Number of RE','FontSize',12,'Interpreter','latex')
    box on;
    xlim([min(axis_lambda) max(axis_lambda)])
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
    
%     saveas(g,[ST.path 'SyntheticCameraFortran_pplane_ss_' num2str(ss)],'fig')
end
end

function Psyn = singleParticleSpectrum(ST,lambda,g,eta)
Psyn = zeros(size(lambda));

q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;
v = c*sqrt(1-1/g^2);
ep = 8.854E-12;% Electric permitivity

k = q*ST.params.fields.Bo*sin(eta)/(g*m*v);
l = lambda; 
lc = 4*pi/(3*k*g^3);

z = lc./l;

BK53 = @(x) besselk(5/3,x);
IntBKv = @(nu,x) (pi/sqrt(2))*(1 - 0.25*(4*nu^2 -1))*erfc(sqrt(x)) + ...
    0.25*(4*nu^2 - 1)*sqrt(0.5*pi./x).*exp(-x);

for ii=1:numel(z)
    if (z(ii) < 0.5)
        a = (2.16/2^(2/3))*z(ii)^(1/3);
        Psyn(ii) = integral(BK53,z(ii),a) + IntBKv(5/3,a);
    elseif (z(ii) >= 0.5) && (z(ii) < 2.5)
        a = 0.72*(z(ii) + 1);
        Psyn(ii) = integral(BK53,z(ii),a) + IntBKv(5/3,a);
    else
        Psyn(ii) = IntBKv(5/3,z(ii));
    end
end

Psyn = c*q^2*Psyn./(sqrt(3)*ep*g^2*l.^3);
end

function Psyn = averagedSpectrum(ST)
Np = 50;
Nchi = 100;

q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;

Ebar = ST.params.avalanche_pdf_params.Epar/ST.params.avalanche_pdf_params.Ec;
Zeff = ST.params.avalanche_pdf_params.Zeff;
Ehat = (Ebar - 1)/(1 + Zeff);
pmax = ST.params.avalanche_pdf_params.max_p;
pmin = ST.params.avalanche_pdf_params.min_p;
pitchmax = deg2rad(ST.params.avalanche_pdf_params.max_pitch_angle);
chimin = cos(pitchmax);
Cz = sqrt(3*(Zeff + 5)/pi)*ST.params.avalanche_pdf_params.Clog;

g = @(p) sqrt(p.^2 + 1);
eta = @(x) acos(x);

f = @(p,x) (Ehat/Cz)*p.*exp( -p.*(x/Cz + 0.5*Ehat*(1 - x.^2)./x) )./x;

% sanityIntegral = integral2(f,pmin,pmax,chimin,1);
% sanityIntegral = integral2(f,0,500,0,1);

p = linspace(pmin,pmax,Np);
pitch = linspace(0,pitchmax,Nchi);
chi = cos(pitch);
% chi = linspace(chimin,1,Nchi);
% pitch = eta(chi);

l = ST.params.synthetic_camera_params.lambda;
Psyn = zeros(size(l));
Psyn_p_chi = zeros(numel(l),Np,Nchi);
Psyn_p = zeros(numel(l),Np);

for ll=1:numel(l)
    for pp=1:Np
        for cc=1:Nchi
            Psyn_p_chi(ll,pp,cc) = f(p(pp),chi(cc))*singleParticleSpectrum(ST,l(ll),g(p(pp)),eta(chi(cc)));
        end
        Psyn_p(ll,pp) = trapz(fliplr(chi),squeeze(Psyn_p_chi(ll,pp,:)));
    end
    Psyn(ll) = trapz(p,squeeze(Psyn_p(ll,:)));
end

E = (g(p)*m*c^2/q)/1E6; % MeV
xAxis = rad2deg(pitch);
lAxis = l/1E-9;

% figure;
% surf(E,lAxis,squeeze(Psyn_p),'LineStyle','none');
% colormap(jet);
% ylabel('$\lambda$ (nm)','Interpreter','latex')
% xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex')


% figure
% plot(lAxis,Psyn)
% xlabel('$\lambda$ (nm)','Interpreter','latex')
% ylabel('$P_{R}$ (Watts)','Interpreter','latex')
end

function generateFigures(ST)
% Psyn_L1 = PTot
% Psyn_L2 = Psyn(lambda) in poloidal plane
% Psyn_L3 = Psyn(lambda) in pixel plane
% Psyn_L4 = Psyn(lambda,psi,chi) in pixel plane
disp('Plotting snapshots...')

lambda = ST.params.synthetic_camera_params.lambda;
xAxis = ST.params.synthetic_camera_params.pixels_nodes_x;
yAxis = ST.params.synthetic_camera_params.pixels_nodes_y;

NX = ST.params.synthetic_camera_params.num_pixels(1);
NY = ST.params.synthetic_camera_params.num_pixels(2);

RAxis = ST.params.poloidal_plane_params.nodes_R;
ZAxis = ST.params.poloidal_plane_params.nodes_Z;
NR = ST.params.poloidal_plane_params.grid_dims(1);
NZ = ST.params.poloidal_plane_params.grid_dims(1);

[~,i1] = min(abs(lambda - ST.lambdas(1)));
[~,i2] = min(abs(lambda - ST.lambdas(2)));
Nl = i2 - i1 + 1;

for ss=1:ST.params.simulation.num_species
    
    num_particles = double(ST.params.simulation.nmpi)*double(ST.params.species.ppp);
    
    disp(['Species: ' num2str(ss)])
    
    if (ST.params.synthetic_camera_params.integrated_opt == 0)
        Psyn_L4_lambda = zeros(NX,NY,Nl);
        np_L4_lambda = zeros(NX,NY,Nl);
        
        Psyn_L3_lambda = zeros(NX,NY,Nl);
        np_L3_lambda = zeros(NX,NY,Nl);
        
        Psyn_L2_lambda = zeros(NR,NZ,Nl);
        np_L2_lambda = zeros(NR,NZ,Nl);
    end
    
    Psyn_L4 = zeros(NX,NY);
    np_L4 = zeros(NX,NY);
    
    Psyn_L3 = zeros(NX,NY);
    np_L3 = zeros(NX,NY);
    
    Psyn_L2 = zeros(NR,NZ);
    np_L2 = zeros(NR,NZ);
    
    Psyn_L1 = zeros(NR,NZ);
    
    if (ST.params.synthetic_camera_params.integrated_opt == 1)
        np_L3 = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel,3);
        np_L4 = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel,3);
        
        Psyn_L3 = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel,3);
        Psyn_L4 = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel,3);
    else
        for ii=1:NX
            for jj=1:NY
                if numel(size(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel)) == 4
                    Psyn_L3_lambda(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel(ii,jj,i1:i2,:),4);
                    Psyn_L4_lambda(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel(ii,jj,i1:i2,:),4);
                    
                    np_L3_lambda(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel(ii,jj,i1:i2,:),4);
                    np_L4_lambda(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel(ii,jj,i1:i2,:),4);
                else
                    Psyn_L3_lambda(ii,jj,:) = ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel(ii,jj,i1:i2);
                    Psyn_L4_lambda(ii,jj,:) = ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel(ii,jj,i1:i2);
                    
                    np_L3_lambda(ii,jj,:) = ST.data.(['sp' num2str(ss)]).np_lambda_pixel(ii,jj,i1:i2);
                    np_L4_lambda(ii,jj,:) = ST.data.(['sp' num2str(ss)]).np_angular_pixel(ii,jj,i1:i2);
                end
                
                np_L3(ii,jj) = sum(np_L3_lambda(ii,jj,:),3);
                np_L4(ii,jj) = sum(np_L4_lambda(ii,jj,:),3);
                
                Psyn_L3(ii,jj) = trapz(lambda(i1:i2),Psyn_L3_lambda(ii,jj,:));
                Psyn_L4(ii,jj) = trapz(lambda(i1:i2),Psyn_L4_lambda(ii,jj,:));
            end
        end
    end
    
    if (ST.params.synthetic_camera_params.integrated_opt == 1)
        Psyn_L1 = sum(abs(ST.data.(['sp' num2str(ss)]).PTot_pplane),3);
        
        np_L2 = sum(ST.data.(['sp' num2str(ss)]).np_pplane,3);
        Psyn_L2 = sum(ST.data.(['sp' num2str(ss)]).Psyn_pplane,3);
    else
        for ii=1:NR
            for jj=1:NZ
                if numel(size(ST.data.(['sp' num2str(ss)]).Psyn_pplane)) == 4
                    Psyn_L2_lambda(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_pplane(ii,jj,i1:i2,:),4);
                    np_L2_lambda(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_pplane(ii,jj,i1:i2,:),4);
                    Psyn_L1(ii,jj) = sum(abs(ST.data.(['sp' num2str(ss)]).PTot_pplane(ii,jj,:)),3);
                else
                    Psyn_L2_lambda(ii,jj,:) = ST.data.(['sp' num2str(ss)]).Psyn_pplane(ii,jj,i1:i2);
                    np_L2_lambda(ii,jj,:) = ST.data.(['sp' num2str(ss)]).np_pplane(ii,jj,i1:i2);
                    Psyn_L1(ii,jj) = abs(ST.data.(['sp' num2str(ss)]).PTot_pplane(ii,jj));
                end
                
                np_L2(ii,jj) = sum(np_L2_lambda(ii,jj,:),3);
                Psyn_L2(ii,jj) = trapz(lambda(i1:i2),Psyn_L2_lambda(ii,jj,:));
            end
        end
    end

    if isfield(ST.params,'avalanche_pdf_params')
%         Psyn_avg = averagedSpectrum(ST);
        Psyn_avg = zeros(1,Nl);
    else
        Psyn_sp = singleParticleSpectrum(ST,lambda(i1:i2),...
        ST.params.species.go(ss),deg2rad(ST.params.species.etao(ss)));
    end
    
    % Convert from m to nm
    axis_lambda = 1E9*lambda(i1:i2);
    if (ST.params.synthetic_camera_params.integrated_opt == 0)
        Psyn_L4_lambda = 1E-9*Psyn_L4_lambda;
    end
    
    
    h = figure;
    h.Position(3:4) = [800 860];
    
    nt = 100;
    t = linspace(0,2*pi,nt);
    x = ST.params.fields.Ro + ST.params.fields.a*cos(t);
    y = ST.params.fields.a*sin(t);
    
    xpixel = zeros(1,nt);
    ypixel = zeros(1,nt);
    
    Rc = ST.params.synthetic_camera_params.position(1);
    Zc = ST.params.synthetic_camera_params.position(2);
    f = ST.params.synthetic_camera_params.focal_length;
    incline = deg2rad(ST.params.synthetic_camera_params.incline);
    
    % Last closed surface
    m = -tan(pi/2 - incline);
    n = tan(incline);
    
    xo = m*Rc/(m - n);
    yo = xo*tan(incline);
    Ro = sqrt(xo^2 + yo^2);
    
    d = sqrt( (xo - Rc)^2 + yo^2 );
    xc = -f*Ro/d;
    yc = -f*Zc/d;  

    xperp = zeros(1,nt);
    yperp = zeros(1,nt);
    
    for tt=1:nt
        xtmp = cos(incline)*x(tt);
        ytmp = sin(incline)*x(tt);
        d = sqrt( (xtmp - Rc)^2 + ytmp^2 );
        
        yperp(tt)= -f*y(tt)/d;
        xperp(tt) = f*(x(tt)-Ro)/d - xc;
    end
    % Last closed surface
    
    % Inner wall
    niw = 25;
    
    tiw = linspace(pi/2,3*pi/2,nt);
    xiw = ST.params.fields.Ro + ST.params.fields.a*cos(tiw);
    yiw = ST.params.fields.a*sin(tiw);
    
    Xiwo = zeros(2,nt);
    Xiw = zeros(2,nt,niw);
       
    for tt=1:nt
        xtmp = cos(incline)*xiw(tt);
        ytmp = sin(incline)*xiw(tt);
        d = sqrt( (xtmp - Rc)^2 + ytmp^2 );
        
        Xiwo(2,tt)= -f*yiw(tt)/d;
        Xiwo(1,tt) = f*(xiw(tt)-Ro)/d;
    end

    xtmp = Xiwo(1,:) - xc;

    rotation_angle = linspace(0,pi,niw);
    
    for ii=1:niw
        Xiw(1,:,ii) = xtmp*cos(rotation_angle(ii));
        Xiw(2,:,ii) = Xiwo(2,:);
    end
    % Inner wall
    
    % Magnetic axis
    xtmp = cos(incline)*ST.params.fields.Ro;
    ytmp = sin(incline)*ST.params.fields.Ro;
    d = sqrt( (xtmp - Rc)^2 + ytmp^2 );
    
    ymag_axis = 0;  
    xmag_axis = f*(ST.params.fields.Ro-Ro)/d - xc;
    % Magnetic axis

    % Initial condition
    xic = ST.params.species.Ro(ss) + ST.params.species.r(ss)*cos(t);
    yic = ST.params.species.r(ss)*sin(t);
        
    for tt=1:nt
        xtmp = cos(incline)*xic(tt);
        ytmp = sin(incline)*xic(tt);
        d = sqrt( (xtmp - Rc)^2 + ytmp^2 );
        
        ypixel(tt)= -f*yic(tt)/d;
        xpixel(tt) = f*(xic(tt)-Ro)/d - xc;
    end
    % Initial condition
    
    A = Psyn_L1';
    minval = min(min(A));
    maxval = 0.9*max(max(A));
    v = linspace(minval,maxval,25);
        
    figure(h);
    subplot(4,2,1)
    contourf(RAxis,ZAxis,A,v,'LineStyle','none')
    hold on;plot(x,y,'w','Linewidth',2);hold off
    hold on;plot(ST.params.fields.Ro,0,'wo','Markersize',3,'LineWidth',1,'MarkerFaceColor',[1,1,1],...
        'MarkerEdgeColor',[0.6,0.6,0.6]);hold off
    colormap(jet); hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    xlabel(hc,'$P_{Tot}$ (Watts)','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$Z$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$R$-axis','FontSize',12,'Interpreter','latex')
    

    figure(h);  
    subplot(4,2,2)
    if (ST.params.synthetic_camera_params.integrated_opt == 0)
        f_L4 = squeeze(sum(sum(Psyn_L4_lambda,1),2));
        f_L4 = f_L4/max(f_L4);
        f_L3 = squeeze(sum(sum(Psyn_L3_lambda,1),2));
        f_L3 = f_L3/max(f_L3);
    end
    f_L2 = ST.data.(['sp' num2str(ss)]).P_lambda(:,ss);
%     figure;plot(axis_lambda,f_L2);
    f_L2 = f_L2/max(f_L2);
    if isfield(ST.params,'avalanche_pdf_params')
        P_theory = Psyn_avg/max(Psyn_avg);
    else
        P_theory = Psyn_sp/max(Psyn_sp);
    end
    figure(h)
    try
        plot(axis_lambda,P_theory,'k',axis_lambda,f_L2,'r',axis_lambda,f_L3,'b',axis_lambda,f_L4,'g','LineWidth',1)
    catch
        plot(axis_lambda,P_theory,'k',axis_lambda,f_L2,'r','LineWidth',1)
    end
    ylabel('$P_{syn}$ (A.U.)','FontSize',12,'Interpreter','latex')
    xlim([min(axis_lambda) max(axis_lambda)])
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
    legend({'$\mathcal{P}_{R}(\lambda)$ (Theory)','$\mathcal{P}_{R}(\lambda)$',...
        '$P_R(\lambda)$','$P_R(\lambda,\psi,\chi)$'},'Interpreter','latex')
    
    
    A = Psyn_L2';
    minval = min(min(A));
    maxval = 0.7*max(max(A));
    v = linspace(minval,maxval,25);
    
    figure(h);
    subplot(4,2,3)
    contourf(RAxis,ZAxis,A,v,'LineStyle','none')
    hold on;plot(x,y,'w','Linewidth',2);hold off
    hold on;plot(ST.params.fields.Ro,0,'wo','Markersize',3,'LineWidth',1,'MarkerFaceColor',[1,1,1],...
        'MarkerEdgeColor',[0.6,0.6,0.6]);hold off
    colormap(jet); hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    xlabel(hc,'$P_{syn}$ (Photon/s)','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$Z$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$R$-axis','FontSize',12,'Interpreter','latex')
    
    
    A = np_L2';
    minval = min(min(A));
    maxval = 0.7*max(max(A));
    v = linspace(minval,maxval,25);
    
    figure(h);
    subplot(4,2,4)
    contourf(RAxis,ZAxis,A,v,'LineStyle','none')
    hold on;plot(x,y,'w','Linewidth',2);hold off
    hold on;plot(ST.params.fields.Ro,0,'wo','Markersize',3,'LineWidth',1,'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0.6,0.6,0.6]);hold off
    colormap(jet);  hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    xlabel(hc,'Number of RE','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$Z$-axis','FontSize',14,'Interpreter','latex')
    xlabel('$R$-axis','FontSize',14,'Interpreter','latex')
    
    
    rescaled_xAxis = ST.params.fields.Ro*(xAxis - xc)/xmag_axis;
    rescaled_yAxis = ST.params.fields.Ro*yAxis/xmag_axis;
    
    rescaled_xpixel = ST.params.fields.Ro*xpixel/xmag_axis;
    rescaled_ypixel = ST.params.fields.Ro*ypixel/xmag_axis;
    
    rescaled_xperp = ST.params.fields.Ro*xperp/xmag_axis;
    rescaled_yperp = ST.params.fields.Ro*yperp/xmag_axis;
    
    rescaled_Xiw = ST.params.fields.Ro*(Xiw)/xmag_axis;
      
    % % % % % Camera diagnostic % % % % % 
    
    A = Psyn_L3';
    minval = min(min(A));
    maxval = 0.7*max(max(A));
    v = linspace(minval,maxval,25);
    
    figure(h);
    subplot(4,2,5)
    contourf(rescaled_xAxis,rescaled_yAxis,A,v,'LineStyle','none')
    hold on;plot(rescaled_xpixel,rescaled_ypixel,'w','Linewidth',0.5);hold off
    hold on;plot(rescaled_xperp,rescaled_yperp,-rescaled_xperp,rescaled_yperp,'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    for ii=1:niw
        hold on;plot(rescaled_Xiw(1,:,ii),rescaled_Xiw(2,:,ii),'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    end
    a = max(rescaled_yperp);angle = pi/4;b = a*sin(angle);xb = a*cos(angle);
    rescaled_wxAxisn = ST.params.fields.Ro-xb;
    rescaled_wxAxisp = ST.params.fields.Ro+xb;
    hold on;plot([rescaled_wxAxisn,-rescaled_wxAxisp],[b,b],[rescaled_wxAxisn,-rescaled_wxAxisp],[-b,-b],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    angle = pi/8;b = a*sin(angle);xb = a*cos(angle);
    rescaled_wxAxisn = ST.params.fields.Ro-xb;
    rescaled_wxAxisp = ST.params.fields.Ro+xb;
    hold on;plot([rescaled_wxAxisn,-rescaled_wxAxisp],[b,b],[rescaled_wxAxisn,-rescaled_wxAxisp],[-b,-b],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    hold on;plot([ST.params.fields.Ro,-ST.params.fields.Ro],[a,a],[ST.params.fields.Ro,-ST.params.fields.Ro],[-a,-a],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    hold on;plot(ST.params.fields.Ro,0,'wo',0,ST.params.fields.Ro*yc/xmag_axis,'wo','Markersize',3,'LineWidth',1,'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0.6,0.6,0.6]);hold off
    cm = colormap(jet);cm(1,:) = [0,0,0];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    ax = gca;ax.Color = cm(1,:);ax.ClippingStyle = 'rectangle';
    xlabel(hc,'$P_{syn}$ (Photon/s)','Interpreter','latex','FontSize',12)
    ymin=min(rescaled_yAxis);ymax=max(rescaled_yAxis);xmin=ST.params.fields.Ro-abs(ymin);xmax=ST.params.fields.Ro+abs(ymax);
    axis([xmin, xmax, ymin, ymax]);
    box on; axis square;axis ij;
    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
    
   
    A = np_L3';
    minval = min(min(A));
    maxval = 0.6*max(max(A));
    v = linspace(minval,maxval,25);
    
    figure(h);
    subplot(4,2,6)
    contourf(rescaled_xAxis,rescaled_yAxis,A,v,'LineStyle','none')
    hold on;plot(rescaled_xpixel,rescaled_ypixel,'w','Linewidth',0.5);hold off
    hold on;plot(rescaled_xperp,rescaled_yperp,-rescaled_xperp,rescaled_yperp,'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    for ii=1:niw
        hold on;plot(rescaled_Xiw(1,:,ii),rescaled_Xiw(2,:,ii),'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    end
    a = max(rescaled_yperp);angle = pi/4;b = a*sin(angle);xb = a*cos(angle);
    rescaled_wxAxisn = ST.params.fields.Ro-xb;
    rescaled_wxAxisp = ST.params.fields.Ro+xb;
    hold on;plot([rescaled_wxAxisn,-rescaled_wxAxisp],[b,b],[rescaled_wxAxisn,-rescaled_wxAxisp],[-b,-b],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    angle = pi/8;b = a*sin(angle);xb = a*cos(angle);
    rescaled_wxAxisn = ST.params.fields.Ro-xb;
    rescaled_wxAxisp = ST.params.fields.Ro+xb;
    hold on;plot([rescaled_wxAxisn,-rescaled_wxAxisp],[b,b],[rescaled_wxAxisn,-rescaled_wxAxisp],[-b,-b],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    hold on;plot([ST.params.fields.Ro,-ST.params.fields.Ro],[a,a],[ST.params.fields.Ro,-ST.params.fields.Ro],[-a,-a],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    hold on;plot(ST.params.fields.Ro,0,'wo',0,ST.params.fields.Ro*yc/xmag_axis,'wo','Markersize',3,'LineWidth',1,'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0.6,0.6,0.6]);hold off
    cm = colormap(jet);cm(1,:) = [0,0,0];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    ax = gca;ax.Color = cm(1,:);ax.ClippingStyle = 'rectangle';
    xlabel(hc,'Number of RE','Interpreter','latex','FontSize',12)
    ymin=min(rescaled_yAxis);ymax=max(rescaled_yAxis);xmin=ST.params.fields.Ro-abs(ymin);xmax=ST.params.fields.Ro+abs(ymax);
    axis([xmin, xmax, ymin, ymax]);
    box on; axis square;axis ij;
    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
    
    A = Psyn_L4';
    minval = min(min(A));
    maxval = 0.7*max(max(A));
    v = linspace(minval,maxval,25);
    
    figure(h);
    subplot(4,2,7)
    contourf(rescaled_xAxis,rescaled_yAxis,A,v,'LineStyle','none')
    hold on;plot(rescaled_xpixel,rescaled_ypixel,'w','Linewidth',0.5);hold off
    hold on;plot(rescaled_xperp,rescaled_yperp,-rescaled_xperp,rescaled_yperp,'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    for ii=1:niw
        hold on;plot(rescaled_Xiw(1,:,ii),rescaled_Xiw(2,:,ii),'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    end
    a = max(rescaled_yperp);angle = pi/4;b = a*sin(angle);xb = a*cos(angle);
    rescaled_wxAxisn = ST.params.fields.Ro-xb;
    rescaled_wxAxisp = ST.params.fields.Ro+xb;
    hold on;plot([rescaled_wxAxisn,-rescaled_wxAxisp],[b,b],[rescaled_wxAxisn,-rescaled_wxAxisp],[-b,-b],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    angle = pi/8;b = a*sin(angle);xb = a*cos(angle);
    rescaled_wxAxisn = ST.params.fields.Ro-xb;
    rescaled_wxAxisp = ST.params.fields.Ro+xb;
    hold on;plot([rescaled_wxAxisn,-rescaled_wxAxisp],[b,b],[rescaled_wxAxisn,-rescaled_wxAxisp],[-b,-b],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    hold on;plot([ST.params.fields.Ro,-ST.params.fields.Ro],[a,a],[ST.params.fields.Ro,-ST.params.fields.Ro],[-a,-a],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    hold on;plot(ST.params.fields.Ro,0,'wo',0,ST.params.fields.Ro*yc/xmag_axis,'wo','Markersize',3,'LineWidth',1,'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0.6,0.6,0.6]);hold off
    cm = colormap(jet);cm(1,:) = [0,0,0];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    ax = gca;ax.Color = cm(1,:);ax.ClippingStyle = 'rectangle';
    xlabel(hc,'$P_{syn}$ (Photon/s)','Interpreter','latex','FontSize',12)
    ymin=min(rescaled_yAxis);ymax=max(rescaled_yAxis);xmin=ST.params.fields.Ro-abs(ymin);xmax=ST.params.fields.Ro+abs(ymax);
    axis([xmin, xmax, ymin, ymax]);
    box on; axis square;axis ij;
    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
    
   
    A = np_L4';
    minval = min(min(A));
    maxval = 0.7*max(max(A));
    v = linspace(minval,maxval,25);
    
    figure(h);
    subplot(4,2,8)
    contourf(rescaled_xAxis,rescaled_yAxis,A,v,'LineStyle','none')
    hold on;plot(rescaled_xpixel,rescaled_ypixel,'w','Linewidth',0.5);hold off
    hold on;plot(rescaled_xperp,rescaled_yperp,-rescaled_xperp,rescaled_yperp,'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    for ii=1:niw
        hold on;plot(rescaled_Xiw(1,:,ii),rescaled_Xiw(2,:,ii),'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    end
    a = max(rescaled_yperp);angle = pi/4;b = a*sin(angle);xb = a*cos(angle);
    rescaled_wxAxisn = ST.params.fields.Ro-xb;
    rescaled_wxAxisp = ST.params.fields.Ro+xb;
    hold on;plot([rescaled_wxAxisn,-rescaled_wxAxisp],[b,b],[rescaled_wxAxisn,-rescaled_wxAxisp],[-b,-b],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    angle = pi/8;b = a*sin(angle);xb = a*cos(angle);
    rescaled_wxAxisn = ST.params.fields.Ro-xb;
    rescaled_wxAxisp = ST.params.fields.Ro+xb;
    hold on;plot([rescaled_wxAxisn,-rescaled_wxAxisp],[b,b],[rescaled_wxAxisn,-rescaled_wxAxisp],[-b,-b],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    hold on;plot([ST.params.fields.Ro,-ST.params.fields.Ro],[a,a],[ST.params.fields.Ro,-ST.params.fields.Ro],[-a,-a],...
        'Color',[0.7,0.7,0.7],'Linewidth',1);hold off
    hold on;plot(ST.params.fields.Ro,0,'wo',0,ST.params.fields.Ro*yc/xmag_axis,'wo','Markersize',3,'LineWidth',1,'MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor',[0.6,0.6,0.6]);hold off
    cm = colormap(jet);cm(1,:) = [0,0,0];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    ax = gca;ax.Color = cm(1,:);ax.ClippingStyle = 'rectangle';
    xlabel(hc,'Number of RE','Interpreter','latex','FontSize',12)
    ymin=min(rescaled_yAxis);ymax=max(rescaled_yAxis);xmin=ST.params.fields.Ro-abs(ymin);xmax=ST.params.fields.Ro+abs(ymax);
    axis([xmin, xmax, ymin, ymax]);
    box on; axis square;axis ij;
    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')    
    
    saveas(h,[ST.path 'SyntheticCamera_ss_' num2str(ss)],'fig')
end
end