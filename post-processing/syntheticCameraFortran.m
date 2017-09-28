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

plotCameraSnapshots(ST);

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

filename = [ST.path 'synthetic_camera_snapshots.h5'];
H = h5info(filename);
it = zeros(1,size(H.Groups,1));
for ii=1:size(H.Groups,1)
    tmpstr = strrep(H.Groups(ii).Name,'/','');
    it(ii) = str2num(tmpstr);
end
it = sort(it,'ascend');

% it = ST.range(1):1:ST.range(2);

NX = ST.params.synthetic_camera_params.num_pixels(1);
NY = ST.params.synthetic_camera_params.num_pixels(2);
Nl = ST.params.synthetic_camera_params.Nlambda;


if ~isfield(ST.params.synthetic_camera_params,'toroidal_sections')
    Ntor = 0;
elseif (ST.params.synthetic_camera_params.toroidal_sections == 1)
    Ntor = ST.params.synthetic_camera_params.ntor_sections;
else
    Ntor = 0;
end

filename = [ST.path 'synthetic_camera_snapshots.h5'];

if (ST.params.synthetic_camera_params.integrated_opt == 0)
    list = {'Psyn_angular_pixel','np_angular_pixel','Psyn_lambda_pixel','np_lambda_pixel','Psyn_pplane','np_pplane'};
    
    for ll=1:length(list)
        disp(['Loading ' list{ll}])
        for ss=1:ST.params.simulation.num_species
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,Nl,ST.num_snapshots);
            for ii=1:ST.num_snapshots
                dataset = ...
                    ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                data.(['sp' num2str(ss)]).(list{ll})(:,:,:,ii) = h5read(filename, dataset);
            end
        end
    end
end


if (ST.params.synthetic_camera_params.integrated_opt == 1)
    list={'Psyn_angular_pixel','np_angular_pixel','Psyn_lambda_pixel','np_lambda_pixel','Psyn_pplane','np_pplane',...
        'PTot_pplane','P_lambda','np_lambda','P_a_pixel','P_l_pixel','np_pixel'};
else
    list={'PTot_pplane'};
end

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    
    for ss=1:ST.params.simulation.num_species
        try
            if strcmp(list{ll},'P_lambda')
                data.(['sp' num2str(ss)]).(list{ll}) = zeros(Nl,ST.num_snapshots);
                for ii=1:ST.num_snapshots % Here
                    dataset = ...
                        ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                        '/' list{ll}];
                    data.(['sp' num2str(ss)]).(list{ll})(:,ii) = h5read(filename, dataset);
                end
            elseif (strcmp(list{ll},'P_lambda') || strcmp(list{ll},'P_a_pixel') || strcmp(list{ll},'P_l_pixel'))
                if (ST.params.synthetic_camera_params.toroidal_sections == 1)
                    data.(['sp' num2str(ss)]).(list{ll}) = zeros(Nl,Ntor,ST.num_snapshots);
                else
                    data.(['sp' num2str(ss)]).(list{ll}) = zeros(Nl,ST.num_snapshots);
                end
                for ii=1:ST.num_snapshots % Here
                    dataset = ...
                        ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                        '/' list{ll}];
                    if (ST.params.synthetic_camera_params.toroidal_sections == 1)
                        data.(['sp' num2str(ss)]).(list{ll})(:,:,ii) = h5read(filename, dataset);
                    else
                        data.(['sp' num2str(ss)]).(list{ll})(:,ii) = h5read(filename, dataset);
                    end
                end
            elseif (strcmp(list{ll},'np_pixel') || strcmp(list{ll},'np_lambda'))
                data.(['sp' num2str(ss)]).(list{ll}) = zeros(1,ST.num_snapshots);
                for ii=1:ST.num_snapshots % Here
                    dataset = ...
                        ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                        '/' list{ll}];
                    data.(['sp' num2str(ss)]).(list{ll})(ii) = h5read(filename, dataset);
                end
            elseif (strcmp(list{ll},'Psyn_pplane') || strcmp(list{ll},'np_pplane') || strcmp(list{ll},'PTot_pplane'))
                data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,ST.num_snapshots);
                for ii=1:ST.num_snapshots % Here
                    dataset = ...
                        ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                        '/' list{ll}];
                    
                    data.(['sp' num2str(ss)]).(list{ll})(:,:,ii) = h5read(filename, dataset);
                end
            else
                if (ST.params.synthetic_camera_params.toroidal_sections == 1)
                    data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,Ntor,ST.num_snapshots);
                else
                    data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,ST.num_snapshots);
                end
                
                for ii=1:ST.num_snapshots % Here
                    dataset = ...
                        ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                        '/' list{ll}];
                    if (ST.params.synthetic_camera_params.toroidal_sections == 1)
                        data.(['sp' num2str(ss)]).(list{ll})(:,:,:,ii) = h5read(filename, dataset);
                    else
                        data.(['sp' num2str(ss)]).(list{ll})(:,:,ii) = h5read(filename, dataset);
                    end
                end
            end
        catch
            data.(['sp' num2str(ss)]).(list{ll}) = [];
        end
    end
    
end


end

function Psyn = singleParticleSpectrum(ST,lambda,g,eta)
Psyn = zeros(size(lambda));

q = abs(ST.params.species.q(1));
m = ST.params.species.m(1);
c = ST.params.scales.v;
v = c*sqrt(1-1/g^2);
ep = 8.854E-12;% Electric permitivity

k = q*ST.params.fields_and_profiles.Bo*sin(eta)/(g*m*v);
l = lambda;
lc = 4*pi/(3*k*g^3);

z = lc./l;

% figure;plot(1E9*l,z)

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

function Psyn = averagedSpectrum(ST,Np,Nchi)

narginchk(1,3);

if (nargin == 1)
    Np = 50;
    Nchi = 100;
end

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

fRE = @(p,x) (Ehat/Cz)*p.*exp( -p.*(x/Cz + 0.5*Ehat*(1 - x.^2)./x) )./x;

Fo = integral2(fRE,pmin,pmax,chimin,1);
% Fo = integral2(fRE,0,500,0,1);

f = @(p,x) fRE(p,x)/Fo;

p = linspace(pmin,pmax,Np);
pitch = linspace(0,pitchmax,Nchi);
chi = cos(pitch);

l = ST.params.synthetic_camera_params.lambda;
Psyn = zeros(size(l));
Psyn_p_chi = zeros(numel(l),Np,Nchi);
Psyn_p = zeros(numel(l),Np);
Psyn_chi = zeros(numel(l),Nchi);

% for ll=1:numel(l)
%     for pp=1:Np
%         for cc=1:Nchi
%             Psp = singleParticleSpectrum(ST,l(ll),g(p(pp)),eta(chi(cc)));
%             Psyn_p_chi(ll,pp,cc) = f(p(pp),chi(cc))*Psp;
%         end
%         Psyn_p(ll,pp) = trapz(fliplr(chi),squeeze(Psyn_p_chi(ll,pp,:)));
%     end
%     Psyn(ll) = trapz(p,squeeze(Psyn_p(ll,:)));
% end


for ll=1:numel(l)
    for cc=1:Nchi
        for pp=1:Np
            Psp = singleParticleSpectrum(ST,l(ll),g(p(pp)),eta(chi(cc)));
            Psyn_p_chi(ll,pp,cc) = f(p(pp),chi(cc))*Psp;
        end
        Psyn_chi(ll,cc) = trapz(p,squeeze(Psyn_p_chi(ll,:,cc)));
    end
    Psyn(ll) = trapz(fliplr(chi),squeeze(Psyn_chi(ll,:)));
end

E = (g(p)*m*c^2/q)/1E6; % MeV
xAxis = rad2deg(pitch);
lAxis = l/1E-9;

% figure;
% surf(E,lAxis,squeeze(Psyn_p),'LineStyle','none');
% colormap(jet(1024));
% ylabel('$\lambda$ (nm)','Interpreter','latex')
% xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex')


% figure
% plot(lAxis,Psyn)
% xlabel('$\lambda$ (nm)','Interpreter','latex')
% ylabel('$P_{R}$ (Watts)','Interpreter','latex')
end

function plotCameraSnapshots(ST)
% Psyn_L1 = PTot
% Psyn_L2 = Psyn(lambda) in poloidal plane
% Psyn_L3 = Psyn(lambda) in pixel plane
% Psyn_L4 = Psyn(lambda,psi,chi) in pixel plane
disp('Plotting snapshots...')

xRectangles = [1.41,1.619; 1.537,1.66; 1.943,2.077];
yRectangles = [-0.1075,0.1161; -0.0626,0.07147; 0.1831,0.363];
colorRectangles = [0,0,0; 1,0,0; 0,0,1];
plotToroidalSections = true;
figuresToShare = false;


lambda = ST.params.synthetic_camera_params.lambda;
xAxis = ST.params.synthetic_camera_params.pixels_nodes_x;
yAxis = ST.params.synthetic_camera_params.pixels_nodes_y;

NX = ST.params.synthetic_camera_params.num_pixels(1);
NY = ST.params.synthetic_camera_params.num_pixels(2);

RAxis = ST.params.poloidal_plane_params.nodes_R;
ZAxis = ST.params.poloidal_plane_params.nodes_Z;
NR = ST.params.poloidal_plane_params.grid_dims(1);
NZ = ST.params.poloidal_plane_params.grid_dims(1);


if ~isfield(ST.params.synthetic_camera_params,'toroidal_sections')
    NT = 0;
elseif (ST.params.synthetic_camera_params.toroidal_sections == 1)
    NT = ST.params.synthetic_camera_params.ntor_sections;
else
    NT = 0;
end


if isempty(ST.lambdas)
    ST.lambdas = [min(lambda) max(lambda)];
end

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
    
    if (NT ~= 0)
        Psyn_L4 = zeros(NX,NY,NT);
        np_L4 = zeros(NX,NY,NT);
        
        Psyn_L3 = zeros(NX,NY,NT);
        np_L3 = zeros(NX,NY,NT);
    else
        Psyn_L4 = zeros(NX,NY);
        np_L4 = zeros(NX,NY);
        
        Psyn_L3 = zeros(NX,NY);
        np_L3 = zeros(NX,NY);
    end
    
    Psyn_L2 = zeros(NR,NZ);
    np_L2 = zeros(NR,NZ);
    
    Psyn_L1 = zeros(NR,NZ);
    
    if (ST.params.synthetic_camera_params.integrated_opt == 0)
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
    else
        Psyn_L1 = sum(abs(ST.data.(['sp' num2str(ss)]).PTot_pplane),3);
        
        np_L2 = sum(ST.data.(['sp' num2str(ss)]).np_pplane,3);
        Psyn_L2 = sum(ST.data.(['sp' num2str(ss)]).Psyn_pplane,3);
    end
    
    if (ST.params.synthetic_camera_params.integrated_opt == 0)
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
    else
        npl = ST.data.(['sp' num2str(ss)]).np_lambda;
        
        P_L2 = squeeze(sum(ST.data.(['sp' num2str(ss)]).P_lambda(i1:i2,:),2))/sum(sum(np_L2));        
        
        if (NT ~= 0)
            np = ST.data.(['sp' num2str(ss)]).np_pixel;
            
            if ndims(ST.data.(['sp' num2str(ss)]).P_l_pixel) == 2
                P_L3 = ST.data.(['sp' num2str(ss)]).P_l_pixel(i1:i2,:);
            else
                P_L3 = squeeze(sum(ST.data.(['sp' num2str(ss)]).P_l_pixel(i1:i2,:,:),3));
            end
            
            if ndims(ST.data.(['sp' num2str(ss)]).P_a_pixel) == 2
                P_L4 = ST.data.(['sp' num2str(ss)]).P_a_pixel(i1:i2,:);
            else
                P_L4 = squeeze(sum(ST.data.(['sp' num2str(ss)]).P_a_pixel(i1:i2,:,:),3));
            end
            
            np_L3 = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel,4);
            np_L4 = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel,4);
            
            Psyn_L3 = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel,4);
            Psyn_L4 = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel,4);
            
        else
            P_L3 = squeeze(sum(ST.data.(['sp' num2str(ss)]).P_l_pixel(i1:i2,:),2));

            P_L4 = squeeze(sum(ST.data.(['sp' num2str(ss)]).P_a_pixel(i1:i2,:),2));
            
            np_L3 = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel,3);
            np_L4 = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel,3);
            
            Psyn_L3 = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel,3);
            Psyn_L4 = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel,3);
        end
    end
    
    if isfield(ST.params,'avalanche_pdf_params')
        Psyn_avg = zeros(size(P_L4));%averagedSpectrum(ST,40,50);
    else
        Psyn_sp = singleParticleSpectrum(ST,lambda(i1:i2),...
            ST.params.species.go(ss),deg2rad(ST.params.species.etao(ss)));
    end
    
    % Convert from m to nm
    axis_lambda = 1E9*lambda(i1:i2);
    if (ST.params.synthetic_camera_params.integrated_opt == 0)
        Psyn_L4_lambda = 1E-9*Psyn_L4_lambda;
    else
        P_L2 = 1E-9*P_L2;
        P_L3 = 1E-9*P_L3;
        P_L4 = 1E-9*P_L4;
    end
    
    nt = 100;
    t = linspace(0,2*pi,nt);
    if isfield(ST.params.fields_and_profiles,'a')
        a95 = ST.params.fields_and_profiles.a;
    else
        a95 = max(ST.params.fields_and_profiles.R) - ST.params.fields_and_profiles.Ro;
    end
    x = ST.params.fields_and_profiles.Ro + a95*cos(t);
    y = a95*sin(t);
    
    xpixel = zeros(1,nt);
    ypixel = zeros(1,nt);
    
    Rc = ST.params.synthetic_camera_params.position(1);
    Zc = ST.params.synthetic_camera_params.position(2);
    f = ST.params.synthetic_camera_params.focal_length;
    incline = deg2rad(ST.params.synthetic_camera_params.incline);
    Rc_zero = Rc*sin(pi/2 - incline);
    angle = atan(abs(xAxis(1))/f);
    Rd = Rc_zero*cot(pi/2-incline);
    Rc_min = Rc_zero - Rd*tan(angle);
    Rc_max = Rc_zero + Rd*tan(angle);
    
    scaling_factor = (Rc_max - Rc_min)/(max(xAxis) - min(xAxis));
    xAxis_rescaled = scaling_factor*xAxis + Rc_zero;
    yAxis_rescaled = scaling_factor*yAxis + Zc;
    
%     xAxis_rescaled = xAxis;
%     yAxis_rescaled = yAxis;
    
    if (NT ~= 0)
        f_L2 = P_L2;
        f_L3 = P_L3;
        f_L4 = P_L4;
        if isfield(ST.params,'avalanche_pdf_params')
            P_theory = 1E-9*Psyn_avg; %P_theory = Psyn_avg/max(Psyn_avg);
        else
            P_theory = 1E-9*Psyn_sp; % P_theory = Psyn_sp/max(Psyn_sp);
        end
        
        fig = figure;
        subplot(3,1,1)
        plot(axis_lambda,f_L2,'k','LineWidth',2)
        hold on;plot(axis_lambda,P_theory,'r');hold off
        ylabel('$P_R(\lambda)$ (Watts)','FontSize',12,'Interpreter','latex')
        xlim([min(axis_lambda) max(axis_lambda)]);box on;grid on;
        xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
        
        subplot(3,1,2)
        plot(axis_lambda,f_L3,'--',axis_lambda,sum(f_L3,2),'LineWidth',2)
        ylabel('$P_R(\lambda)$ (Watts)','FontSize',12,'Interpreter','latex')
        xlim([min(axis_lambda) max(axis_lambda)]);box on;grid on;
        xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
        
        subplot(3,1,3)
        plot(axis_lambda,f_L4,'--',axis_lambda,sum(f_L4,2),'LineWidth',2)
        ylabel('$P_R(\lambda)$ (Watts)','FontSize',12,'Interpreter','latex')
        xlim([min(axis_lambda) max(axis_lambda)]);box on;grid on;
        xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
        
%         saveas(fig,[ST.path 'Spectra_ss_' num2str(ss)],'fig')
        
        Dtor = 360/double(NT);
        if plotToroidalSections
            for tt=1:NT
                h = figure;
                h.Position(3:4) = [800 860];
                annotation('textbox',[0.03 0.03 0.5 0.05],'String',['$\phi\in$[' num2str((tt-1)*Dtor) '$^\circ$,' num2str(tt*Dtor) '$^\circ$]'],...
                    'FitBoxToText','on','Interpreter','latex');
                
                % % % % % Camera diagnostic % % % % %
                
                A = Psyn_L3(:,:,tt)';
                minval = min(min(A(A~=0)));
                maxval = 0.8*max(max(A));
                
                if ~isempty(minval)
                    v = linspace(minval,maxval,25);                    
                    figure(h);
                    subplot(2,2,1)
                    contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
                    ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
                    axis([xmin, xmax, ymin, ymax]);
                    box on; axis equal
                    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
                    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
                    
                    cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([minval,maxval]);
                    ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
                    if ST.params.synthetic_camera_params.photon_count
                        xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
                    else
                        xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
                    end
                    
                    
                    
                    A = np_L3(:,:,tt)';
                    minval = min(min(A(A~=0)));
                    maxval = 0.8*max(max(A));
                    v = linspace(minval,maxval,25);
                    
                    figure(h);
                    subplot(2,2,2)
                    contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
                    ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
                    axis([xmin, xmax, ymin, ymax]);
                    box on; axis equal;
                    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
                    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
                    
                    cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([minval,maxval]);
                    ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
                    xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
                    
                    
                    A = Psyn_L4(:,:,tt)';
                    minval = min(min(A(A~=0)));
                    maxval = 0.8*max(max(A));
                    v = linspace(minval,maxval,25);
                    
                    figure(h);
                    subplot(2,2,3)
                    contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
                    ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
                    axis([xmin, xmax, ymin, ymax]);
                    box on; axis equal;
                    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
                    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
                    
                    cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([minval,maxval]);
                    ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
                    if ST.params.synthetic_camera_params.photon_count
                        xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
                    else
                        xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
                    end
                    
                    
                    A = np_L4(:,:,tt)';
                    minval = min(min(A(A~=0)));
                    maxval = 0.8*max(max(A));
                    v = linspace(minval,maxval,25);
                    
                    figure(h);
                    subplot(2,2,4)
                    contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
                    ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
                    axis([xmin, xmax, ymin, ymax]);
                    box on; axis equal;
                    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
                    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
                    
                    cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([minval,maxval]);
                    ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
                    xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
                    
                                saveas(h,[ST.path 'SyntheticCamera_ss_' num2str(ss)],'fig')
                else
                    close(h)
                end
            end
        end
        
        h = figure;
        h.Position(3:4) = [1350 815];
        annotation('textbox',[0.03 0.03 0.5 0.05],'String','Total','FitBoxToText','on','Interpreter','latex');
        
        % % % % % Camera diagnostic % % % % %
        
        if ~figuresToShare

            A = Psyn_L2';
            B = reshape(A,[numel(A),1]);
            B(B==0) = [];
            if ST.params.synthetic_camera_params.photon_count
                B(B<1) = [];
            end
            minval = min(B);
            maxval = 3*std(B);
            v = linspace(minval,maxval,50);
            
            figure(h);
            subplot(2,4,1)
            contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
            ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
            if ~figuresToShare
                axis([xmin, xmax, ymin, ymax]);
            else
                axis([0.95 2.45 -0.5 1])
            end
            box on; axis equal;
            ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
            xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
            
            cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
            ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
            if ST.params.synthetic_camera_params.photon_count
                xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
            else
                xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
            end
            
            
            A = np_L2';
            minval = min(min(A(A~=0)));
            maxval = 0.8*max(max(A));
            v = linspace(minval,maxval,25);
            
            figure(h);
            subplot(2,4,5)
            contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
            ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
            if ~figuresToShare
                axis([xmin, xmax, ymin, ymax]);
            else
                axis([0.95 2.45 -0.5 1])
            end
            box on; axis equal;
            ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
            xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
            
            cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
            ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
            xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
            
            
            A = sum(Psyn_L3,3)';
            B = reshape(A,[numel(A),1]);
            B(B==0) = [];
            if ST.params.synthetic_camera_params.photon_count
                B(B<1) = [];
            end
            minval = min(B);
            maxval = 3*std(B);
            v = linspace(minval,maxval,50);
            
            figure(h);
            subplot(2,4,2)
            contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
            ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
            if ~figuresToShare
                axis([xmin, xmax, ymin, ymax]);
            else
                axis([0.95 2.45 -0.5 1])
            end
            box on; axis equal;
            ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
            xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
            
            cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
            ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
            if ST.params.synthetic_camera_params.photon_count
                xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
            else
                xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
            end
            
            
            A = sum(np_L3,3)';
            minval = min(min(A(A~=0)));
            maxval = 0.8*max(max(A));
            v = linspace(minval,maxval,25);
            
            figure(h);
            subplot(2,4,6)
            contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
            ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
            if ~figuresToShare
                axis([xmin, xmax, ymin, ymax]);
            else
                axis([0.95 2.45 -0.5 1])
            end
            box on; axis equal;
            ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
            xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
            
            cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
            ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
            xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        end
        
        A = sum(Psyn_L4,3)';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        if ST.params.synthetic_camera_params.photon_count
            B(B<1) = [];
        end
        minval = min(B);
        maxval = 3*std(B);
        v = linspace(minval,maxval,50);
        
        figure(h);
        if ~figuresToShare
            subplot(2,4,3)
        else
            subplot(2,2,1)
        end
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        box on; axis equal;
        if ~figuresToShare
            axis([xmin, xmax, ymin, ymax]);
        else
            axis([0.95 2.45 -0.5 1])
        end
        ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
        xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
        
        iX = zeros(2,size(xRectangles,1));
        iY = zeros(2,size(xRectangles,1));
        cameraCount = zeros(1,size(xRectangles,1));
        
        figure(h)
        if ~figuresToShare
            subplot(2,4,3)
        else
            subplot(2,2,1)
        end
%         hold on
%         for ii=1:size(xRectangles,1)
%             iX(ii,1) = find(xAxis_rescaled > xRectangles(ii,1),1,'first');
%             iX(ii,2) = find(xAxis_rescaled > xRectangles(ii,2),1,'first') - 1;
%             
%             iY(ii,1) = find(yAxis_rescaled > yRectangles(ii,1),1,'first');
%             iY(ii,2) = find(yAxis_rescaled > yRectangles(ii,2),1,'first') - 1;
%             
%             cameraCount(ii) = sum(sum(A(iY(ii,1):iY(ii,2),iX(ii,1):iX(ii,2))));
%             
%             xplot = [xRectangles(ii,1),xRectangles(ii,2),xRectangles(ii,2),xRectangles(ii,1),xRectangles(ii,1)];
%             yplot = [yRectangles(ii,1),yRectangles(ii,1),yRectangles(ii,2),yRectangles(ii,2),yRectangles(ii,1)];
%             plot(xplot,yplot,'LineWidth',2,'Color',colorRectangles(ii,:))
%         end
%         hold off
               
        A = sum(np_L4,3)';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        figure(h);
        if ~figuresToShare
            subplot(2,4,7)
        else
            subplot(2,2,2)
        end
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        box on; axis equal;
        if ~figuresToShare
            axis([xmin, xmax, ymin, ymax]);
        else
            axis([0.95 2.45 -0.5 1])
        end
        ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
        xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        
        figure(h)
        if ~figuresToShare
            subplot(2,4,[4 8])
        else
            subplot(2,2,[3 4])
        end        
        hold on
        for ii=1:size(xRectangles,1)
            plot(ii,cameraCount(ii),'s','MarkerSize',6,'Color',colorRectangles(ii,:),'MarkerFaceColor',colorRectangles(ii,:))
        end
        hold off
        box on;xlim([0 size(xRectangles,1)+1])
        ylabel('Box count','FontSize',12,'Interpreter','latex')
        xlabel('Box number','FontSize',12,'Interpreter','latex')
        
        saveas(h,[ST.path 'Total_SE_ss_' num2str(ss)],'fig')
        
    else
        h = figure;
        h.Position(3:4) = [800 860];
        
        A = Psyn_L1';
        minval = min(min(A));
        maxval = 0.8*max(max(A));
        v = linspace(minval,maxval,25);
        
        figure(h);
        subplot(4,2,1)
        contourf(RAxis,ZAxis,A,v,'LineStyle','none')
        ax = gca; ax.YDir = 'reverse';
        colormap(jet(1024)); hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
        box on; axis square;view([0 -90])
        ylabel('$Z$-axis','FontSize',12,'Interpreter','latex')
        xlabel('$R$-axis','FontSize',12,'Interpreter','latex')
        
        
        if (ST.params.synthetic_camera_params.integrated_opt == 0)
            % something here
        else
            f_L3 = P_L3;
            f_L4 = P_L4;
            f_L2 = P_L2;
        end
        
        if isfield(ST.params,'avalanche_pdf_params')
            P_theory = 1E-9*Psyn_avg; %P_theory = Psyn_avg/max(Psyn_avg);
        else
            P_theory = 1E-9*Psyn_sp; % P_theory = Psyn_sp/max(Psyn_sp);
        end
        
        figure(h);
        subplot(4,2,2)
        try
            plot(axis_lambda,P_theory,'k',axis_lambda,f_L2,'r',axis_lambda,f_L3,'b',axis_lambda,f_L4,'g','LineWidth',1)
        catch
            figure(h)
            plot(axis_lambda,P_theory(i1:i2),'k',axis_lambda,f_L2,'r','LineWidth',1)
            figure;plot(axis_lambda,f_L2./P_theory(i1:i2),'r','LineWidth',1)
        end
        ylabel('$P_R(\lambda)$ (Watts)','FontSize',12,'Interpreter','latex')
        xlim([min(axis_lambda) max(axis_lambda)])
        xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
        legend({'$\mathcal{P}_{R}(\lambda)$ (Theory)','$\mathcal{P}_{R}(\lambda)$',...
            '$P_R(\lambda)$','$P_R(\lambda,\psi,\chi)$'},'Interpreter','latex')
        
        
        A = Psyn_L2';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        if ST.params.synthetic_camera_params.photon_count
            B(B<1) = [];
        end
        minval = min(B);
        maxval = 3*std(B);
        v = linspace(minval,maxval,50);
        
        figure(h);
        subplot(4,2,3)
        contourf(RAxis,ZAxis,A,v,'LineStyle','none')
        ax = gca; ax.YDir = 'reverse';
        colormap(jet(1024)); hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
        box on; axis square;view([0 -90])
        ylabel('$Z$-axis','FontSize',12,'Interpreter','latex')
        xlabel('$R$-axis','FontSize',12,'Interpreter','latex')
        
        
        A = np_L2';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        figure(h);
        subplot(4,2,4)
        contourf(RAxis,ZAxis,A,v,'LineStyle','none')
        ax = gca; ax.YDir = 'reverse';
        colormap(jet(1024));  hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        box on; axis square;view([0 -90])
        ylabel('$Z$-axis','FontSize',14,'Interpreter','latex')
        xlabel('$R$-axis','FontSize',14,'Interpreter','latex')
        
        
        % % % % % Camera diagnostic % % % % %
        
        A = Psyn_L3';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        if ST.params.synthetic_camera_params.photon_count
            B(B<1) = [];
        end
        minval = min(B);
        maxval = 3*std(B);
        v = linspace(minval,maxval,50);
        
        figure(h);
        subplot(4,2,5)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        axis([xmin, xmax, ymin, ymax]);
        box on; axis equal
        ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
        xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
        
        A = np_L3';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        figure(h);
        subplot(4,2,6)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        axis([xmin, xmax, ymin, ymax]);
        box on; axis equal;
        ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
        xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        
        A = Psyn_L4';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        if ST.params.synthetic_camera_params.photon_count
            B(B<1) = [];
        end
        minval = min(B);
        maxval = 3*std(B);
        v = linspace(minval,maxval,50);
        
        figure(h);
        subplot(4,2,7)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        axis([xmin, xmax, ymin, ymax]);
        box on; axis equal;
        ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
        xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
        
        
        A = np_L4';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        figure(h);
        subplot(4,2,8)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        axis([xmin, xmax, ymin, ymax]);
        box on; axis equal;
        ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
        xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        
%         saveas(h,[ST.path 'SyntheticCamera_ss_' num2str(ss)],'fig')
    end
end
end