function ST = cameraPicsAnalysis(path,lambdas,visible,range)
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

% plotCameraSnapshots(ST);

picsAnalysis(ST);

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


cadence = double(ST.params.simulation.output_cadence);
[~,I1] = min(abs(cadence*ST.range(1)-it));
[~,I2] = min(abs(cadence*ST.range(2)-it));

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
    list = {'Psyn_angular_pixel','np_angular_pixel','Psyn_lambda_pixel','np_lambda_pixel','Psyn_p_pplane','Psyn_t_pplane','np_p_pplane','np_t_pplane'};
    
    for ll=1:length(list)
        disp(['Loading ' list{ll}])
        for ss=1:ST.params.simulation.num_species
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,Nl,ST.num_snapshots);
            for ii=I1:I2
                dataset = ...
                    ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                data.(['sp' num2str(ss)]).(list{ll})(:,:,:,ii) = h5read(filename, dataset);
            end
        end
    end
end


if (ST.params.synthetic_camera_params.integrated_opt == 1)
    list={'Psyn_angular_pixel','np_angular_pixel','Psyn_lambda_pixel','np_lambda_pixel','Psyn_p_pplane','Psyn_t_pplane','np_p_pplane','np_t_pplane'};
end

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    
    for ss=1:ST.params.simulation.num_species
        try
            if (strcmp(list{ll},'Psyn_p_pplane') || strcmp(list{ll},'np_p_pplane') || strcmp(list{ll},'PTot_p_pplane') ||...
                    strcmp(list{ll},'Psyn_t_pplane') || strcmp(list{ll},'np_t_pplane') || strcmp(list{ll},'PTot_t_pplane'))
                if (ST.params.synthetic_camera_params.toroidal_sections == 1)
                    data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,Ntor,ST.num_snapshots);
                else
                    data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,ST.num_snapshots);
                end
                               
                for ii=I1:I2
                    dataset = ...
                        ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                        '/' list{ll}];
                    if (ST.params.synthetic_camera_params.toroidal_sections == 1)
                        data.(['sp' num2str(ss)]).(list{ll})(:,:,:,ii) = h5read(filename, dataset);
                    else
                        data.(['sp' num2str(ss)]).(list{ll})(:,:,ii) = h5read(filename, dataset);
                    end
                end
            else
                if (ST.params.synthetic_camera_params.toroidal_sections == 1)
                    data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,Ntor,ST.num_snapshots);
                else
                    data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,ST.num_snapshots);
                end
                
                for ii=I1:I2
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

function plotCameraSnapshots(ST)
disp('Plotting snapshots...')

plotToroidalSections = true; % <-------------------------------------------
threshold = 4; % <---------------------------------------------------------

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

    if (NT ~= 0)
        PR_f = zeros(NX,NY,NT);
        np_f = zeros(NX,NY,NT);
        
        PR_c = zeros(NX,NY,NT);
        np_c = zeros(NX,NY,NT);
    else
        PR_f = zeros(NX,NY);
        np_f = zeros(NX,NY);
        
        PR_c = zeros(NX,NY);
        np_c = zeros(NX,NY);
    end
    
    if (ST.params.synthetic_camera_params.integrated_opt == 0)
        for ii=1:NX
            for jj=1:NY
                if numel(size(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel)) == 4
                    P_cone(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel(ii,jj,i1:i2,:),4);
                    P_full(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel(ii,jj,i1:i2,:),4);
                    
                    np_cone(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel(ii,jj,i1:i2,:),4);
                    np_full(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel(ii,jj,i1:i2,:),4);
                else
                    P_cone(ii,jj,:) = ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel(ii,jj,i1:i2);
                    P_full(ii,jj,:) = ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel(ii,jj,i1:i2);
                    
                    np_cone(ii,jj,:) = ST.data.(['sp' num2str(ss)]).np_lambda_pixel(ii,jj,i1:i2);
                    np_full(ii,jj,:) = ST.data.(['sp' num2str(ss)]).np_angular_pixel(ii,jj,i1:i2);
                end
                
                np_c(ii,jj) = sum(np_cone(ii,jj,:),3);
                np_f(ii,jj) = sum(np_full(ii,jj,:),3);
                
                PR_c(ii,jj) = trapz(lambda(i1:i2),P_cone(ii,jj,:));
                PR_f(ii,jj) = trapz(lambda(i1:i2),P_full(ii,jj,:));
            end
        end
    else
        if (NT ~= 0)  
            np_c = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel,4);
            np_f = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel,4);
            
            PR_c = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel,4);
            PR_f = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel,4);
        else
            np_c = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel,3);
            np_f = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel,3);
            
            PR_c = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel,3);
            PR_f = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel,3);
        end
    end

    % Convert from m to nm
    axis_lambda = 1E9*lambda(i1:i2);

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
    
    I1 = find(yAxis_rescaled>0.25,1)-1;z1 = yAxis_rescaled(I1);
    I0 = find(yAxis_rescaled>0,1);z0 = yAxis_rescaled(I0);
    I2 = find(yAxis_rescaled>-0.25,1);z2 = yAxis_rescaled(I2);
    
    J0 = find(xAxis_rescaled>1.3,1)-1;r0 = xAxis_rescaled(J0);
    J1 = find(xAxis_rescaled>1.5,1)-1;r1 = xAxis_rescaled(J1);
    J2 = find(xAxis_rescaled>1.7,1)-1;r2 = xAxis_rescaled(J2);
    
    if (NT ~= 0)
        Dtor = 360/double(NT);
        
        if plotToroidalSections
            fig_tor_c = figure;
            fig_tor_f = figure;
            
            nsubs = ceil(sqrt(NT));
            
            for tt=1:NT               
                Ac = PR_c(:,:,tt)';
                Bc = reshape(Ac,[numel(Ac),1]);
                Bc(Bc==0) = [];
                if ST.params.synthetic_camera_params.photon_count
                    Bc(Bc<1) = [];
                end
                
                Af = PR_f(:,:,tt)';
                Bf = reshape(Af,[numel(Af),1]);
                Bf(Bf==0) = [];
                if ST.params.synthetic_camera_params.photon_count
                    Bf(Bf<1) = [];
                end
                
                if (~isempty(Bc) && ~isempty(Bf))
                    minval = min(Bc);
                    maxval = threshold*std(Bc);
                    v = linspace(minval,maxval,50);
                    
                    figure(fig_tor_c);subplot(nsubs,nsubs,tt)
                    subplot(3,2,3)
                    contourf(xAxis_rescaled,yAxis_rescaled,Ac,v,'LineStyle','none')
                    ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
                    axis([xmin, xmax, ymin, ymax]);
                    box on; axis equal
                    ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
                    xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
                    
                    cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([minval,maxval]);
                    ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
                    if ST.params.synthetic_camera_params.photon_count
                        xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
                    else
                        xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
                    end
                    
                    minval = min(Bf);
                    maxval = threshold*std(Bf);
                    v = linspace(minval,maxval,50);
                    
                    figure(fig_tor_t);subplot(nsubs,nsubs,tt)
                    subplot(3,2,5)
                    contourf(xAxis_rescaled,yAxis_rescaled,Af,v,'LineStyle','none')
                    ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
                    axis([xmin, xmax, ymin, ymax]);
                    box on; axis equal;
                    ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
                    xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
                    
                    cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([minval,maxval]);
                    ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
                    if ST.params.synthetic_camera_params.photon_count
                        xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
                    else
                        xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
                    end
                    
                    saveas(fig_tor_c,[ST.path 'toroidal_section_cone_' num2str(tt)],'fig')
                    saveas(fig_tor_f,[ST.path 'toroidal_section_full_' num2str(tt)],'fig')
                else
                    close(h)
                end
            end
        end
        
        fig = figure;
        
        A = sum(PR_c,3)';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        if ST.params.synthetic_camera_params.photon_count
            B(B<1) = [];
        end
        minval = min(B);
        maxval = threshold*std(B);
        v = linspace(minval,maxval,50);
        
        figure(fig);subplot(2,2,1)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        box on; axis equal;
        axis([xmin, xmax, ymin, ymax]);
        ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
               
        A = sum(np_c,3)';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        figure(fig);subplot(2,2,2)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        box on; axis equal;
        axis([xmin, xmax, ymin, ymax]);
        ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        
        
        A = sum(PR_f,3)';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        if ST.params.synthetic_camera_params.photon_count
            B(B<1) = [];
        end
        minval = min(B);
        maxval = threshold*std(B);
        v = linspace(minval,maxval,50);
        
        figure(fig);subplot(2,2,3)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        box on; axis equal;
        axis([xmin, xmax, ymin, ymax]);
        ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
               
        A = sum(np_f,3)';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        figure(fig);subplot(2,2,4)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        box on; axis equal;
        axis([xmin, xmax, ymin, ymax]);
        ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','northoutside');caxis([minval,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        
        saveas(h,[ST.path 'Total_PR_' num2str(ss)],'fig')
    else
        fig = figure;
        
        if (any(xAxis_rescaled<0))
            xAxis_rescaled = -xAxis_rescaled;
        end
        
        A = PR_c';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        if ST.params.synthetic_camera_params.photon_count
            B(B<1) = [];
        end
        minval = min(B);
        maxval = threshold*std(B);
        v = linspace(minval,maxval,50);
        
        figure(fig);subplot(2,2,1)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        axis([xmin, xmax, ymin, ymax]);
        box on; axis equal
        ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
        
        C = max(A(I0,:));
        A0 = A(I0,:)/C;A1 = A(I1,:)/C;A2 = A(I2,:)/C;
        CC = max(A(:,J1));
        AA0 = A(:,J0)/CC;AA1 = A(:,J1)/CC;AA2 = A(:,J2)/CC;
        
        fig_slices=figure;
        subplot(3,2,1)
        plot(xAxis_rescaled,A1,'g--','LineWidth',2)
        axis([1 2.0 0 1]);box on;grid minor;
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        subplot(3,2,3)
        plot(xAxis_rescaled,A0,'r--','LineWidth',2)
        axis([1 2.0 0 1]);box on;grid minor;
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        subplot(3,2,5)
        plot(xAxis_rescaled,A2,'b--','LineWidth',2)
        axis([1 2.0 0 1]);box on;grid minor;
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        subplot(3,2,2)
        plot(yAxis_rescaled,AA0,'g--','LineWidth',2)
        axis([-0.8 0.8 0 1]);box on;grid minor;
        xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        subplot(3,2,4)
        plot(yAxis_rescaled,AA1,'r--','LineWidth',2)
        axis([-0.8 0.8 0 1]);box on;grid minor;
        xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        subplot(3,2,6)
        plot(yAxis_rescaled,AA2,'b--','LineWidth',2)
        axis([-0.8 0.8 0 1]);box on;grid minor;
        xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        
        A = np_c';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        figure(fig);subplot(2,2,2)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        axis([xmin, xmax, ymin, ymax]);
        box on; axis equal;
        ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        
        A = PR_f';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        if ST.params.synthetic_camera_params.photon_count
            B(B<1) = [];
        end
        minval = min(B);
        maxval = threshold*std(B);
        v = linspace(minval,maxval,50);
        
        figure(fig);subplot(2,2,3)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        hold on;plot(xAxis_rescaled,z0*ones(size(xAxis_rescaled)),'r',...
            xAxis_rescaled,z1*ones(size(xAxis_rescaled)),'g',...,
            xAxis_rescaled,z2*ones(size(xAxis_rescaled)),'b','LineWidth',1);hold off
        hold on;plot(r1*ones(size(yAxis_rescaled)),yAxis_rescaled,'r',...
            r0*ones(size(yAxis_rescaled)),yAxis_rescaled,'g',...
            r2*ones(size(yAxis_rescaled)),yAxis_rescaled,'b','LineWidth',1);hold off
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        axis([xmin, xmax, ymin, ymax]);
        box on; axis equal;
        ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        if ST.params.synthetic_camera_params.photon_count
            xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
        else
            xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
        end
        
        C = max(A(I0,:));
        A0 = A(I0,:)/C;A1 = A(I1,:)/C;A2 = A(I2,:)/C;
        CC = max(A(:,J1));
        AA0 = A(:,J0)/CC;AA1 = A(:,J1)/CC;AA2 = A(:,J2)/CC;
        
        figure(fig_slices);
        subplot(3,2,1)
        hold on;plot(xAxis_rescaled,A1,'g-','LineWidth',2);hold off
%         axis([1 2.0 0 1]);box on;grid minor;
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        subplot(3,2,3)
        hold on;plot(xAxis_rescaled,A0,'r-','LineWidth',2);hold off
%         axis([1 2.0 0 1]);box on;grid minor;
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        subplot(3,2,5)
        hold on;plot(xAxis_rescaled,A2,'b-','LineWidth',2);hold off
%         axis([1 2.0 0 1]);box on;grid minor;
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        subplot(3,2,2)
        hold on;plot(yAxis_rescaled,AA0,'g-','LineWidth',2);hold off
%         axis([-0.8 0.8 0 1]);box on;grid minor;
        xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        subplot(3,2,4)
        hold on;plot(yAxis_rescaled,AA1,'r-','LineWidth',2);hold off
%         axis([-0.8 0.8 0 1]);box on;grid minor;
        xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        subplot(3,2,6)
        hold on;plot(yAxis_rescaled,AA2,'b-','LineWidth',2);hold off
%         axis([-0.8 0.8 0 1]);box on;grid minor;
        xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        saveas(fig_slices,[ST.path 'camera_slices_ss_' num2str(ss)],'fig')
        
        A = np_f';
        B = reshape(A,[numel(A),1]);
        B(B==0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        figure(fig);subplot(2,2,4)
        contourf(xAxis_rescaled,yAxis_rescaled,A,v,'LineStyle','none')
        ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
        axis([xmin, xmax, ymin, ymax]);
        box on; axis equal;
        ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
        xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
        
        cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
        ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
        xlabel(hc,'$\rho_{RE}(R,Z)$ (No. particles)','Interpreter','latex','FontSize',12)
        
        saveas(fig,[ST.path 'cameraPicsAnalysis_' num2str(ss)],'fig')
    end
    
end
end

function D = svdAnalysis(A,Eo)
% A has to be a matrix with L2 norm equal to 1

[U,S,V] = svd(A);

diagS = diag(S);
s = sum(diagS.^2);

rankS = numel(diagS);

E = zeros(1,rankS);

for ii=1:rankS
    E(ii) = 100*sum(diagS(1:ii).^2)/s;
end

[~,I] = find(E>Eo,1,'first');

SS = zeros(size(S));
SS(1:I,1:I) = S(1:I,1:I);

D = U*SS*V';

figure;
plot(1:1:rankS,E,'ko--')
xlabel('Mode number','Interpreter','latex')
ylabel('Energy (\%)','Interpreter','latex')
end

function picsAnalysis(ST)
disp('Plotting snapshots...')

plotToroidalSections = true; % <-------------------------------------------
threshold = 4; % <---------------------------------------------------------

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
    
    if (NT ~= 0)
        PR_f = zeros(NX,NY,NT);
        np_f = zeros(NX,NY,NT);
        
        PR_c = zeros(NX,NY,NT);
        np_c = zeros(NX,NY,NT);
    else
        PR_f = zeros(NX,NY);
        np_f = zeros(NX,NY);
        
        PR_c = zeros(NX,NY);
        np_c = zeros(NX,NY);
    end
    
    if (ST.params.synthetic_camera_params.integrated_opt == 0)
        for ii=1:NX
            for jj=1:NY
                if numel(size(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel)) == 4
                    P_cone(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel(ii,jj,i1:i2,:),4);
                    P_full(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel(ii,jj,i1:i2,:),4);
                    
                    np_cone(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel(ii,jj,i1:i2,:),4);
                    np_full(ii,jj,:) = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel(ii,jj,i1:i2,:),4);
                else
                    P_cone(ii,jj,:) = ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel(ii,jj,i1:i2);
                    P_full(ii,jj,:) = ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel(ii,jj,i1:i2);
                    
                    np_cone(ii,jj,:) = ST.data.(['sp' num2str(ss)]).np_lambda_pixel(ii,jj,i1:i2);
                    np_full(ii,jj,:) = ST.data.(['sp' num2str(ss)]).np_angular_pixel(ii,jj,i1:i2);
                end
                
                np_c(ii,jj) = sum(np_cone(ii,jj,:),3);
                np_f(ii,jj) = sum(np_full(ii,jj,:),3);
                
                PR_c(ii,jj) = trapz(lambda(i1:i2),P_cone(ii,jj,:));
                PR_f(ii,jj) = trapz(lambda(i1:i2),P_full(ii,jj,:));
            end
        end
    else
        if (NT ~= 0)
            np_c = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel,4);
            np_f = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel,4);
            
            PR_c = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel,4);
            PR_f = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel,4);
        else
            np_c = sum(ST.data.(['sp' num2str(ss)]).np_lambda_pixel,3);
            np_f = sum(ST.data.(['sp' num2str(ss)]).np_angular_pixel,3);
            
            PR_c = sum(ST.data.(['sp' num2str(ss)]).Psyn_lambda_pixel,3);
            PR_f = sum(ST.data.(['sp' num2str(ss)]).Psyn_angular_pixel,3);
        end
    end
    
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
    
    I0 = find(yAxis_rescaled>0.25,1)-1; z0 = yAxis_rescaled(I0);
    I1 = find(yAxis_rescaled>0,1); z1 = yAxis_rescaled(I1);
    I2 = find(yAxis_rescaled>-0.25,1); z2 = yAxis_rescaled(I2);
    
    J0 = find(xAxis_rescaled>1.3,1)-1; r0 = xAxis_rescaled(J0);
    J1 = find(xAxis_rescaled>1.5,1)-1; r1 = xAxis_rescaled(J1);
    J2 = find(xAxis_rescaled>1.7,1)-1; r2 = xAxis_rescaled(J2);
    
    DX = mean(diff(xAxis_rescaled));
    DY = mean(diff(yAxis_rescaled));
    
    %  % % % % % % % PRE-PROCESSING % % % % % % % % % % % % % % % % % % % %
    
    fig = figure;
    
    if (any(xAxis_rescaled<0))
        xAxis_rescaled = -xAxis_rescaled;
    end
    
    Ac = sum(PR_c,3)';
    Ac = Ac/(sum(sum(Ac))*DX*DY);
    Ac = Ac/sqrt(sum(sum(Ac.^2)));
    Ac = svdAnalysis(Ac,90.0);
    
    B = reshape(Ac,[numel(Ac),1]);
    B(B==0) = [];
    if ST.params.synthetic_camera_params.photon_count
        B(B<1) = [];
    end
    minval = min(B);
    maxval = threshold*std(B);
    v = linspace(minval,maxval,50);
    
    figure(fig);subplot(1,2,1)
    contourf(xAxis_rescaled,yAxis_rescaled,Ac,v,'LineStyle','none')
    ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
    axis([xmin, xmax, ymin, ymax]);
    box on; axis equal;grid minor
    ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
    xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
    
    cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
    if ST.params.synthetic_camera_params.photon_count
        xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
    else
        xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
    end  
    
    Af = sum(PR_f,3)';
    Af = Af/(sum(sum(Af))*DX*DY);
    Af = Af/sqrt(sum(sum(Af.^2)));
    Af = svdAnalysis(Af,90.0);
    
    B = reshape(Af,[numel(Af),1]);
    B(B==0) = [];
    if ST.params.synthetic_camera_params.photon_count
        B(B<1) = [];
    end
    minval = min(B);
    maxval = threshold*std(B);
    v = linspace(minval,maxval,50);
    
    figure(fig);subplot(1,2,2)
    contourf(xAxis_rescaled,yAxis_rescaled,Af,v,'LineStyle','none')
    ymin=min(yAxis_rescaled);ymax=max(yAxis_rescaled);xmin=min(xAxis_rescaled);xmax=max(xAxis_rescaled);
    axis([xmin, xmax, ymin, ymax]);
    box on; axis equal;grid minor
    ylabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
    xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
    
    cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);hc = colorbar('Location','eastoutside');caxis([0,maxval]);
    ax = gca;ax.Color = [1,1,1];ax.ClippingStyle = 'rectangle';
    if ST.params.synthetic_camera_params.photon_count
        xlabel(hc,'$P_R$ (Photons)','Interpreter','latex','FontSize',12)
    else
        xlabel(hc,'$\int P_R(\lambda) d\lambda$ (Watts)','Interpreter','latex','FontSize',12)
    end
    
    
    hs0c = Ac(I0,:);hs1c = Ac(I1,:);hs2c = Ac(I2,:);
    vs0c = Ac(:,J0);vs1c = Ac(:,J1);vs2c = Ac(:,J2);
    
    hs0f = Af(I0,:);hs1f = Af(I1,:);hs2f = Af(I2,:);
    vs0f = Af(:,J0);vs1f = Af(:,J1);vs2f = Af(:,J2);
    
    ymax = max([max(max(Ac)) max(max(Af))]);
    
    fig_slices=figure;
    subplot(3,2,1)
    plot(xAxis_rescaled,hs0c,'g--',xAxis_rescaled,hs0f,'g-','LineWidth',2)
    axis([1 2.0 0 ymax]);box on;grid minor;
    xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
    subplot(3,2,3)
    plot(xAxis_rescaled,hs1c,'r--',xAxis_rescaled,hs1f,'r-','LineWidth',2)
    axis([1 2.0 0 ymax]);box on;grid minor;
    xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
    subplot(3,2,5)
    plot(xAxis_rescaled,hs2c,'b--',xAxis_rescaled,hs2f,'b-','LineWidth',2)
    axis([1 2.0 0 ymax]);box on;grid minor;
    xlabel('$R$ (m)','FontSize',12,'Interpreter','latex')
    
    subplot(3,2,2)
    plot(yAxis_rescaled,vs0c,'g--',yAxis_rescaled,vs0f,'g-','LineWidth',2)
    axis([-0.8 0.8 0 ymax]);box on;grid minor;
    xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
    subplot(3,2,4)
    plot(yAxis_rescaled,vs1c,'r--',yAxis_rescaled,vs1f,'r-','LineWidth',2)
    axis([-0.8 0.8 0 ymax]);box on;grid minor;
    xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
    subplot(3,2,6)
    plot(yAxis_rescaled,vs2c,'b--',yAxis_rescaled,vs2f,'b-','LineWidth',2)
    axis([-0.8 0.8 0 ymax]);box on;grid minor;
    xlabel('$Z$ (m)','FontSize',12,'Interpreter','latex')
    
    saveas(fig_slices,[ST.path 'slices_ss_' num2str(ss)],'fig')
    
    for ii=1:2
        figure(fig);subplot(1,2,ii);
        hold on;plot(xAxis_rescaled,z0*ones(size(xAxis_rescaled)),'r',...
            xAxis_rescaled,z1*ones(size(xAxis_rescaled)),'g',...,
            xAxis_rescaled,z2*ones(size(xAxis_rescaled)),'b','LineWidth',1);hold off
        hold on;plot(r1*ones(size(yAxis_rescaled)),yAxis_rescaled,'r',...
            r0*ones(size(yAxis_rescaled)),yAxis_rescaled,'g',...
            r2*ones(size(yAxis_rescaled)),yAxis_rescaled,'b','LineWidth',1);hold off
    end
    
    
    saveas(fig,[ST.path 'cameraPicsAnalysis_' num2str(ss)],'fig')
end
end