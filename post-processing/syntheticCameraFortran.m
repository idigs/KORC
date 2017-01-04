function ST = syntheticCameraFortran(path,lambdas,visible,range)
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

plotSyntheticCameraAnalysis(ST)

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
end

function data = loadData(ST)
data = struct;

list = {'Psyn_pixel','part_pixel'};

it = ST.range(1):1:ST.range(2);

NX = ST.params.synthetic_camera_params.num_pixels(1);
NY = ST.params.synthetic_camera_params.num_pixels(2);
Nl = ST.params.synthetic_camera_params.Nlambda;

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    for ss=1:ST.params.simulation.num_species
        data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,Nl,ST.num_snapshots);
        filename = [ST.path 'synthetic_camera_snapshots.h5'];
        for ii=1:numel(it)
            dataset = ...
                ['/' num2str(it(ii)*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                '/' list{ll}];
            data.(['sp' num2str(ss)]).(list{ll})(:,:,:,ii) = h5read(filename, dataset);
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
% Nl = ST.params.synthetic_camera_params.Nlambda;

[~,i1] = min(abs(lambda - ST.lambdas(1)));
[~,i2] = min(abs(lambda - ST.lambdas(2)));
Nl = i2 - i1 + 1;

for ss=1:ST.params.simulation.num_species
    disp(['Species: ' num2str(ss)])
    Psyn_lambda = zeros(NX,NY,Nl);
    Npart_lambda = zeros(NX,NY,Nl);
    Psyn = zeros(NX,NY);
    Npart = zeros(NX,NY);
    
    for it=1:ST.num_snapshots
        for ii=1:NX
            for jj=1:NY
                Psyn_lambda(ii,jj,:) = Psyn_lambda(ii,jj,:) + ...
                    ST.data.(['sp' num2str(ss)]).Psyn_pixel(ii,jj,i1:i2,it);
                Npart_lambda(ii,jj,:) = Npart_lambda(ii,jj,:) + ...
                    ST.data.(['sp' num2str(ss)]).part_pixel(ii,jj,i1:i2,it);
                
                Psyn(ii,jj) = Psyn(ii,jj) + trapz(lambda(i1:i2),ST.data.(['sp' num2str(ss)]).Psyn_pixel(ii,jj,i1:i2,it));
                Npart(ii,jj) = Npart(ii,jj) + sum(ST.data.(['sp' num2str(ss)]).part_pixel(ii,jj,i1:i2,it));
            end
        end
    end
    
    axis_lambda = 1E9*lambda(i1:i2);
    Psyn_lambda = 1E-9*Psyn_lambda;
%     Psyn_mean = squeeze(mean(mean(Psyn_lambda,1),2));
    
    h = figure;
    subplot(3,2,[1 3])
    surfc(xAxis,yAxis,Psyn','LineStyle','none')
    colormap(jet); hc = colorbar('Location','southoutside');
    xlabel(hc,'$P_{syn}$ (W/sr)','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$y$-axis','FontSize',12,'Interpreter','latex')
    xlabel('$x$-axis','FontSize',12,'Interpreter','latex')
    
    figure(h);
    subplot(3,2,5)
    hold on
    for ii=1:NX
        for jj=1:NY
            plot(axis_lambda,squeeze(Psyn_lambda(ii,jj,:)))
        end
    end
%     plot(lambda,Psyn_mean,'k','LineWidth',3)
    hold off
    box on;
    ylabel('$P_{syn}$ (W/(nm$\cdot$sr))','FontSize',12,'Interpreter','latex')
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
    
    figure(h);
    subplot(3,2,[2 4])
    surfc(xAxis,yAxis,Npart','LineStyle','none')
    colormap(jet);  hc = colorbar('Location','southoutside');
    xlabel(hc,'Number of RE','Interpreter','latex','FontSize',12)
    box on; axis square;view([0 -90])
    ylabel('$y$-axis','FontSize',14,'Interpreter','latex')
    xlabel('$x$-axis','FontSize',14,'Interpreter','latex')
    
    figure(h);
    subplot(3,2,6)
    hold on
    for ii=1:NX
        for jj=1:NY
            plot(axis_lambda,squeeze(Npart_lambda(ii,jj,:)))
        end
    end
    hold off
    box on;
    ylabel('Number of RE','FontSize',12,'Interpreter','latex')
    xlabel('$\lambda$ (nm)','FontSize',12,'Interpreter','latex')
    
    saveas(h,[ST.path 'SyntheticCameraFortran_ss_' num2str(ss)],'fig')
end
end

