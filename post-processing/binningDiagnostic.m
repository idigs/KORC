function ST = binningDiagnostic(path,range)
% ST = binningDiagnostic('../KORC-FO/outputFiles/',[99,100])
close all

ST = struct;
ST.path = path;

ST.range = range;
ST.num_snapshots = ST.range(2) - ST.range(1) + 1;

ST.params = loadSimulationParameters(ST);

ST.time = ...
    ST.params.simulation.dt*double(ST.params.simulation.output_cadence)*double(ST.range(1):1:ST.range(2));

ST.data = loadData(ST);

plotBinningSnapshots(ST);

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

info = h5info([ST.path 'binning_diagnostic.h5']);

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

filename = [ST.path 'binning_diagnostic_snapshots.h5'];
H = h5info(filename);
it = zeros(1,size(H.Groups,1));
for ii=1:size(H.Groups,1)
    tmpstr = strrep(H.Groups(ii).Name,'/','');
    it(ii) = str2num(tmpstr);
end
it = sort(it,'ascend');

NX = ST.params.binning_diagnostic_params.num_bins(1);
NY = ST.params.binning_diagnostic_params.num_bins(2);

if (ST.params.binning_diagnostic_params.toroidal_sections == 1)
    Ntor = ST.params.binning_diagnostic_params.ntor_sections;
else
    Ntor = 0;
end

filename = [ST.path 'binning_diagnostic_snapshots.h5'];


list = {'eta','g','N'};

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    for ss=1:ST.params.simulation.num_species
        if (ST.params.binning_diagnostic_params.toroidal_sections == 1)
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,Ntor,ST.num_snapshots);
        else
            data.(['sp' num2str(ss)]).(list{ll}) = zeros(NX,NY,ST.num_snapshots);
        end
        for ii=ST.range(1)+1:ST.range(2)+1
            dataset = ...
                ['/' num2str(it(ii)) '/spp_' num2str(ss)...
                '/' list{ll}];
            if (ST.params.binning_diagnostic_params.toroidal_sections == 1)
                data.(['sp' num2str(ss)]).(list{ll})(:,:,:,ii) = h5read(filename, dataset);
            else
                data.(['sp' num2str(ss)]).(list{ll})(:,:,ii) = h5read(filename, dataset);
            end
            
        end
    end
end


end

function plotBinningSnapshots(ST)
disp('Plotting snapshots...')
plotToroidalSections = false;
figuresToShare = false;

xAxis = ST.params.binning_diagnostic_params.rnodes;
yAxis = ST.params.binning_diagnostic_params.znodes;

NX = ST.params.binning_diagnostic_params.num_bins(1);
NY = ST.params.binning_diagnostic_params.num_bins(2);

if (ST.params.binning_diagnostic_params.toroidal_sections == 1)
    NT = ST.params.binning_diagnostic_params.ntor_sections;
else
    NT = 0;
end

for ss=1:ST.params.simulation.num_species
    num_dims = ndims(ST.data.sp1.eta);
    
    if num_dims > 2
        eta = sum(ST.data.(['sp' num2str(ss)]).eta,num_dims);
        g = sum(ST.data.(['sp' num2str(ss)]).g,num_dims);
        N = sum(ST.data.(['sp' num2str(ss)]).N,num_dims);
    else
        eta = ST.data.(['sp' num2str(ss)]).eta;
        g = ST.data.(['sp' num2str(ss)]).g;
        N = ST.data.(['sp' num2str(ss)]).N;
    end
    
    eta(N~=0) = eta(N~=0)./N(N~=0);
    g(N~=0) = g(N~=0)./N(N~=0);
    
    if (ST.params.binning_diagnostic_params.toroidal_sections && ~plotToroidalSections)
        eta = squeeze(sum(eta,3));
        g = squeeze(sum(g,3));
    end
    
    fig1 = figure;
    subplot(1,2,1)
    surfc(xAxis,yAxis,eta,'LineStyle','none')
    axis equal;view([0 90])
    colormap(jet);cb = colorbar;
    xlabel('$R$','Interpreter','latex')
    ylabel('$Z$','Interpreter','latex')
    ylabel(cb,'Mean pitch angle ($^\circ$)','Interpreter','latex')
    
    subplot(1,2,2)
    surfc(xAxis,yAxis,g,'LineStyle','none')
    axis equal;view([0 90])
    colormap(jet);cb = colorbar;
    xlabel('$R$','Interpreter','latex')
    ylabel('$Z$','Interpreter','latex')
    ylabel(cb,'Mean energy (MeV)','Interpreter','latex')
    
end
end