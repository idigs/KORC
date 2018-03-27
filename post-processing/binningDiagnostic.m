function ST = binningDiagnostic(path,range)
% ST = binningDiagnostic('../KORC-FO/outputFiles/',[99,100])
close all

ST = struct;
ST.path = path;

ST.range = range;
% ST.num_snapshots = ST.range(2) - ST.range(1) + 1;

ST.params = loadSimulationParameters(ST);

% ST.time = ...
%     ST.params.simulation.dt*double(ST.params.simulation.output_cadence)*double(ST.range(1):1:ST.range(2));

ST.data = loadData(ST);

plotBinningSnapshots(ST);

plotRadialDensity(ST);

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

cadence = double(ST.params.simulation.output_cadence);
[~,I1] = min(abs(cadence*ST.range(1)-it));
[~,I2] = min(abs(cadence*ST.range(2)-it));

indices = it(I1:I2);

ST.num_snapshots = I2 - I1 + 1;

% ST.params = loadSimulationParameters(ST);

ST.time = zeros(1,ST.num_snapshots);

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
        for ii=1:ST.num_snapshots
            dataset = ...
                ['/' num2str(indices(ii)) '/spp_' num2str(ss)...
                '/' list{ll}];
            if (ST.params.binning_diagnostic_params.toroidal_sections == 1)
                data.(['sp' num2str(ss)]).(list{ll})(:,:,:,ii) = h5read(filename, dataset);
            else
                data.(['sp' num2str(ss)]).(list{ll})(:,:,ii) = h5read(filename, dataset);
            end
            
        end
    end
end

% try
%     list = ST.params.simulation.outputs_list;
%     
%     for ll=1:length(list)
%         disp(['Loading ' list{ll}])
%         for ss=1:ST.params.simulation.num_species
%             tnp = double(ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
%             
%             if (strcmp(list{ll},'X') || strcmp(list{ll},'V') || strcmp(list{ll},'B') || strcmp(list{ll},'E'))
%                 data.(['sp' num2str(ss)]).('raw').(list{ll}) = zeros(3,tnp,ST.num_snapshots);
%             else
%                 data.(['sp' num2str(ss)]).('raw').(list{ll}) = zeros(tnp,ST.num_snapshots);
%             end
%             
%             for ff=1:ST.params.simulation.nmpi
%                 filename = [ST.path 'file_' num2str(ff-1) '.h5'];
%                 indi = (ff - 1)*double(ST.params.species.ppp(ss)) + 1;
%                 indf = ff*double(ST.params.species.ppp(ss));
%                 for ii=1:ST.num_snapshots
%                     dataset = ...
%                         ['/' num2str(indices(ii)) '/spp_' num2str(ss)...
%                         '/' list{ll}];
%                     if (strcmp(list{ll},'X') || strcmp(list{ll},'V') || strcmp(list{ll},'B') || strcmp(list{ll},'E'))
%                         data.(['sp' num2str(ss)]).('raw').(list{ll})(:,indi:indf,ii) = ...
%                             h5read(filename, dataset);
%                     else
%                         data.(['sp' num2str(ss)]).('raw').(list{ll})(indi:indf,ii) = ...
%                             h5read(filename, dataset);
%                     end
%                 end
%             end
%         end
%     end
% catch
% 	disp('No RAW data found...')
% end

end

function plotBinningSnapshots(ST)
disp('Plotting snapshots...')
plotToroidalSections = false;
figuresToShare = false;

c = ST.params.scales.v;

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
    m = ST.params.species.m(ss);
    q = abs(ST.params.species.q(ss));
    
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
    
    E = (g-1)*m*c^2/(1E6*q);
    
    cmax = max(max(eta));
    cmin = min(min(eta));
    
    fig1 = figure;
    subplot(1,3,1)
    surfc(xAxis,yAxis,eta','LineStyle','none')
    axis equal;view([0 90])
    colormap(jet);cb = colorbar;caxis([cmin cmax])
    xlabel('$R$','Interpreter','latex')
    ylabel('$Z$','Interpreter','latex')
    ylabel(cb,'Mean pitch angle ($^\circ$)','Interpreter','latex')
    
    cmax = max(max(E(g~=0)));
    cmin = min(min(E(g~=0)));
    
    subplot(1,3,2)
    surfc(xAxis,yAxis,E','LineStyle','none')
    axis equal;view([0 90])
    colormap(jet);cb = colorbar;caxis([cmin cmax])
    xlabel('$R$','Interpreter','latex')
    ylabel('$Z$','Interpreter','latex')
    ylabel(cb,'Mean energy (MeV)','Interpreter','latex')
    
    
    cmin = min(min(N));
    cmax = max(max(N));
    
    subplot(1,3,3)
    surfc(xAxis,yAxis,N','LineStyle','none')
    axis equal;view([0 90])
    colormap(jet);cb = colorbar;caxis([cmin cmax])
    xlabel('$R$','Interpreter','latex')
    ylabel('$Z$','Interpreter','latex')
    ylabel(cb,'No. particles','Interpreter','latex')
    
end
end

function plotRadialDensity(ST)
disp('Plotting radial density...')
Nq = 20;
angle = [30 45 60];
falloff_rate = [1 2 3 4];

xAxis = ST.params.binning_diagnostic_params.rnodes;
yAxis = ST.params.binning_diagnostic_params.znodes;

NX = ST.params.binning_diagnostic_params.num_bins(1);
NY = ST.params.binning_diagnostic_params.num_bins(2);

if (ST.params.binning_diagnostic_params.toroidal_sections == 1)
    NT = ST.params.binning_diagnostic_params.ntor_sections;
else
    NT = 0;
end

Ro = ST.params.fields.Ro;
Zo = ST.params.fields.Zo;

fig = figure;
nrows = numel(angle);
for ss=1:ST.params.simulation.num_species
    num_dims = ndims(ST.data.sp1.eta);
    
    if num_dims > 2
        N = sum(ST.data.(['sp' num2str(ss)]).N,num_dims);
    else
        N = ST.data.(['sp' num2str(ss)]).N;
    end
    
    cmin = min(min(N));
    cmax = max(max(N));
    
    N = N';

    for aa=1:nrows
        Rq1 = linspace(Ro,max(xAxis),Nq);
        Zq1 = tan(deg2rad(angle(aa)))*(Rq1 - Ro) + Zo;
        rq1 = 100*sqrt((Rq1 - Ro).^2 + (Zq1 - Zo).^2);
        
        Rq2 = linspace(min(xAxis),Ro,Nq);
        Zq2 = tan(deg2rad(angle(aa)))*(Rq2 - Ro) + Zo;
        rq2 = 100*sqrt((Rq2 - Ro).^2 + (Zq2 - Zo).^2);
        
        
        subplot(nrows,2,(aa-1)*2 + 1)
        contourf(xAxis,yAxis,N,10,'LineStyle','none')
        axis equal;view([0 90]);axis([1 2.5 -1 1])
        colormap(jet);cb = colorbar;caxis([cmin cmax])
        xlabel('$R$','Interpreter','latex')
        ylabel('$Z$','Interpreter','latex')
        ylabel(cb,'No. particles','Interpreter','latex')
        
        [R,Z] = meshgrid(xAxis,yAxis);
        
        Nq1 = interp2(R,Z,N,Rq1,Zq1);
        Nq2 = interp2(R,Z,N,Rq2,Zq2);
        
        figure(fig)
        subplot(nrows,2,(aa-1)*2 + 1)
        hold on;
        plot(Rq1(Nq1~=0),Zq1(Nq1~=0),'go',Rq2(Nq2~=0),Zq2(Nq2~=0),'mo','MarkerFaceColor',[0.6,0.6,0.6],'MarkerSize',4)
        hold off
               
        figure(fig)
        subplot(nrows,2,(aa-1)*2 + 2)
        plot(rq1,log10(Nq1),'go-',rq2,log10(Nq2),'mo-','MarkerFaceColor',[0,0,0])
        grid minor; box on;xlim([0 60])
        xlabel('$r$ (cm)','Interpreter','latex')
        ylabel('$N$','Interpreter','latex')
        
        for ff=1:numel(falloff_rate)
            r = linspace(0,max([max(rq1(Nq1~=0)) max(rq2(Nq2~=0))]),Nq);
            log10No = max([max(log10(Nq1(Nq1~=0))) max(log10(Nq2(Nq2~=0)))]);
            log10Nq = -0.01*falloff_rate(ff)*r + log10No;
            
            figure(fig)
            subplot(nrows,2,(aa-1)*2 + 2)
            hold on
            plot(r,log10Nq,'k--')
            hold off
            text(1.05*max(r),min(log10Nq),num2str(falloff_rate(ff)),'Interpreter','latex')
        end
    end
end
end