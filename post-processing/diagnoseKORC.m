function ST = diagnoseKORC(path)
close all

ST = struct;
ST.path = path;

ST.params = loadSimulationParameters(ST);

ST.data = loadData(ST);

% energyConservation(ST);

ST.PD = pitchAngleDiagnostic(ST,30);

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

end

function data = loadData(ST)
data = struct;

list = {'X'};%,'V','Rgc'};

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    for ss=1:ST.params.simulation.num_species
        tnp = double(ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
        
        data.(['sp' num2str(ss)]).(list{ll}) = ...
            zeros(3,tnp,ST.params.simulation.num_snapshots);
        
        for ff=1:ST.params.simulation.nmpi
            filename = [ST.path 'file_' num2str(ff-1) '.h5'];
            indi = (ff - 1)*double(ST.params.species.ppp(ss)) + 1;
            indf = ff*double(ST.params.species.ppp(ss));
            for ii=1:ST.params.simulation.num_snapshots
                dataset = ...
                    ['/' num2str(ii*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                
                data.(['sp' num2str(ss)]).(list{ll})(:,indi:indf,ii) = ...
                    h5read(filename, dataset);
            end
            
        end
        
    end
end


list = {'eta','gamma'};%,'mu','kappa','tau'};

for ll=1:length(list)
    for ss=1:ST.params.simulation.num_species
        tnp = double(ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
        
        data.(['sp' num2str(ss)]).(list{ll}) = ...
            zeros(tnp,ST.params.simulation.num_snapshots);
        
        for ii=1:ST.params.simulation.num_snapshots
            disp(['Loading: ' list{ll} ' Snapshot: ' num2str(ii)])
            
            for ff=1:ST.params.simulation.nmpi
                
                filename = [ST.path 'file_' num2str(ff-1) '.h5'];
                indi = (ff - 1)*double(ST.params.species.ppp(ss)) + 1;
                indf = ff*double(ST.params.species.ppp(ss));
                
                dataset = ...
                    ['/' num2str(ii*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                
                data.(['sp' num2str(ss)]).(list{ll})(indi:indf,ii) = ...
                    h5read(filename, dataset);
            end
            
        end
        
    end
end

end

function energyConservation(ST)

try
    
catch
    
end

end

function PD = pitchAngleDiagnostic(ST,numBins)
PD = struct;

tmp = [];

for ss=1:ST.params.simulation.num_species
    tmp = [tmp;ST.data.(['sp' num2str(ss)]).eta];
end

minVal = min(min( tmp ));
maxVal = max(max( tmp ));
vals = linspace(minVal,maxVal,numBins);

f = zeros(numBins,ST.params.simulation.num_snapshots);


for ii=1:ST.params.simulation.num_snapshots
    [f(:,ii),~] = hist(tmp(:,ii),vals);
end

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(cad:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

figure
surf(time,vals,log10(f),'LineStyle','none')
axis([tmin tmax minVal maxVal])
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('Pitch angle $\theta$ (degrees)','Interpreter','latex','FontSize',16)
colormap(jet)




end