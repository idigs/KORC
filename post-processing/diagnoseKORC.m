function ST = diagnoseKORC(path)
close all

ST = struct;
ST.path = path;

ST.params = loadSimulationParameters(ST);

ST.data = loadData(ST);

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

dataCell = {};

for ss=1:ST.params.simulation.num_species
    for ff=1:ST.params.simulation.nmpi
        
    end
    
end
end