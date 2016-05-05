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

for ss=1:ST.params.simulation.num_species
    data.(['sp' num2str(ss)]) = struct;
    for ff=1:ST.params.simulation.nmpi 
        filename = [ST.path 'file_' num2str(ff-1) '.h5'];
        for ii=1:ST.params.simulation.num_snapshots
            position = ...
                ['/' num2str(ii*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)];
            h5disp(filename,position,'min')
            % y = h5read(filename)
        end
    end
end

end