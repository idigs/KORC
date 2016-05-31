function ST = diagnoseKORC(path)
close all

ST = struct;
ST.path = path;

ST.params = loadSimulationParameters(ST);

ST.data = loadData(ST);

energyConservation(ST);

ST.PD = pitchAngleDiagnostic(ST,50);

ST.PD = magneticMomentDiagnostic(ST,50);

sp = 'sp1'

R = ST.params.scales.l*squeeze( sqrt( ST.data.(sp).X(1,:,:).^2 + ST.data.(sp).X(2,:,:).^2 ) );

[V,I] = max(R');
[~,II] = max(V);

% [V,I] = min(R');
% [~,II] = min(V);

X = squeeze(ST.data.(sp).X(:,II,:))*ST.params.scales.l;

figure
plot3(X(1,:),X(2,:),X(3,:))
axis equal; box on

R = sqrt( X(1,:).^2 + X(2,:).^2 );
Z = X(3,:);
figure;plot(R,Z);axis equal

X = squeeze(ST.data.(sp).Rgc(:,II,:))*ST.params.scales.l;
R = sqrt( X(1,:).^2 + X(2,:).^2 );
Z = X(3,:);
hold on;plot(R,Z);hold off
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

% params.simulation.num_snapshots = 3828;
% params.simulation.t_steps = 19140000;

end

function data = loadData(ST)
data = struct;

list = {'X','Rgc'};%,'V','Rgc'};

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
                    ['/' num2str((ii-1)*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                
                data.(['sp' num2str(ss)]).(list{ll})(:,indi:indf,ii) = ...
                    h5read(filename, dataset);
            end
            
        end
        
    end
end


% list = {'eta','gamma'};%,'mu','kappa','tau'};
list = {'eta','gamma','kappa','mu'};

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    for ss=1:ST.params.simulation.num_species
        tnp = double(ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
        
        data.(['sp' num2str(ss)]).(list{ll}) = ...
            zeros(tnp,ST.params.simulation.num_snapshots);
        
        for ii=1:ST.params.simulation.num_snapshots
            % disp(['Loading: ' list{ll} ' Snapshot: ' num2str((ii-1))])
            
            for ff=1:ST.params.simulation.nmpi
                
                filename = [ST.path 'file_' num2str(ff-1) '.h5'];
                indi = (ff - 1)*double(ST.params.species.ppp(ss)) + 1;
                indf = ff*double(ST.params.species.ppp(ss));
                
                dataset = ...
                    ['/' num2str((ii-1)*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                
                data.(['sp' num2str(ss)]).(list{ll})(indi:indf,ii) = ...
                    h5read(filename, dataset);
            end
            
        end
        
    end
end

end

function energyConservation(ST)

err = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(cad:cad:ST.params.simulation.t_steps);

h=figure;

try
    for ss=1:ST.params.simulation.num_species
        tmp = zeros(size(ST.data.(['sp' num2str(ss)]).gamma));
        for ii=1:ST.params.species.ppp(ss)
            tmp(ii,:) = ...
                100*( ST.data.(['sp' num2str(ss)]).gamma(ii,1) - ...
                ST.data.(['sp' num2str(ss)]).gamma(ii,:) )./ST.data.(['sp' num2str(ss)]).gamma(ii,1);
        end
        err(:,ss) = mean(tmp,1);
        figure(h)
        hold on
        plot(time,err(:,ss))
        hold off
    end
catch
    error('Something went wrong: energyConservation')
end

figure(h)
box on
grid on
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('Energy conservation (\%)','Interpreter','latex','FontSize',16)

end

function PD = pitchAngleDiagnostic(ST,numBins)
PD = struct;

mean_pitch = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
std_pitch = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
skewness_pitch = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);

f_tot = zeros(ST.params.simulation.num_species,numBins);

tmp = [];

for ss=1:ST.params.simulation.num_species
    mean_pitch(ss,:) = mean(ST.data.(['sp' num2str(ss)]).eta,1);
    std_pitch(ss,:) = std(ST.data.(['sp' num2str(ss)]).eta,0,1);
    skewness_pitch(ss,:) = skewness(ST.data.(['sp' num2str(ss)]).eta,0,1);
    
    tmp = [tmp; ST.data.(['sp' num2str(ss)]).eta];
    
    aux = reshape(ST.data.(['sp' num2str(ss)]).eta,1,numel(ST.data.(['sp' num2str(ss)]).eta));
    
    minVal = min(min( aux ));
    maxVal = max(max( aux ));
    vals = linspace(minVal,maxVal,numBins);
    
    f_tot(ss,:) =  hist(aux,vals);
    figure
    plot(vals,f_tot(ss,:))
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
% surf(time,vals,f,'LineStyle','none')
axis([tmin tmax minVal maxVal])
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('Pitch angle $\theta$ (degrees)','Interpreter','latex','FontSize',16)
colormap(jet)

figure
for ii=1:ST.params.simulation.num_species
    subplot(3,1,1)
    hold on
    plot(time,mean_pitch(ii,:))
    hold off
    subplot(3,1,2)
    hold on
    plot(time,std_pitch(ii,:))
    hold off
    subplot(3,1,3)
    hold on
    plot(time,skewness_pitch(ii,:))
    hold off
end

PD.f = f;
PD.vals = vals;
end

function PD = magneticMomentDiagnostic(ST,numBins)
PD = struct;

tmp = [];

for ss=1:ST.params.simulation.num_species
    tmp = [tmp;ST.data.(['sp' num2str(ss)]).mu];
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
% surf(time,vals,f,'LineStyle','none')
axis([tmin tmax minVal maxVal])
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('Magnetic moment $\mu$ (arbitrary units)','Interpreter','latex','FontSize',16)
colormap(jet)

PD.f = f;
PD.vals = vals;
end