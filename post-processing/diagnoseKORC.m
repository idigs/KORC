function ST = diagnoseKORC(path)
close all

ST = struct;
ST.path = path;

ST.params = loadSimulationParameters(ST);

ST.data = loadData(ST);

energyConservation(ST);

pitchAngleDiagnostic(ST,100);

magneticMomentDiagnostic(ST,150);



% % % % % % % % % % % % % 
sp = 'sp1'
R = squeeze( sqrt( ST.data.(sp).X(1,:,:).^2 + ST.data.(sp).X(2,:,:).^2 ) );

[V,I] = max(R');
[~,II] = max(V);
% [V,I] = min(R');
% [~,II] = min(V);

X = squeeze(ST.data.(sp).X(:,II,:));

figure
plot3(X(1,:),X(2,:),X(3,:))
axis equal; box on

R = sqrt( X(1,:).^2 + X(2,:).^2 );
Z = X(3,:);
figure;plot(R,Z);axis equal

X = squeeze(ST.data.(sp).Rgc(:,II,:));
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


% list = {'eta','gamma'};%,'mu','Prad','tau'};
list = {'eta','gamma','mu'};

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
set(h,'name','Energy conservation','numbertitle','off')
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

function pitchAngleDiagnostic(ST,numBins)
N = 10;
nbins = 30;

mean_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
std_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
skewness_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
kurtosis_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);

fx = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots);
x = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots);

f_tot = zeros(numBins,ST.params.simulation.num_snapshots);

data = [];

for ss=1:ST.params.simulation.num_species
    tmp = cos(pi*ST.data.(['sp' num2str(ss)]).eta/180);
%     tmp = ST.data.(['sp' num2str(ss)]).eta;
    mean_f(ss,:) = mean(tmp,1);
    std_f(ss,:) = std(tmp,0,1);
    skewness_f(ss,:) = skewness(tmp,0,1);
    kurtosis_f(ss,:) = kurtosis(tmp,1,1);
    
    data = [data; ST.data.(['sp' num2str(ss)]).eta];
    
    for ii=1:ST.params.simulation.num_snapshots
        [fx(ss,:,ii),x(ss,:,ii)] = ...
            hist(tmp(:,ii),nbins);
        
        dx = mean(diff(x(ss,:,ii)));
        fx(ss,:,ii) = fx(ss,:,ii)/(sum(fx(ss,:,ii))*dx);
        
        x(ss,:,ii) = ( x(ss,:,ii) - mean_f(ss,ii) )/std_f(ss,ii);
        fx(ss,:,ii) = std_f(ss,ii)*fx(ss,:,ii);
    end
end

minVal = min(min( data ));
maxVal = max(max( data ));
vals = linspace(minVal,maxVal,numBins);

for ii=1:ST.params.simulation.num_snapshots
    [f_tot(:,ii),~] = hist(data(:,ii),vals);
end

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(cad:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

h1 = figure;
set(h1,'name','PDF time evolution','numbertitle','off')
surf(time,vals,log10(f_tot),'LineStyle','none')
% surf(time,vals,f,'LineStyle','none')
axis([tmin tmax minVal maxVal])
box on
axis on
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('Pitch angle $\theta$ (degrees)','Interpreter','latex','FontSize',16)
colormap(jet)

h2 = figure
set(h2,'name','Statistical moments','numbertitle','off')
for ii=1:ST.params.simulation.num_species
    figure(h2)
    subplot(4,1,1)
    hold on
    plot(time,mean_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,2)
    hold on
    plot(time,std_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,3)
    hold on
    plot(time,skewness_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,4)
    hold on
    plot(time,kurtosis_f(ii,:))
    hold off
end

figure(h2)
subplot(4,1,1)
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('mean($\theta$)','Interpreter','latex','FontSize',16)
box on
grid on
figure(h2)
subplot(4,1,2)
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('std($\theta$)','Interpreter','latex','FontSize',16)
box on
grid on
figure(h2)
subplot(4,1,3)
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('skewness($\theta$)','Interpreter','latex','FontSize',16)
box on
grid on
figure(h2)
subplot(4,1,4)
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('kurtosis($\theta$)','Interpreter','latex','FontSize',16)
box on
grid on


offset = floor(double(ST.params.simulation.num_snapshots)/N);

z = linspace(-4,4,100);
fz = exp( -0.5*z.^2 )/sqrt(2*pi);

h3 = figure;
set(h3,'name','PDF pitch angle','numbertitle','off')
for ii=1:N
    it = ii*offset;
%     subplot(N,1,ii)
    nc = floor(N/2);
    nr = floor(N/nc);
    subplot(nr,nc,ii)
    plot(z,log10(fz),'k')
    hold on
    for ss=1:ST.params.simulation.num_species
        figure(h3)
        plot(x(ss,:,it),log10(fx(ss,:,it)),'o:')
        
        figure(h2)
        subplot(4,1,1)
        hold on
        plot(time(it),mean_f(ss,it),'rs')
        hold off
        figure(h2)
        subplot(4,1,2)
        hold on
        plot(time(it),std_f(ss,it),'rs')
        hold off
        figure(h2)
        subplot(4,1,3)
        hold on
        plot(time(it),skewness_f(ss,it),'rs')
        hold off
        figure(h2)
        subplot(4,1,4)
        hold on
        plot(time(it),kurtosis_f(ss,it),'rs')
        hold off
    end
    figure(h3)
    title(['Time: ' num2str(time(it))],'Interpreter','latex','FontSize',11)
    hold off
    box on
    grid on
end


ft = zeros(ST.params.simulation.num_species,nbins,N);
t = zeros(ST.params.simulation.num_species,nbins,N);

for ii=1:N
    it = ii*offset;
    for ss=1:ST.params.simulation.num_species
        I = find( any(ST.data.sp3.eta > 90, 2) == 0 );
        X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,I,it));
        theta = atan2(X(3,:), sqrt(X(1,:).^2 + X(2,:).^2) - ST.params.fields.Ro);
        theta(theta < 0) = theta(theta < 0) + 2*pi;
        theta = (180/pi)*theta;
        [ft(ss,:,ii),t(ss,:,ii)] = hist(theta,nbins);
        dt = mean(diff(t(ss,:,ii)));
        ft(ss,:,ii) = ft(ss,:,ii)/(sum(ft(ss,:,ii))*dt);
    end
end

barcolor = [1,0,0;0,1,0;0,0,1];

h4 = figure;
set(h4,'name','PDF poloidal angle','numbertitle','off')
for ii=1:N
    it = ii*offset;
    nc = floor(N/2);
    nr = floor(N/nc);
    figure(h4)
    subplot(nr,nc,ii)
    hold on
    for ss=1:ST.params.simulation.num_species
%         plot(t(ss,:,ii),ft(ss,:,ii),'o:')
        bar(t(ss,:,ii),ft(ss,:,ii),'FaceColor',barcolor(ss,:))
    end
    currentAxis = gca;
    set(currentAxis.Children(:),'FaceAlpha',0.2);
    hold off
    xlim([0 360])
    title(['Time: ' num2str(time(it))],'Interpreter','latex','FontSize',11)
    xlabel('$\theta$','Interpreter','latex','FontSize',16)
    box on
    grid on
end

end

function magneticMomentDiagnostic(ST,numBins)
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
set(gcf,'name','Magnetic moment','numbertitle','off')
surf(time,vals,log10(f),'LineStyle','none')
% surf(time,vals,f,'LineStyle','none')
axis([tmin tmax minVal maxVal])
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('Magnetic moment $\mu$ (arbitrary units)','Interpreter','latex','FontSize',16)
colormap(jet)
end

function poloidalPlots(ST)



end