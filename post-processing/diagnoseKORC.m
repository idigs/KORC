function ST = diagnoseKORC(path)
close all

ST = struct;
ST.path = path;

ST.params = loadSimulationParameters(ST);

ST.data = loadData(ST);

% energyConservation(ST);

% ST.RT = radialTransport(ST);

% confined_particles(ST);

% pitchAngleDiagnostic(ST,100);

% magneticMomentDiagnostic(ST,50);

% poloidalPlaneDistributions(ST,25);

% angularMomentum(ST);

changeOfMagneticField(ST)

% energyLimit(ST);

% LarmorVsLL(ST);


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

% params.simulation.num_snapshots = 250;
% params.simulation.t_steps = params.simulation.output_cadence*params.simulation.num_snapshots;

end

function data = loadData(ST)
data = struct;

list = {'X','V','B'};

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    for ss=1:ST.params.simulation.num_species
        tnp = double(ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
        
        data.(['sp' num2str(ss)]).(list{ll}) = ...
            zeros(3,tnp,ST.params.simulation.num_snapshots+1);
        
        for ff=1:ST.params.simulation.nmpi
            filename = [ST.path 'file_' num2str(ff-1) '.h5'];
            indi = (ff - 1)*double(ST.params.species.ppp(ss)) + 1;
            indf = ff*double(ST.params.species.ppp(ss));
            for ii=1:(ST.params.simulation.num_snapshots+1)
                dataset = ...
                    ['/' num2str((ii-1)*double(ST.params.simulation.output_cadence)) '/spp_' num2str(ss)...
                    '/' list{ll}];
                
                data.(['sp' num2str(ss)]).(list{ll})(:,indi:indf,ii) = ...
                    h5read(filename, dataset);
            end
            
        end
        
    end
end


% list = {'eta','gamma','Prad','Pin','flag','mu'};
list = {'eta','gamma','flag','mu'};

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    for ss=1:ST.params.simulation.num_species
        tnp = double(ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
        
        data.(['sp' num2str(ss)]).(list{ll}) = ...
            zeros(tnp,ST.params.simulation.num_snapshots);
        
        for ii=1:ST.params.simulation.num_snapshots            
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
maxerr = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);
minerr = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);

st1 = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);
st2 = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);
st3 = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);
st4 = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);

h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
h5 = figure;
h6 = figure;

set(h1,'name','Energy conservation','numbertitle','off')
set(h2,'name','Velocity components','numbertitle','off')
set(h3,'name','Radiated power','numbertitle','off')
set(h4,'name','Input power','numbertitle','off')
set(h5,'name','Energy gain/loss','numbertitle','off')
set(h6,'name','Energy statistics','numbertitle','off')
try
    for ss=1:ST.params.simulation.num_species
        pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
%         passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
%         bool = pin & passing;
        gamma = ST.data.(['sp' num2str(ss)]).gamma(pin,:);
        tmp = zeros(size(gamma));
        for ii=1:size(tmp,1)
            tmp(ii,:) = ...
                100*( gamma(ii,:) - gamma(ii,1) )./gamma(ii,1);
%             tmp(ii,:) = ST.data.(['sp' num2str(ss)]).gamma(ii,:)./ST.data.(['sp' num2str(ss)]).gamma(ii,1);
        end
        err(:,ss) = mean(tmp,1);
%         maxerr(:,ss) = max(tmp,[],1);
%         minerr(:,ss) = min(tmp,[],1);
        maxerr(:,ss) = err(:,ss) + std(tmp,0,1)';
        minerr(:,ss) = err(:,ss) - std(tmp,0,1)';
        
        if (~isempty(tmp))
            st1(:,ss) = mean(tmp,1);
            st2(:,ss) = std(tmp,0,1);
            st3(:,ss) = skewness(tmp,1,1);
            st4(:,ss) = kurtosis(tmp,1,1);
        end
        
        Prad = mean(abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:)),1);
        minPrad = Prad + std(abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:)),0,1);
        maxPrad = Prad - std(abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:)),0,1);
        
        Pin = mean(abs(ST.data.(['sp' num2str(ss)]).Pin(pin,:)),1);
        minPin = Pin + std(abs(ST.data.(['sp' num2str(ss)]).Pin(pin,:)),0,1);
        maxPin = Pin - std(abs(ST.data.(['sp' num2str(ss)]).Pin(pin,:)),0,1);
        
        eta = ST.data.(['sp' num2str(ss)]).eta(pin,:);
        V = sqrt(1 - 1./gamma.^2); % in units of the speed of light
        Vpar = mean((V.*cos(eta)).^2,1);
        Vperp = mean((V.*sin(eta)).^2,1);
        
        figure(h1)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(time,err(:,ss),'k-',time,minerr(:,ss),'r:',time,maxerr(:,ss),'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\Delta \mathcal{E}/\mathcal{E}_0$ (\%)','Interpreter','latex','FontSize',16)
        
        figure(h2)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(time,Vpar-Vperp,'k-')
%         plot(time,Vpar,'k-',time,Vperp,'r-')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$v_\parallel$, $v_\perp$ ($c$)','Interpreter','latex','FontSize',16)

        figure(h6)
        subplot(4,1,1)
        hold on; plot(time,st1(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\mu$','Interpreter','latex','FontSize',16)
        subplot(4,1,2)
        hold on; plot(time,st2(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\sigma$','Interpreter','latex','FontSize',16)
        subplot(4,1,3)
        hold on; plot(time,st3(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$s$','Interpreter','latex','FontSize',16)
        subplot(4,1,4)
        hold on; plot(time,st4(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$k$','Interpreter','latex','FontSize',16)
        
        
        figure(h3)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(time,Prad,'k',time,minPrad,'r:',time,maxPrad,'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{rad} \rangle$','Interpreter','latex','FontSize',16)
        
        figure(h4)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(time,Pin,'k',time,minPin,'r:',time,maxPin,'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{in} \rangle$','Interpreter','latex','FontSize',16)
        
        figure(h5)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(time,Prad./Pin,'k')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{rad} \rangle / \langle P_{in} \rangle$','Interpreter','latex','FontSize',16)
    end
catch
    error('Something went wrong: energyConservation')
end


end

function pitchAngleDiagnostic(ST,numBins)
N = 10;
nbins = 30;

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

mean_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
std_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
skewness_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
kurtosis_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);

fx = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots);
x = zeros(ST.params.simulation.num_species,nbins);

f_tot = zeros(numBins,ST.params.simulation.num_snapshots);

data = [];

h1 = figure;
set(h1,'name','PDF time evolution','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
%     tmp = cos(pi*ST.data.(['sp' num2str(ss)]).eta/180);
    tmp = ST.data.(['sp' num2str(ss)]).eta(pin,:);
    
    mean_f(ss,:) = mean(tmp,1);
    std_f(ss,:) = std(tmp,0,1);
    skewness_f(ss,:) = skewness(tmp,0,1);
    kurtosis_f(ss,:) = kurtosis(tmp,1,1);
    
    data = [data; ST.data.(['sp' num2str(ss)]).eta];
    
    
    minVal = min(min( data ));
    maxVal = max(max( data ));
    x(ss,:) = linspace(minVal,maxVal,nbins);
    
    for ii=1:ST.params.simulation.num_snapshots
        [fx(ss,:,ii),~] = hist(tmp(:,ii),x(ss,:));
    end
    
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    surf(time,squeeze(x(ss,:)),log10(squeeze(fx(ss,:,:))),'LineStyle','none')
    axis([tmin tmax minVal maxVal])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('Pitch angle $\theta$ (degrees)','Interpreter','latex','FontSize',16)
    colormap(jet)
end

minVal = min(min( data ));
maxVal = max(max( data ));
vals = linspace(minVal,maxVal,numBins);

for ii=1:ST.params.simulation.num_snapshots
    [f_tot(:,ii),~] = hist(data(:,ii),vals);
end


% h1 = figure;
% set(h1,'name','PDF time evolution','numbertitle','off')
% surf(time,vals,log10(f_tot),'LineStyle','none')
% % surf(time,vals,f,'LineStyle','none')
% axis([tmin tmax minVal maxVal])
% box on
% axis on
% xlabel('Time (s)','Interpreter','latex','FontSize',16)
% ylabel('Pitch angle $\theta$ (degrees)','Interpreter','latex','FontSize',16)
% colormap(jet)

h2 = figure
set(h2,'name','Statistical moments pitch angle','numbertitle','off')
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
    nc = floor(N/2);
    nr = floor(N/nc);
    subplot(nr,nc,ii)
    plot(z,log10(fz),'k')
    hold on
    for ss=1:ST.params.simulation.num_species
        dx = mean(diff(x(ss,:)));
        xAxis = ( x(ss,:) - mean_f(ss,it) )/std_f(ss,it);
        f = std_f(ss,it)*fx(ss,:,it)/(sum(fx(ss,:,it))*dx);
        
        figure(h3)
        plot(xAxis,log10(f),'o:')
        
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


% % Poloidal angle distribution function % %

ft = zeros(ST.params.simulation.num_species,nbins,N);
t = zeros(ST.params.simulation.num_species,nbins,N);

for ii=1:N
    it = ii*offset;
    for ss=1:ST.params.simulation.num_species
%         I = find( any(ST.data.sp3.eta > 90, 2) == 0 );
%         X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,I,it));
        X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,:,it));
        theta = atan2(X(3,:), sqrt(X(1,:).^2 + X(2,:).^2) - ST.params.fields.Ro);
        theta(theta < 0) = theta(theta < 0) + 2*pi;
        theta = (180/pi)*theta;
        [ft(ss,:,ii),t(ss,:,ii)] = hist(theta,nbins);
        dt = mean(diff(t(ss,:,ii)));
        ft(ss,:,ii) = ft(ss,:,ii)/(sum(ft(ss,:,ii))*dt);
    end
end

barcolor = [1,0,0;0,1,0;0,0,1;0.5,0.5,1.0;1.0,0.2,0.2];

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
cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

mean_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
std_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
skewness_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);
kurtosis_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots);

fx = zeros(ST.params.simulation.num_species,numBins,ST.params.simulation.num_snapshots);
x = zeros(ST.params.simulation.num_species,numBins);

h1 = figure;
set(h1,'name','PDF of magnetic moment','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    tmp = ST.data.(['sp' num2str(ss)]).mu(pin,:);
    for ii=1:size(tmp,1)
        tmp(ii,:) = 100*(tmp(ii,:) - tmp(ii,1))./tmp(ii,1);
    end
    
    mean_f(ss,:) = mean(tmp,1);
    std_f(ss,:) = std(tmp,0,1);
    skewness_f(ss,:) = skewness(tmp,0,1);
    kurtosis_f(ss,:) = kurtosis(tmp,1,1);
      
    minVal = min(min( tmp ));
    maxVal = max(max( tmp ));
    x(ss,:) = linspace(minVal,maxVal,numBins);
    
    for ii=1:ST.params.simulation.num_snapshots
        try
            [fx(ss,:,ii),~] = hist(tmp(:,ii),x(ss,:));
        catch
        end
    end
    
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    surf(time,squeeze(x(ss,:)),log10(squeeze(fx(ss,:,:))),'LineStyle','none')
%     axis([tmin tmax minVal maxVal])
    view([0,90])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\Delta \mu/\mu_0$ ($\%$)','Interpreter','latex','FontSize',16)
    colormap(jet)
end

h2 = figure;
set(h2,'name','Statistical moments magnetic moment','numbertitle','off')
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

end

function poloidalPlaneDistributions(ST,nbins)
N = 3;
offset = floor(double(ST.params.simulation.num_snapshots)/N);

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);

% m = figure;


for ss=1:ST.params.simulation.num_species
    for ii=1:N
        it = ii*offset;
        
        R = squeeze( sqrt( ST.data.(['sp' num2str(ss)]).X(1,:,it).^2 + ...
            ST.data.(['sp' num2str(ss)]).X(2,:,it).^2 ) );
        Z = squeeze( ST.data.(['sp' num2str(ss)]).X(3,:,it) );
        Prad = squeeze( ST.data.(['sp' num2str(ss)]).Prad(:,it) );

        % Poloidal distribution of particles
        h = figure;
        subplot(3,1,1)
        n = histogram2(R,Z,[nbins,nbins],'FaceColor','flat','Normalization','probability');
        colorbar
        xAxis = n.XBinEdges;
        dx = mean( diff(xAxis) );
        yAxis = n.YBinEdges;
        dy = mean( diff(yAxis) );
        axis([min(xAxis) max(xAxis) min(yAxis) max(yAxis)])
        axis equal
        view([0 90])
        colormap(jet)
        title(['Species: ' num2str(ss) 'Time: ' num2str(time(it))],'Interpreter','latex','FontSize',11)
        xlabel('$R$','Interpreter','latex','FontSize',16)
        ylabel('$Z$','Interpreter','latex','FontSize',16)
        
        x = ST.params.fields.Ro + ...
            ST.params.species.r(ss)*cos(linspace(0,2*pi,100));
        y = ST.params.species.r(ss)*sin(linspace(0,2*pi,100));
        hold on; plot3(x,y,1*ones(size(x)),'k');hold off
        
        
        % Poloidal distribution of radiated synchroton power
%        prad = histogram2(R,Z,[nbins,nbins],'FaceColor','flat')       
        prad = zeros(nbins,nbins);
        
        R = R - min(xAxis);
        Z = Z + abs(min(yAxis));
        indx = floor( (R + 0.5*dx)/dx ) + 1;
        indy = floor( (Z + 0.5*dy)/dy ) + 1;
        
        I = find(indx > nbins);
        indx(I) = indx(I) - 1;
        
        I = find(indy > nbins);
        indy(I) = indy(I) - 1;
        
        for jj=1:ST.params.species.ppp(ss)
            prad(indy(jj),indx(jj)) = prad(indx(jj),indy(jj)) + Prad(jj);
        end
        
        figure(h)
        subplot(3,1,2)
        surf(xAxis(1:end-1)',yAxis(1:end-1)',prad,'LineStyle','none')
        axis equal
        view([0 90])        
        colormap(jet)
        axis([min(xAxis) max(xAxis) min(yAxis) max(yAxis)])
        xlabel('$R$','Interpreter','latex','FontSize',16)
        ylabel('$Z$','Interpreter','latex','FontSize',16)
        
        A = n.Values'.*prad;
%         I = isinf(A);
%         A(I) = 0;
        
        figure(h)
        subplot(3,1,3)
        surf(xAxis(1:end-1)',yAxis(1:end-1)',A,'LineStyle','none')
        axis equal
        view([0 90])        
        colormap(jet)
        axis([min(xAxis) max(xAxis) min(yAxis) max(yAxis)])
        xlabel('$R$','Interpreter','latex','FontSize',16)
        ylabel('$Z$','Interpreter','latex','FontSize',16)

    end

% for ii=1:N
%     it = ii*offset;
%     
%     R = squeeze( sqrt( ST.data.(['sp' num2str(ss)]).X(1,:,it).^2 + ...
%         ST.data.(['sp' num2str(ss)]).X(2,:,it).^2 ) );
%     Z = squeeze( ST.data.(['sp' num2str(ss)]).X(3,:,it) );
%     Prad = squeeze( ST.data.(['sp' num2str(ss)]).Prad(:,it) );
%     
%     % Poloidal distribution of particles
%     figure(m)
%     n = histogram2(R,Z,[nbins,nbins],'FaceColor','flat','Normalization','probability');
%     colorbar
%     xAxis = n.XBinEdges;
%     dx = mean( diff(xAxis) );
%     yAxis = n.YBinEdges;
%     dy = mean( diff(yAxis) );
%     axis([min(xAxis) max(xAxis) min(yAxis) max(yAxis)])
%     axis equal
%     view([0 90])
%     colormap(jet)
%     title(['Species: ' num2str(ss) 'Time: ' num2str(time(it))],'Interpreter','latex','FontSize',11)
%     xlabel('$R$','Interpreter','latex','FontSize',16)
%     ylabel('$Z$','Interpreter','latex','FontSize',16)
%     
%     x = ST.params.fields.Ro + ...
%         ST.params.species.r(ss)*cos(linspace(0,2*pi,100));
%     y = ST.params.species.r(ss)*sin(linspace(0,2*pi,100));
%     hold on; plot3(x,y,1*ones(size(x)),'k');hold off
%     
%     F(ii) = getframe(gcf);
% end
end

end

function angularMomentum(ST)
c = 299792458.0;
c = 1E2*c;

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);

Bo = 1E4*ST.params.fields.Bo;
Ro = 1E2*ST.params.fields.Ro; % Major radius in meters.
a = 1E2*ST.params.fields.a;% Minor radius in meters.
co = 0.5; % Extra parameter
lambda = a/co;
Bpo = 1E4*ST.params.fields.Bpo;

st1 = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);
st2 = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);
st3 = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);
st4 = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);

h = figure;
set(h,'name','Angular momentum conservation','numbertitle','off')
h1 = figure;
set(h1,'name','Angular momentum statistics','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    num_part = numel(find(pin==1));
%     num_part = ST.params.species.ppp(ss)*ST.params.simulation.nmpi;
    
    m = 1E3*ST.params.species.m(ss);
    q = 3E9*ST.params.species.q(ss);
    
    invariant = ...
        zeros(num_part,ST.params.simulation.num_snapshots);
    err = zeros(1,ST.params.simulation.num_snapshots);
    minerr = zeros(1,ST.params.simulation.num_snapshots);
    maxerr = zeros(1,ST.params.simulation.num_snapshots);
    
    for ii=1:ST.params.simulation.num_snapshots
        X = 1E2*squeeze( ST.data.(['sp' num2str(ss)]).X(:,pin,ii) );
        V = 1E2*squeeze( ST.data.(['sp' num2str(ss)]).V(:,pin,ii) );

        % Toroidal coordinates
        % r = radius, theta = poloidal angle, phi = toroidal angle
        r = sqrt( (sqrt(X(1,:).^2 + X(2,:).^2) - Ro).^2 + X(3,:).^2 );
        theta = atan2(X(3,:),sqrt(X(1,:).^2 + X(2,:).^2) - Ro);
        theta(theta<0) = theta(theta<0) + 2*pi;
        zeta = atan2(X(1,:),X(2,:));
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        % Toroidal coordinates

        gamma = squeeze( ST.data.(['sp' num2str(ss)]).gamma(pin,ii) )';

        eta = r/Ro;
        psi = 0.5*lambda*Bpo*log(1 + r.^2/lambda^2);
        wo = q*Bo./(m*c*gamma);

        dzeta = ...
        (X(2,:).*V(1,:) - X(1,:).*V(2,:))./( sum(X(1:2,:).^2,1) );
    
        invariant(:,ii) = dzeta.*( 1 + eta.*cos(theta) ).^2 - wo.*psi/(Ro*Bo);
        tmp_vec = 100*(invariant(:,ii) - invariant(:,1))./invariant(:,1);
        err(ii) = mean( tmp_vec );
%         minerr(ii) = min( tmp_vec );
%         maxerr(ii) = max( tmp_vec );
        minerr(ii) = err(ii) - std( tmp_vec );
        maxerr(ii) = err(ii) + std( tmp_vec );
        
        if (~isempty(tmp_vec))
            st1(ii,ss) = mean(tmp_vec);
            st2(ii,ss) = std(tmp_vec);
            st3(ii,ss) = skewness(tmp_vec);
            st4(ii,ss) = kurtosis(tmp_vec);
        end
    end
    
    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(time,err,'k-',time,minerr,'r:',time,maxerr,'r:')
    box on
    grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',12)
    ylabel('Angular momentum conservation (\%)','Interpreter','latex','FontSize',12)
    
    figure(h1)
    subplot(4,1,1)
    hold on; plot(time,st1(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\mu$','Interpreter','latex','FontSize',16)
    subplot(4,1,2)
    hold on; plot(time,st2(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\sigma$','Interpreter','latex','FontSize',16)
    subplot(4,1,3)
    hold on; plot(time,st3(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$s$','Interpreter','latex','FontSize',16)
    subplot(4,1,4)
    hold on; plot(time,st4(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$k$','Interpreter','latex','FontSize',16)
    
end


end

function changeOfMagneticField(ST)
cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

for ss=1:ST.params.simulation.num_species   
    q = ST.params.species.q(ss);
    m = ST.params.species.m(ss);
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    aux = find(pin == 1);
    
    Bo = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,pin,1).^2,1) ) );
    gamma = ST.data.(['sp' num2str(ss)]).gamma(pin,1)';
    wc = abs(q)*Bo./(gamma*m);
    Tc = 2*pi./wc;
    I = zeros(size(Bo));
    B = zeros(size(Bo));
    DB = zeros(size(Bo));
    for ii=1:numel(Bo)
        [~,I(ii)] = min( abs(Tc(ii) - time) );
        tmp = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),1:I(ii)).^2,1) ) );
        DB(ii) = max( 100*abs(tmp - Bo(ii))/Bo(ii) );
        
%         B = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),I(ii)).^2,1) ) );
%         DB(ii) = 100*abs(B - Bo(ii))/Bo(ii);
        
%         B1 = ST.data.(['sp' num2str(ss)]).B(:,aux(ii),1);
%         B2 = ST.data.(['sp' num2str(ss)]).B(:,aux(ii),I(ii));
%         DB(ii) = 100*sqrt(sum((B2-B1).^2))/sqrt(sum(B1.^2));
        
%         B(ii) = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),I(ii)).^2,1) ) );
    end
%     DB = 100*abs( B - Bo )./Bo;

    minVal = min( DB );
    maxVal = max( DB );
    x = linspace(minVal,maxVal,50);

    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,pin,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);

    S = 12*ones(size(gamma));
    h=figure;
    set(h,'name',['Change of B-field: sp' num2str(ss)],'numbertitle','off')
%     scatter3(R,Z,DB,S,DB,'square','filled')
    subplot(1,2,1)
    histogram(DB,x)
    subplot(1,2,2)
    scatter3(R,Z,DB,S,DB,'square','filled')
    colormap(jet(256))
    colorbar
    view([0,90])
    box on; axis on; axis square
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
end



end

function RT = radialTransport(ST)
RT = struct;

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
Ro = ST.params.fields.Ro;
rc = zeros(1,ST.params.simulation.num_species);

C = colormap(jet(512));
offset = floor(512/ST.params.simulation.num_species);
colour = C(1:offset:end,:);

h1=figure;
set(h1,'name','IC','numbertitle','off')
legends = cell(1,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species
%     h1=figure;
%     set(h1,'name',['Species ' num2str(ss)],'numbertitle','off')
    pin = zeros(1,ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
    for pp=1:ST.params.species.ppp(ss)*ST.params.simulation.nmpi
        X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,pp,:));
        
        % Toroidal coordinates
        % r = radius, theta = poloidal angle, zeta = toroidal angle
        r = sqrt( (sqrt(sum(X(1:2,:).^2,1)) - Ro).^2 + X(3,:).^2 );
        theta = atan2(X(3,:),sqrt(sum(X(1:2,:).^2,1)) - Ro);
        theta(theta<0) = theta(theta<0) + 2*pi;
        theta = 180*theta/pi;
        zeta = atan2(X(1),X(2));
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        % Toroidal coordinates
        
        bool = all(ST.data.(['sp' num2str(ss)]).eta(pp,:) < 90);
        
        if all(r < ST.params.fields.a) && bool && all
%             rc(ss) = r(1);
            pin(pp) = 1;
        end
%        figure(h1)
%        hold on
%        plot(theta,r,'.')
%        hold off
    end
    
    t = linspace(0,2*pi,200);
    Rs = ST.params.fields.Ro + ST.params.fields.a*cos(t);
    Zs = ST.params.fields.a*sin(t);

    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,logical(pin),1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);
    
    figure(h1)
%     subplot(2,1,1)
    hold on
    plot(R,Z,'k.','MarkerSize',12,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
    hold off
%     hold on
%     plot(Rs,Zs,'r')
%     box on
%     axis on
%     axis equal
%     xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
%     ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

    
%     X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,logical(pin),end));
%     R = sqrt( sum(X(1:2,:).^2,1) );
%     Z = X(3,:);
%     
%     subplot(2,1,2)
%     plot(R,Z,'b.',Rs,Zs,'r')
%     box on
%     axis on
%     axis equal
%     xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
%     ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
% subplot(2,1,1)
legend(legends)
hold on
plot(Rs,Zs,'r')
box on
axis on
axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

RT.rc = rc;

end

function confined_particles(ST)
RT = struct;

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

C = colormap(jet(512));
offset = floor(512/ST.params.simulation.num_species);
colour = C(1:offset:end,:);

t = linspace(0,2*pi,200);
Rs = ST.params.fields.Ro + ST.params.fields.a*cos(t);
Zs = ST.params.fields.a*sin(t);

h0 = figure;
set(h0,'name','Particle loss','numbertitle','off')
legends = cell(1,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species
    confinedParticles = ...
        100*sum(ST.data.(['sp' num2str(ss)]).flag,1)/size(ST.data.(['sp' num2str(ss)]).flag,1);
    figure(h0)
    hold on
    plot(time,confinedParticles)
    hold off
    
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h0)
legend(legends,'interpreter','latex','FontSize',12)
box on
axis on
axis square
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
ylabel('$\%$ of confined RE','Interpreter','latex','FontSize',16)


% h = figure;
% set(h,'name','Spatial distribution particle loss','numbertitle','off')
% legends = cell(1,ST.params.simulation.num_species);
% for ss=1:ST.params.simulation.num_species
%     pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
%     I = find(pin==0);
%     for ii=1:numel(I)
%         it = find(ST.data.(['sp' num2str(ss)]).flag(I(ii),:)==0,1,'first');
%         
%         X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,I(ii),it));
%         R = sqrt( sum(X(1:2).^2,1) );
%         Z = X(3);
%         
%         figure(h)
%         subplot(1,double(ST.params.simulation.num_species),double(ss))
%         hold on
%         plot(R,Z,'o','MarkerSize',6,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
%         hold off
%     end
%     figure(h)
%     subplot(1,double(ST.params.simulation.num_species),double(ss))
%     hold on
%     plot(Rs,Zs,'r')
%     hold off
%     box on
%     axis on
%     axis equal
%     xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
%     ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
% end

% % Toroidal coordinates
%         % r = radius, theta = poloidal angle, zeta = toroidal angle
%         r = sqrt( (sqrt(sum(X(1:2,:).^2,1)) - Ro).^2 + X(3,:).^2 );
%         theta = atan2(X(3,:),sqrt(sum(X(1:2,:).^2,1)) - Ro);
%         theta(theta<0) = theta(theta<0) + 2*pi;
%         theta = 180*theta/pi;
%         zeta = atan2(X(1),X(2));
%         zeta(zeta<0) = zeta(zeta<0) + 2*pi;
%         % Toroidal coordinates
%         
%         
% h = figure;
% set(h,'name','Spatial distribution particle loss','numbertitle','off')
% legends = cell(1,ST.params.simulation.num_species);
% for ss=1:ST.params.simulation.num_species
%     pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
%     I = find(pin==0);
%     for ii=1:numel(I)
%         it = find(ST.data.(['sp' num2str(ss)]).flag(I(ii),:)==0,1,'first');
%         
%         X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,I(ii),it));
%         R = sqrt( sum(X(1:2).^2,1) );
%         Z = X(3);
%         
%         figure(h)
%         subplot(1,double(ST.params.simulation.num_species),double(ss))
%         hold on
%         plot(R,Z,'o','MarkerSize',6,'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0])
%         hold off
%     end
%     figure(h)
%     subplot(1,double(ST.params.simulation.num_species),double(ss))
%     hold on
%     plot(Rs,Zs,'r')
%     hold off
%     box on
%     axis on
%     axis equal
%     xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
%     ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
% end


h1=figure;
set(h1,'name','IC','numbertitle','off')
legends = cell(1,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species   
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin & passing;
    
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);

    figure(h1)
    subplot(1,2,1)
    hold on
    plot(R,Z,'s','MarkerSize',4,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
%     plot(R,Z,'.','MarkerSize',10,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
    hold off
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
subplot(1,2,1)
legend(legends,'interpreter','latex','FontSize',12)
hold on
plot(Rs,Zs,'r')
hold off
box on
axis on
axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

legends = cell(1,ST.params.simulation.num_species);
for ss=ST.params.simulation.num_species:-1:1
    t = linspace(0,2*pi,200);
    Rs = ST.params.fields.Ro + ST.params.fields.a*cos(t);
    Zs = ST.params.fields.a*sin(t);

    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    trapped = logical( any(ST.data.(['sp' num2str(ss)]).eta > 90,2) );

    bool = pin & trapped;
    
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);
    
    figure(h1)
    subplot(1,2,2)
    hold on
    plot(R,Z,'s','MarkerSize',6,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
%     plot(R,Z,'.','MarkerSize',10,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
    hold off
    legends{ST.params.simulation.num_species + 1 -ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
subplot(1,2,2)
legend(legends,'interpreter','latex','FontSize',12)
hold on
plot(Rs,Zs,'r')
box on
axis on
axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

end

function energyLimit(ST)
cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

h1=figure;
set(h1,'name','Energy vs. pitch angle','numbertitle','off')
h2=figure;
set(h2,'name','Energy vs. time','numbertitle','off')
legends = cell(1,ST.params.simulation.num_species);

for ss=1:ST.params.simulation.num_species   
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    
    energy = mean(ST.data.(['sp' num2str(ss)]).gamma(pin,:),1);
    energy = ST.params.species.m(ss)*ST.params.scales.v^2*energy/ST.params.scales.q;
    energy = energy/1E6;
    pitch = mean(ST.data.(['sp' num2str(ss)]).eta(pin,:),1);
    
    figure(h1)
    hold on
%     plot(pitch(1),energy(1),'ko',pitch,energy)
    plot(pitch,energy)
    hold off
        
    figure(h2)
    hold on
    plot(time,energy)
    hold off
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
legend(legends,'interpreter','latex','FontSize',12)
grid on; box on; axis on; axis square
xlabel('Pitch angle $\eta$ ($^\circ$)','Interpreter','latex','FontSize',16)
ylabel('$\mathcal{E}_0$ (MeV)','Interpreter','latex','FontSize',16)

figure(h2)
legend(legends,'interpreter','latex','FontSize',12)
grid on; box on; axis on; axis square
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
ylabel('$\mathcal{E}_0$ (MeV)','Interpreter','latex','FontSize',16)
end

function B = analyticalB(ST,X)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

narginchk(1,2);

% Parameters of the analytical magnetic field
Bo = ST.params.fields.Bo;
a = ST.params.fields.a;
Ro = ST.params.fields.Ro;
qa = ST.params.fields.qa;
co = ST.params.fields.co;
lamb = ST.params.fields.lambda;
Bpo = ST.params.fields.Bpo;
% Parameters of the analytical magnetic field


% Toroidal coordinates
% r = radius, theta = poloidal angle, zeta = toroidal angle
r = squeeze(sqrt( (sqrt(sum(X(1:2,:,:).^2,1)) - Ro).^2 + X(3,:,:).^2));
theta = atan2(squeeze(X(3,:,:)),squeeze(sqrt(sum(X(1:2,:,:).^2,1)) - Ro));
theta(theta<0) = theta(theta<0) + 2*pi;
zeta = atan2(squeeze(X(1,:,:)),squeeze(X(2,:,:)));
zeta(zeta<0) = zeta(zeta<0) + 2*pi;
% Toroidal coordinates

% Poloidal magnetic field
% Minus sign = TEXTOR
% Plus sign = default
Bp = Bpo*(r/lamb)./( 1 + (r/lamb).^2 );

eta = r/Ro;
Br = 1./( 1 + eta.*cos(theta) );

Bx = Br.*( Bo*cos(zeta) - Bp.*sin(theta).*sin(zeta) );
By = -Br.*( Bo*sin(zeta) + Bp.*sin(theta).*cos(zeta) );
Bz = Br.*Bp.*cos(theta);

B = zeros(3,size(r,1),size(r,2));
B(1,:,:) = Bx;
B(2,:,:) = By;
B(3,:,:) = Bz;

end

function E = analyticalE(ST,X)
% Analytical magnetic field
% X is a vector X(1)=x, X(2)=y, X(3)=z.

narginchk(1,2);

% Parameters of the analytical magnetic field
Eo = ST.params.fields.Eo;
Ro = ST.params.fields.Ro;
% Parameters of the analytical magnetic field


% Toroidal coordinates
% r = radius, theta = poloidal angle, zeta = toroidal angle
r = squeeze(sqrt( (sqrt(sum(X(1:2,:,:).^2,1)) - Ro).^2 + X(3,:,:).^2));
theta = atan2(squeeze(X(3,:,:)),squeeze(sqrt(sum(X(1:2,:,:).^2,1)) - Ro));
theta(theta<0) = theta(theta<0) + 2*pi;
zeta = atan2(squeeze(X(1,:,:)),squeeze(X(2,:,:)));
zeta(zeta<0) = zeta(zeta<0) + 2*pi;
% Toroidal coordinates

% Poloidal magnetic field
eta = r/Ro;
Ezeta = Eo./( 1 + eta.*cos(theta) );

Ex = Ezeta.*cos(zeta);
Ey = -Ezeta.*sin(zeta);
Ez = 0;

E = zeros(3,size(r,1),size(r,2));
E(1,:,:) = Ex;
E(2,:,:) = Ey;
E(3,:,:) = Ez;

end

function LarmorVsLL(ST)
cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
% % % % % % % % % % %


% 
%         % % % % % % % % % % % % % % %
%         % Radiation losses operator
%         vmag = sqrt( V*V' );
%         aux =  cross(V,E) + V*(V*B') - B*vmag^2;
%         curv = abs(q)*sqrt( aux*aux' )/(gamma*m*vmag^3);
%         % Radiation losses operator
%         % % % % 

h=figure;
set(h,'name','Radiation: Larmor approx','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    q = ST.params.species.q(ss);
    m = ST.params.species.m(ss);
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    
    V = ST.data.(['sp' num2str(ss)]).V(:,pin,:);
    gamma = ST.data.(['sp' num2str(ss)]).gamma(pin,:);
    X = ST.data.(['sp' num2str(ss)]).X(:,pin,:);
    B = analyticalB(ST,X);
    E = analyticalE(ST,X);
    
    vmag = squeeze(sqrt(sum(V.^2,1)));
    VxE = zeros(size(X));
    VxVxB = zeros(size(X));
    for it=1:size(X,3)
        VxE(:,:,it) = cross(V(:,:,it),E(:,:,it));
        VxVxB(:,:,it) = cross(V(:,:,it),cross(V(:,:,it),B(:,:,it)));
    end
    clear E B
    
    kappa = q*sqrt( squeeze(sum((VxE + VxVxB).^2,1)) )./(m*gamma.*vmag.^3);
    clear VxE VxVxB
    
    PR = 2*Kc*q^2*mean((gamma.*vmag).^4.*kappa.^2,1)/(3*c^3);
    
    eta = ST.data.(['sp' num2str(ss)]).eta(pin,:);
    v = c*sqrt(1 - 1./gamma.^2); % in units of the speed of light
    vperp = v.*sin(eta);
    wc = q*ST.params.fields.Bo./(gamma*m);
    rg = vperp./wc;
    kappa2 = 1/ST.params.fields.Ro.^2 + ...
        sin(eta).^4./rg.^2;
    
    PR_approx = 2*Kc*q^2*mean((gamma.*vmag).^4.*kappa2,1)/(3*c^3);
    
    
    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
%     semilogy(time,PR,'k',time,PR_approx,'r')
    plot(time,PR./PR_approx,'k')
    box on; grid on;
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\langle P_{R} \rangle$/$\langle P_{RL} \rangle$','Interpreter','latex','FontSize',16)
end

end