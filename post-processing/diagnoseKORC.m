function ST = diagnoseKORC(path)
% close all

ST = struct;
ST.path = path;

ST.params = loadSimulationParameters(ST);

ST.data = loadData(ST);

energyConservation(ST);

pitchAngleDiagnostic(ST,100);

% magneticMomentDiagnostic(ST,100);

% poloidalPlaneDistributions(ST,25);

% angularMomentum(ST);

% changeOfMagneticField(ST)


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

% params.simulation.num_snapshots = 577;
% params.simulation.t_steps = 18464000;

end

function data = loadData(ST)
data = struct;

list = {'X','V'};%,'Rgc'};

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
list = {'eta','gamma','mu','Prad','Pin'};

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
maxerr = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);
minerr = zeros(ST.params.simulation.num_snapshots,ST.params.simulation.num_species);

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(cad:cad:ST.params.simulation.t_steps);

h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
h5 = figure;
set(h1,'name','Energy conservation','numbertitle','off')
% set(h2,'name','Relativistic energy','numbertitle','off')
set(h3,'name','Radiated power','numbertitle','off')
set(h4,'name','Input power','numbertitle','off')
set(h5,'name','Energy gain/loss','numbertitle','off')
try
    for ss=1:ST.params.simulation.num_species
        tmp = zeros(size(ST.data.(['sp' num2str(ss)]).gamma));
        for ii=1:ST.params.species.ppp(ss)*ST.params.simulation.nmpi
            tmp(ii,:) = ...
                100*( ST.data.(['sp' num2str(ss)]).gamma(ii,:) - ...
                ST.data.(['sp' num2str(ss)]).gamma(ii,1) )./ST.data.(['sp' num2str(ss)]).gamma(ii,1);
        end
        err(:,ss) = mean(tmp,1);
        maxerr(:,ss) = max(tmp,[],1);
        minerr(:,ss) = min(tmp,[],1);
        figure(h1)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(time,err(:,ss),'k-',time,minerr(:,ss),'r:',time,maxerr(:,ss),'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\Delta E/E_0$ (\%)','Interpreter','latex','FontSize',16)
        
%         figure(h2)
%         subplot(double(ST.params.simulation.num_species),1,double(ss))
%         Erel = (299792458.0)^2*ST.params.species.m(ss)*mean(ST.data.(['sp' num2str(ss)]).gamma,1)/abs(ST.params.species.q(ss));
%         plot(time,Erel)
%         box on
%         grid on
%         xlabel('Time (s)','Interpreter','latex','FontSize',16)
%         ylabel('$\langle \gamma m_0 c^2 \rangle$','Interpreter','latex','FontSize',16)
        
        figure(h3)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        Prad = mean(abs(ST.data.(['sp' num2str(ss)]).Prad),1);
        minPrad = min( abs(ST.data.(['sp' num2str(ss)]).Prad), [], 1 );
        maxPrad = max( abs(ST.data.(['sp' num2str(ss)]).Prad), [], 1 );
        plot(time,Prad,'k',time,minPrad,'r:',time,maxPrad,'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{rad} \rangle$','Interpreter','latex','FontSize',16)
        
        figure(h4)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        Pin = mean(abs(ST.data.(['sp' num2str(ss)]).Pin),1);
        minPin = min( abs(ST.data.(['sp' num2str(ss)]).Pin), [], 1 );
        maxPin = max( abs(ST.data.(['sp' num2str(ss)]).Pin), [], 1 );
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
time = ST.params.simulation.dt*double(cad:cad:ST.params.simulation.t_steps);
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
%     tmp = cos(pi*ST.data.(['sp' num2str(ss)]).eta/180);
    tmp = ST.data.(['sp' num2str(ss)]).eta;
    
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
tmp = [];

for ss=1:ST.params.simulation.num_species
    tmp = [tmp;ST.data.(['sp' num2str(ss)]).mu];
end


f = zeros(numBins,ST.params.simulation.num_snapshots);

for ii=2:ST.params.simulation.num_snapshots
    tmp(:,ii) = 100*(tmp(:,1) - tmp(:,ii))./tmp(:,1);
    minVal = min(min( tmp(:,ii) ));
    maxVal = max(max( tmp(:,ii) ));
    vals = linspace(0,20,numBins);
    [f(:,ii),~] = hist(tmp(:,ii),vals);
end

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(cad:cad:ST.params.simulation.t_steps);
tmax = max(time);
tmin = min(time);

figure
set(gcf,'name','Magnetic moment','numbertitle','off')
surf(time(2:end),vals,f(:,2:end),'LineStyle','none')
% surf(time,vals,f,'LineStyle','none')
% axis([tmin tmax minVal maxVal])
xlabel('Time (s)','Interpreter','latex','FontSize',16)
ylabel('Magnetic moment $\mu$ (arbitrary units)','Interpreter','latex','FontSize',16)
colormap(jet)
end

function poloidalPlaneDistributions(ST,nbins)
N = 3;
offset = floor(double(ST.params.simulation.num_snapshots)/N);

cad = ST.params.simulation.output_cadence;
time = ST.params.simulation.dt*double(cad:cad:ST.params.simulation.t_steps);

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
time = ST.params.simulation.dt*double(cad:cad:ST.params.simulation.t_steps);

Bo = 1E4*ST.params.fields.Bo;
Ro = 1E2*ST.params.fields.Ro; % Major radius in meters.
a = 1E2*ST.params.fields.a;% Minor radius in meters.
co = 0.5; % Extra parameter
lambda = a/co;
Bpo = 1E4*ST.params.fields.Bpo;

h = figure;
set(h,'name','Angular momentum conservation','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    m = 1E3*ST.params.species.m(ss);
    q = 3E9*ST.params.species.q(ss);
    
    invariant = ...
        zeros(ST.params.species.ppp(ss)*ST.params.simulation.nmpi,ST.params.simulation.num_snapshots);
    err = zeros(1,ST.params.simulation.num_snapshots);
    minerr = zeros(1,ST.params.simulation.num_snapshots);
    maxerr = zeros(1,ST.params.simulation.num_snapshots);
    
    for ii=1:ST.params.simulation.num_snapshots
        X = 1E2*squeeze( ST.data.(['sp' num2str(ss)]).X(:,:,ii) );
        V = 1E2*squeeze( ST.data.(['sp' num2str(ss)]).V(:,:,ii) );

        % Toroidal coordinates
        % r = radius, theta = poloidal angle, phi = toroidal angle
        r = sqrt( (sqrt(X(1,:).^2 + X(2,:).^2) - Ro).^2 + X(3,:).^2 );
        theta = atan2(X(3,:),sqrt(X(1,:).^2 + X(2,:).^2) - Ro);
        theta(theta<0) = theta(theta<0) + 2*pi;
        zeta = atan2(X(1,:),X(2,:));
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        % Toroidal coordinates

        gamma = squeeze( ST.data.(['sp' num2str(ss)]).gamma(:,ii) )';

        eta = r/Ro;
        psi = 0.5*lambda*Bpo*log(1 + r.^2/lambda^2);
        wo = q*Bo./(m*c*gamma);

        dzeta = ...
        (X(2,:).*V(1,:) - X(1,:).*V(2,:))./( sum(X(1:2,:).^2,1) );
    
        invariant(:,ii) = dzeta.*( 1 + eta.*cos(theta) ).^2 - wo.*psi/(Ro*Bo);
        tmp_vec = 100*(invariant(:,1) - invariant(:,ii))./invariant(:,1);
        err(ii) = mean( tmp_vec );
        minerr(ii) = min( tmp_vec );
        maxerr(ii) = max( tmp_vec );
    end
    
    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(time,err,'k-',time,minerr,'r:',time,maxerr,'r:')
    box on
    grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',12)
    ylabel('Angular momentum conservation (\%)','Interpreter','latex','FontSize',12)
end


end

function changeOfMagneticField(ST)

for ss=1:ST.params.simulation.num_species
    m = ST.params.species.m(ss);
    q = ST.params.species.q(ss);
    gamma = ST.data.(['sp' num2str(ss)]).gamma(:,1); % initial relativistic factor
    
end

end