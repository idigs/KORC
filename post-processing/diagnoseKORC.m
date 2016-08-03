function ST = diagnoseKORC(path,visible)
close all

ST = struct;
ST.path = path;
ST.visible = visible;

ST.params = loadSimulationParameters(ST);

ST.time = ...
    ST.params.simulation.dt*double(ST.params.simulation.output_cadence)*...
    double(0:1:ST.params.simulation.num_snapshots);

ST.data = loadData(ST);

energyConservation(ST);

% ST.RT = radialTransport(ST);

% ST.CP = confined_particles(ST);

% ST.PAD = pitchAngleDiagnostic(ST,30);

% ST.MMD = magneticMomentDiagnostic(ST,70);

% poloidalPlaneDistributions(ST,25);

% angularMomentum(ST);

% ST.CMF = changeOfMagneticField(ST)

% energyLimit(ST);

% ST.PR = LarmorVsLL(ST);

% stackedPlots(ST,40);

% scatterPlots(ST);





% save('energy_limit','ST')
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

params.simulation.num_snapshots = 470;
params.simulation.t_steps = params.simulation.output_cadence*params.simulation.num_snapshots;
end

function data = loadData(ST)
data = struct;

list = {'X','V'};%,'B'};

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
list = {'eta','gamma','Prad','Pin','flag','mu'};

for ll=1:length(list)
    disp(['Loading ' list{ll}])
    for ss=1:ST.params.simulation.num_species
        tnp = double(ST.params.species.ppp(ss)*ST.params.simulation.nmpi);
        
        data.(['sp' num2str(ss)]).(list{ll}) = ...
            zeros(tnp,ST.params.simulation.num_snapshots);
        
        for ii=1:(ST.params.simulation.num_snapshots+1)           
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

err = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);
maxerr = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);
minerr = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);

st1 = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);
st2 = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);
st3 = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);
st4 = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);

h1 = figure('Visible',ST.visible);
h2 = figure('Visible',ST.visible);
h3 = figure('Visible',ST.visible);
h4 = figure('Visible',ST.visible);
h5 = figure('Visible',ST.visible);
h6 = figure('Visible',ST.visible);

set(h1,'name','Energy conservation','numbertitle','off')
set(h2,'name','Velocity components','numbertitle','off')
set(h3,'name','Radiated power','numbertitle','off')
set(h4,'name','Input power','numbertitle','off')
set(h5,'name','Energy gain/loss','numbertitle','off')
set(h6,'name','Energy statistics','numbertitle','off')
% try
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
        
        eta = 2*pi*ST.data.(['sp' num2str(ss)]).eta(pin,:)/180;
        V = sqrt(1 - 1./gamma.^2); % in units of the speed of light
        Vpar = mean(V.*cos(eta),1);
        Vperp = mean(V.*sin(eta),1);
        
        figure(h1)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,err(:,ss),'k-',ST.time,minerr(:,ss),'r:',ST.time,maxerr(:,ss),'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\Delta \mathcal{E}/\mathcal{E}_0$ (\%)','Interpreter','latex','FontSize',16)
        saveas(h1,[ST.path 'energy_conservation'],'fig')
        
        figure(h2)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Vpar,'k',ST.time,Vperp,'r')
%         plot(time,Vpar,'k-',time,Vperp,'r-')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$v_\parallel$, $v_\perp$ ($c$)','Interpreter','latex','FontSize',16)
        saveas(h2,[ST.path 'velocity_components'],'fig')

        figure(h6)
        subplot(4,1,1)
        hold on; plot(ST.time,st1(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\mu$','Interpreter','latex','FontSize',16)
        subplot(4,1,2)
        hold on; plot(ST.time,st2(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\sigma$','Interpreter','latex','FontSize',16)
        subplot(4,1,3)
        hold on; plot(ST.time,st3(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$s$','Interpreter','latex','FontSize',16)
        subplot(4,1,4)
        hold on; plot(ST.time,st4(:,ss)); hold off
        box on; grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$k$','Interpreter','latex','FontSize',16)
        saveas(h6,[ST.path 'energy_statistics'],'fig')
        
        
        figure(h3)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Prad,'k',ST.time,minPrad,'r:',ST.time,maxPrad,'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{rad} \rangle$','Interpreter','latex','FontSize',16)
        saveas(h3,[ST.path 'radiated_power'],'fig')
        
        figure(h4)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Pin,'k',ST.time,minPin,'r:',ST.time,maxPin,'r:')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{in} \rangle$','Interpreter','latex','FontSize',16)
        saveas(h4,[ST.path 'input_power'],'fig')
        
        figure(h5)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Prad./Pin,'k')
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\langle P_{rad} \rangle / \langle P_{in} \rangle$','Interpreter','latex','FontSize',16)
        saveas(h5,[ST.path 'energy_gain-lost'],'fig')
    end
% catch
%     error('Something went wrong: energyConservation')
% end


end

function PAD = pitchAngleDiagnostic(ST,nbins)
PAD = struct;
N = 10;

tmax = max(ST.time);
tmin = min(ST.time);

mean_fx = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
std_fx = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
skewness_fx = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
kurtosis_fx = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);

mean_fz = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
std_fz = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
skewness_fz = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
kurtosis_fz = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);

stats = zeros(2,ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);

fx = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots+1);
x = zeros(ST.params.simulation.num_species,nbins);
fz = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots+1);
z = zeros(ST.params.simulation.num_species,nbins);

h1 = figure('Visible',ST.visible);
set(h1,'name','Pitch angle PDF vs. time','numbertitle','off')
h0 = figure('Visible',ST.visible);
set(h0,'name','Energy PDF vs. time','numbertitle','off')
h = figure('Visible',ST.visible);
set(h,'name','Pitch angle: variability','numbertitle','off')
h3 = figure('Visible',ST.visible);
set(h3,'name','PDF pitch angle','numbertitle','off')
hh = figure('Visible',ST.visible);
set(hh,'name','PDF energy','NumberTitle','off')
for ss=1:ST.params.simulation.num_species
    q = ST.params.species.q(ss);
    m = ST.params.species.m(ss);
    c = ST.params.scales.v;
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin & passing;
    eta = ST.data.(['sp' num2str(ss)]).eta(bool,:);
    Eo = ST.data.(['sp' num2str(ss)]).gamma(bool,:)*m*c^2/abs(q);
    Eo = Eo/1E6;
    
    if ~isempty(eta)
        mean_fx(ss,:) = mean(eta,1);
        std_fx(ss,:) = std(eta,0,1);
        skewness_fx(ss,:) = skewness(eta,0,1);
        kurtosis_fx(ss,:) = kurtosis(eta,1,1);
        
        mean_fz(ss,:) = mean(Eo,1);
        std_fz(ss,:) = std(Eo,0,1);
        skewness_fz(ss,:) = skewness(Eo,0,1);
        kurtosis_fz(ss,:) = kurtosis(Eo,1,1);
        
        minVal = min(min( eta ));
        maxVal = max(max( eta ));
        x(ss,:) = linspace(minVal,maxVal,nbins);
        
        minVal = min(min( Eo ));
        maxVal = max(max( Eo ));
        z(ss,:) = linspace(minVal,maxVal,nbins);
        
        for ii=1:ST.params.simulation.num_snapshots+1
            stats(1,ss,ii) = 100*mean((eta(:,ii) - eta(:,1))./eta(:,1));
            stats(2,ss,ii) = 100*std((eta(:,ii) - eta(:,1))./eta(:,1));
            [fx(ss,:,ii),~] = hist(eta(:,ii),x(ss,:));
            [fz(ss,:,ii),~] = hist(Eo(:,ii),z(ss,:));
        end
    end
    
    offset = floor(double(ST.params.simulation.num_snapshots+1)/N);
    
    for ii=1:N
        it = ii*offset;
        nc = floor(N/2);
        nr = floor(N/nc);
        figure(h3)
        subplot(nr,nc,ii)
        hold on
        figure(hh)
        subplot(nr,nc,ii)
        hold on
        
        minVal = min(eta(:,it));
        maxVal = max(eta(:,it));
        xAxis = linspace(minVal,maxVal,nbins);
        f = hist(eta(:,it),xAxis);
        f = f/(mean(diff(xAxis))*sum(f));
        
        figure(h3)
        plot(xAxis,f,'o:')
        
        minVal = min(Eo(:,it));
        maxVal = max(Eo(:,it));
        xAxis = linspace(minVal,maxVal,nbins);
        f = hist(Eo(:,it),xAxis);
        f = f/(mean(diff(xAxis))*sum(f));
        
        figure(hh)
        plot(xAxis,f,'o:')
        
        
        figure(h3)
        title(['Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
        hold off
        box on
        grid on
        figure(hh)
        title(['Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
        hold off
        box on
        grid on
    end
    
    figure(h1)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    surf(ST.time,squeeze(x(ss,:)),log10(squeeze(fx(ss,:,:))),'LineStyle','none')
    %     surf(ST.time,squeeze(x(ss,:)),squeeze(fx(ss,:,:)),'LineStyle','none')
    axis([tmin tmax minVal maxVal])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\eta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    colormap(jet(256))
    
    figure(h)
    subplot(2,1,1)
    hold on
    plot(ST.time,squeeze(stats(1,ss,:)))
    hold off
    xlim([tmin tmax])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('mean($\Delta |\theta|/\theta_0$) ($\%$)','Interpreter','latex','FontSize',16)
    subplot(2,1,2)
    hold on
    plot(ST.time,squeeze(stats(2,ss,:)))
    hold off
    xlim([tmin tmax])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('std($\Delta |\theta|/\theta_0$) ($\%$)','Interpreter','latex','FontSize',16)
    
    
    tmp = ...
        ST.data.(['sp' num2str(ss)]).gamma(bool,:)*ST.params.species.m(ss)*ST.params.scales.v^2/abs(ST.params.species.q(ss));
    tmp = tmp/1E6;
    if ~isempty(tmp)
        minVal = min(min( tmp ));
        maxVal = max(max( tmp ));
        y = linspace(minVal,maxVal,nbins);
        fy = zeros(nbins,ST.params.simulation.num_snapshots+1);
        for ii=1:ST.params.simulation.num_snapshots+1
            [fy(:,ii),~] = hist(tmp(:,ii),y);
        end
    end
    
    figure(h0)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    if ~isempty(tmp)
        surf(ST.time,y,log10(fy),'LineStyle','none')
        axis([tmin tmax minVal maxVal])
    end
    box on
    axis on
    xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
    ylabel('$\mathcal{E}(t)$ (MeV)','Interpreter','latex','FontSize',16)
    colormap(jet(256))
end

h2 = figure('Visible',ST.visible);
set(h2,'name','Statistical moments pitch angle','numbertitle','off')
for ii=1:ST.params.simulation.num_species
    figure(h2)
    subplot(4,1,1)
    hold on
    plot(ST.time,mean_fx(ii,:))
    hold off
    figure(h2)
    subplot(4,1,2)
    hold on
    plot(ST.time,std_fx(ii,:))
    hold off
    figure(h2)
    subplot(4,1,3)
    hold on
    plot(ST.time,skewness_fx(ii,:))
    hold off
    figure(h2)
    subplot(4,1,4)
    hold on
    plot(ST.time,kurtosis_fx(ii,:))
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


% % % Poloidal angle distribution function % %
%
% ft = zeros(ST.params.simulation.num_species,nbins,N);
% t = zeros(ST.params.simulation.num_species,nbins,N);
%
% for ii=1:N
%     it = ii*offset;
%     for ss=1:ST.params.simulation.num_species
% %         I = find( any(ST.data.sp3.eta > 90, 2) == 0 );
% %         X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,I,it));
%         X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,:,it));
%         theta = atan2(X(3,:), sqrt(X(1,:).^2 + X(2,:).^2) - ST.params.fields.Ro);
%         theta(theta < 0) = theta(theta < 0) + 2*pi;
%         theta = (180/pi)*theta;
%         [ft(ss,:,ii),t(ss,:,ii)] = hist(theta,nbins);
%         dt = mean(diff(t(ss,:,ii)));
%         ft(ss,:,ii) = ft(ss,:,ii)/(sum(ft(ss,:,ii))*dt);
%     end
% end
%
% barcolor = [1,0,0;0,1,0;0,0,1;0.5,0.5,1.0;1.0,0.2,0.2];
%
% h4 = figure('Visible',ST.visible);
% set(h4,'name','PDF poloidal angle','numbertitle','off')
% for ii=1:N
%     it = ii*offset;
%     nc = floor(N/2);
%     nr = floor(N/nc);
%     figure(h4)
%     subplot(nr,nc,ii)
%     hold on
%     for ss=1:ST.params.simulation.num_species
%         plot(t(ss,:,ii),ft(ss,:,ii),'o:')
% %         bar(t(ss,:,ii),ft(ss,:,ii),'FaceColor',barcolor(ss,:))
%     end
%     currentAxis = gca;
% %     set(currentAxis.Children(:),'FaceAlpha',0.2);
%     hold off
%     xlim([0 360])
%     title(['Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
%     xlabel('$\theta$','Interpreter','latex','FontSize',16)
%     box on
%     grid on
% end

saveas(h,[ST.path 'variability'],'fig')
saveas(h0,[ST.path 'Energy_PDF_vs_time'],'fig')
saveas(h1,[ST.path 'pitch_vs_time'],'fig')
saveas(h2,[ST.path 'pitch_stats'],'fig')
saveas(h3,[ST.path 'pitch_pdfs'],'fig')
saveas(hh,[ST.path 'energy_pdfs'],'fig')

PAD.mean = mean_fx;
PAD.std = std_fx;

end

function MMD = magneticMomentDiagnostic(ST,numBins)
MMD = struct;
tmax = max(ST.time);
tmin = min(ST.time);

mean_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
std_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
skewness_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
kurtosis_f = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);

fx = zeros(ST.params.simulation.num_species,numBins,ST.params.simulation.num_snapshots+1);
x = zeros(ST.params.simulation.num_species,numBins);

h1 = figure('Visible',ST.visible);
set(h1,'name','PDF of magnetic moment','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    tmp = ST.data.(['sp' num2str(ss)]).mu(pin,:);
    
    if ~isempty(tmp)
        for ii=1:size(tmp,1)
            tmp(ii,:) = 100*abs(tmp(ii,:) - tmp(ii,1))./tmp(ii,1);
        end
        
        mean_f(ss,:) = mean(tmp,1);
        std_f(ss,:) = std(tmp,0,1);
        skewness_f(ss,:) = skewness(tmp,0,1);
        kurtosis_f(ss,:) = kurtosis(tmp,1,1);
        
        minVal = min(min( tmp ));
        maxVal = max(max( tmp ));
        x(ss,:) = linspace(minVal,maxVal,numBins);
        
        for ii=1:ST.params.simulation.num_snapshots+1
            try
                [fx(ss,:,ii),~] = hist(tmp(:,ii),x(ss,:));
%                 fx(ss,:,ii) = fx(ss,:,ii)/( mean(diff(x(ss,:)))*sum(fx(ss,:,ii)) );
                fx(ss,:,ii) = fx(ss,:,ii)/max(squeeze(fx(ss,:,ii)));
            catch
            end
        end
    end
    
    subplot(double(ST.params.simulation.num_species),1,double(ss))
%     surf(ST.time,squeeze(x(ss,:)),log10(squeeze(fx(ss,:,2:end))),'LineStyle','none')
    surf(ST.time(2:end),squeeze(x(ss,:)),squeeze(fx(ss,:,2:end)),'LineStyle','none')
%     axis([tmin tmax minVal maxVal])
    view([0,90])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\Delta \mu/\mu_0$ ($\%$)','Interpreter','latex','FontSize',16)
    colormap(jet)
end

h2 = figure('Visible',ST.visible);
set(h2,'name','Statistical moments magnetic moment','numbertitle','off')
for ii=1:ST.params.simulation.num_species
    figure(h2)
    subplot(4,1,1)
    hold on
    plot(ST.time,mean_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,2)
    hold on
    plot(ST.time,std_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,3)
    hold on
    plot(ST.time,skewness_f(ii,:))
    hold off
    figure(h2)
    subplot(4,1,4)
    hold on
    plot(ST.time,kurtosis_f(ii,:))
    hold off
end

saveas(h1,[ST.path 'magnetic_moment_pdfs'],'fig')
saveas(h2,[ST.path 'magnetic_moment_stats'],'fig')

MMD.stat1 = mean_f;
MMD.stat2 = std_f;
MMD.stat3 = skewness_f;
MMD.stat4 = kurtosis_f;
end

function poloidalPlaneDistributions(ST,nbins)
N = 3;
offset = floor(double(ST.params.simulation.num_snapshots+1)/N);

% cad = ST.params.simulation.output_cadence;
% time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);

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
        title(['Species: ' num2str(ss) 'Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
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
%     title(['Species: ' num2str(ss) 'Time: ' num2str(ST.time(it))],'Interpreter','latex','FontSize',11)
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

% cad = ST.params.simulation.output_cadence;
% time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);

Bo = 1E4*ST.params.fields.Bo;
Ro = 1E2*ST.params.fields.Ro; % Major radius in meters.
a = 1E2*ST.params.fields.a;% Minor radius in meters.
lambda = a/sqrt(ST.params.fields.qa/ST.params.fields.qo - 1.0);
Bpo = 1E4*ST.params.fields.Bpo;

st1 = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);
st2 = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);
st3 = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);
st4 = zeros(ST.params.simulation.num_snapshots+1,ST.params.simulation.num_species);

h = figure('Visible',ST.visible);
set(h,'name','Angular momentum conservation','numbertitle','off')
h1 = figure('Visible',ST.visible);
set(h1,'name','Angular momentum statistics','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    num_part = numel(find(pin==1));
%     num_part = ST.params.species.ppp(ss)*ST.params.simulation.nmpi;
    
    m = 1E3*ST.params.species.m(ss);
    q = 3E9*ST.params.species.q(ss);
    
    invariant = ...
        zeros(num_part,ST.params.simulation.num_snapshots+1);
    err = zeros(1,ST.params.simulation.num_snapshots+1);
    minerr = zeros(1,ST.params.simulation.num_snapshots+1);
    maxerr = zeros(1,ST.params.simulation.num_snapshots+1);
    
    for ii=1:ST.params.simulation.num_snapshots+1
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
    plot(ST.time,err,'k-',ST.time,minerr,'r:',ST.time,maxerr,'r:')
    box on
    grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',12)
    ylabel('$\Delta p_\zeta/p_{\zeta 0}$ (\%)','Interpreter','latex','FontSize',12)
    saveas(h,[ST.path 'ang_mom_conservation'],'fig')
    
    figure(h1)
    subplot(4,1,1)
    hold on; plot(ST.time,st1(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\mu$','Interpreter','latex','FontSize',16)
    subplot(4,1,2)
    hold on; plot(ST.time,st2(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\sigma$','Interpreter','latex','FontSize',16)
    subplot(4,1,3)
    hold on; plot(ST.time,st3(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$s$','Interpreter','latex','FontSize',16)
    subplot(4,1,4)
    hold on; plot(ST.time,st4(:,ss)); hold off
    box on; grid on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$k$','Interpreter','latex','FontSize',16)
    saveas(h1,[ST.path 'ang_mom_stat'],'fig')
end


end

function CMF = changeOfMagneticField(ST)
CMF = struct;

CMF.stat.mag = zeros(2,ST.params.simulation.num_species);
CMF.stat.vec = zeros(2,ST.params.simulation.num_species);

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
    DBvec = zeros(size(Bo));
    DBmag = zeros(size(Bo));
    for ii=1:numel(Bo)
        [~,I(ii)] = min( abs(Tc(ii) - ST.time) );
        tmp = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),1:I(ii)).^2,1) ) );
        DBmag(ii) = max( 100*abs(tmp - Bo(ii))/Bo(ii) );
        
        B1 = squeeze(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),1));
        for jj=1:I(ii)
            B2 = squeeze(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),jj));
            tmp = 100*sqrt(sum((B2-B1).^2))/sqrt(sum(B1.^2));
            if(abs(tmp) > DBvec(ii))
                DBvec(ii) = tmp;
            end
        end
        
%         B(ii) = squeeze( sqrt( sum(ST.data.(['sp' num2str(ss)]).B(:,aux(ii),I(ii)).^2,1) ) );
    end
%     DBmag = 100*abs( B - Bo )./Bo;

    CMF.stat.mag(1,ss) = mean(DBmag);
    CMF.stat.mag(2,ss) = std(DBmag);
    CMF.stat.vec(1,ss) = mean(DBvec);
    CMF.stat.vec(2,ss) = std(DBvec);

    xmag = linspace(min(DBmag),max(DBmag),35);
    xvec = linspace(min(DBvec),max(DBvec),35);
    xAxisMin = min([min(DBmag) min(DBvec)]);
    xAxisMax = max([max(DBmag) max(DBvec)]);

    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,pin,1));
    R = sqrt( sum(X(1:2,:).^2,1) );
    Z = X(3,:);

    S = 12*ones(size(gamma));
    h=figure('Visible',ST.visible,'units','normalized','OuterPosition',[0.1,0.25,0.75,0.4]);
    set(h,'name',['Change of B-field: sp' num2str(ss)],'numbertitle','off')
    subplot(1,3,1)
    histogram(DBvec,xvec)
    hold on
    histogram(DBmag,xmag)
    hold off
    xlim([xAxisMin xAxisMax])
    currentAxis = gca;
    set(currentAxis.Children(:),'FaceAlpha',0.2);
    ylabel('$f(\Delta B/B_0)$','Interpreter','latex','FontSize',16)
    xlabel('$\Delta B/B_0$ ($\%$)','Interpreter','latex','FontSize',16)
    subplot(1,3,2)
    scatter3(R,Z,DBvec,S,DBvec,'square','filled')
    title('$\Delta \mathbf{B}$','Interpreter','latex','FontSize',16)
    colormap(jet(256))
    colorbar
    view([0,90])
    box on; axis on; axis square
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    subplot(1,3,3)
    scatter3(R,Z,DBmag,S,DBmag,'square','filled')
    title('$\Delta B$','Interpreter','latex','FontSize',16)
    colormap(jet(256))
    colorbar
    view([0,90])
    box on; axis on; axis square
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    saveas(h,[ST.path 'B-field_sp' num2str(ss)],'fig')
%     close(h)
end



end

function RT = radialTransport(ST)
RT = struct;

% cad = ST.params.simulation.output_cadence;
% time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
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

function CP = confined_particles(ST)
CP = struct;

tmax = max(ST.time);
tmin = min(ST.time);

t = linspace(0,2*pi,200);
Rs = ST.params.fields.Ro + ST.params.fields.a*cos(t);
Zs = ST.params.fields.a*sin(t);

CP.confined = zeros(1,ST.params.simulation.num_species);

h0 = figure('Visible',ST.visible);
set(h0,'name','Particle loss','numbertitle','off');

C = colormap(h0,jet(512));
offset = floor(512/ST.params.simulation.num_species);
colour = C(1:offset:end,:);

legends = cell(1,ST.params.simulation.num_species);
for ss=1:ST.params.simulation.num_species
    numConfPart = sum(ST.data.(['sp' num2str(ss)]).flag(:,1),1);
    confinedParticles = ...
        100*sum(ST.data.(['sp' num2str(ss)]).flag,1)/numConfPart;
    figure(h0)
    hold on
    plot(ST.time,confinedParticles)
    hold off
    
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
    
    CP.confined(ss) = confinedParticles(end);
end


figure(h0)
legend(legends,'interpreter','latex','FontSize',12)
box on
axis on
axis square
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
ylabel('$\%$ of confined RE','Interpreter','latex','FontSize',16)
saveas(h0,[ST.path 'particle_loss'],'fig')
% close(h0)

h1=figure;
set(h1,'Visible',ST.visible,'name','IC','numbertitle','off')
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
saveas(h1,[ST.path 'IC'],'fig')
% close(h1)
end

function energyLimit(ST)
% cad = ST.params.simulation.output_cadence;
% time = ST.params.simulation.dt*double(0:cad:ST.params.simulation.t_steps);
tmax = max(ST.time);
tmin = min(ST.time);

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
    plot(ST.time,energy)
    hold off
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
legend(legends,'interpreter','latex','FontSize',12)
grid on; box on; axis on; axis square
xlabel('Pitch angle $\eta$ ($^\circ$)','Interpreter','latex','FontSize',16)
ylabel('$\mathcal{E}_0$ (MeV)','Interpreter','latex','FontSize',16)
saveas(h1,[ST.path 'Energy_vs_pitch'],'fig')


figure(h2)
legend(legends,'interpreter','latex','FontSize',12)
grid on; box on; axis on; axis square
xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
ylabel('$\mathcal{E}_0$ (MeV)','Interpreter','latex','FontSize',16)
saveas(h2,[ST.path 'mean_energy_vs_time'],'fig')
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

function PR = LarmorVsLL(ST)
PR = struct;
tmax = max(ST.time);
tmin = min(ST.time);

kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
% % % % % % % % % % %

PR.A = zeros(3,ST.params.simulation.num_species);

h = figure('Visible',ST.visible);
set(h,'name','Radiation: Larmor approx','numbertitle','off')
h1 = figure('Visible',ST.visible);
set(h1,'name','Scatter plot','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    aux = find(pin == 1);
    
    V = ST.data.(['sp' num2str(ss)]).V(:,pin,:);
    v = squeeze( sqrt( sum(V.^2,1) ) );
    gamma = ST.data.(['sp' num2str(ss)]).gamma(pin,:);
    eta = pi*ST.data.(['sp' num2str(ss)]).eta(pin,:)/180;
%     eta = pi*repmat(ST.params.species.etao(ss),size(gamma,1),size(gamma,2))/180;
    X = ST.data.(['sp' num2str(ss)]).X(:,pin,:);
    try
        B = ST.data.(['sp' num2str(ss)]).B(:,pin,:);
    catch
        B = analyticalB(ST,X);
    end
    E = analyticalE(ST,X);
    
    vec_mag = zeros(size(gamma));
    for it=1:size(X,3)
        for ii=1:size(gamma,1)
            VxE = cross(squeeze(V(:,ii,it)),squeeze(E(:,ii,it)));
            VxB = cross(squeeze(V(:,ii,it)),squeeze(B(:,ii,it)));
            VxVxB = cross(squeeze(V(:,ii,it)),VxB);
            vec = VxE + VxVxB;
            vec_mag(ii,it) = sqrt( vec'*vec );
        end
    end

    % Approximation of <1/R^2> ~ sin^4(eta)/rg^2
    vperp = v.*sin(eta);
    wc = q*ST.params.fields.Bo./(gamma*m);
    rg = vperp./wc;
    kappa2 = sin(eta).^4./rg.^2;
    
    % Actual curvature
    kappa = q*vec_mag./(m*gamma.*v.^3);    

    % Landau-Lifshiftz radiation formula
%     PR_LL = mean(abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:)),1);
    PR_LL = abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:));
    
    % Larmor approximation using actual curvature
%     PR_L = 2*Kc*q^2*mean((gamma.*v).^4.*kappa.^2,1)/(3*c^3);
    PR_L = 2*Kc*q^2*(gamma.*v).^4.*kappa.^2/(3*c^3);
    
    % Larmor approximation using approximation for curvature
%     PR_approx = 2*Kc*q^2*mean((gamma.*v).^4.*kappa2,1)/(3*c^3);
    PR_approx = 2*Kc*q^2*(gamma.*v).^4.*kappa2/(3*c^3);
    
    figure(h1)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(squeeze(180*eta(:,end)/pi),squeeze(PR_LL(:,end)),'b.',...
        squeeze(180*eta(:,end)/pi),squeeze(PR_approx(:,end)),'r.')
    box on; grid on;
    xlabel('$\eta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    ylabel('$P_{rad}$','Interpreter','latex','FontSize',16)
    saveas(h1,[ST.path 'radiation_scatter_plot'],'fig')
    
    PR_LL = mean(PR_LL,1);
    PR_L = mean(PR_L,1);
    PR_approx = mean(PR_approx,1);
    
    RATIO1 = PR_L./PR_LL;
    RATIO2 = PR_approx./PR_LL;   
    PR.A(:,ss) = [PR_LL(end), PR_L(end), PR_approx(end)];
    
    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(ST.time,RATIO1,'k',ST.time,RATIO2,'r')
    box on; grid on;
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$P_R$','Interpreter','latex','FontSize',12)
    saveas(h,[ST.path 'radiation_ratios'],'fig')
end

end

function stackedPlots(ST,nbins)
tmax = max(ST.time);
tmin = min(ST.time);

fx = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots+1);
Prad = zeros(ST.params.simulation.num_species,ST.params.simulation.num_snapshots+1);
x = zeros(ST.params.simulation.num_species,nbins);

h2 = figure;
set(h2,'name','PDF of radiated power','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
%     data = ST.data.(['sp' num2str(ss)]).eta(pin,:);
    data = abs(ST.data.(['sp' num2str(ss)]).Prad(pin,:));

    minVal = min(min( data ));
    maxVal = max(max( data ));
    x(ss,:) = linspace(minVal,maxVal,nbins);
    
    for ii=1:ST.params.simulation.num_snapshots+1
        [fx(ss,:,ii),~] = hist(data(:,ii),x(ss,:));
        Prad(ss,ii) = mean(diff(x(ss,:)))*sum(x(ss,:).*fx(ss,:,ii));
    end
    
    figure(h2)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
%     surf(ST.time,squeeze(x(ss,:)),log10(squeeze(fx(ss,:,:))),'LineStyle','none')
    plot(ST.time,Prad(ss,:))
    axis([tmin tmax min(Prad(ss,:)) max(Prad(ss,:))])
%     axis([tmin tmax minVal maxVal])
    box on
    axis on
    xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
    ylabel('$P_{rad}$ (Watts/electron)','Interpreter','latex','FontSize',16)
    colormap(jet)
    colorbar
end


end

function scatterPlots(ST)
kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
% % % % % % % % % % %

h=figure;
set(h,'name','Scatter: E vs pitch','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    eta = ST.data.(['sp' num2str(ss)]).eta(pin,:);
    Et = ST.data.(['sp' num2str(ss)]).gamma(pin,:)*m*ST.params.scales.v^2/(abs(q)*1E6);
    Eo = ST.params.species.Eo(ss)/1E6;
    etao = ST.params.species.etao(ss);

    eta_mean = mean(eta,1);
    Et_mean = mean(Et,1);

    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(squeeze(eta(:,end)),squeeze(Et(:,end)),'.','MarkerFaceColor',[0,0.45,0.74],...
        'MarkerEdgeColor',[0,0.45,0.74])
    hold on
    plot(eta_mean,Et_mean,'k',...
        etao,Eo,'rs','MarkerFaceColor','r','MarkerSize',10)
    hold off
    box on; grid on;
    xlabel('$\eta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
end

end


