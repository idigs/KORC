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

% energyConservation(ST);

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

% ST.P = synchrotronSpectrum(ST,true);

% ST.VS = identifyVisibleParticles(ST);

ST.SD = syntheticDiagnosticSynchrotron(ST);

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

% params.simulation.num_snapshots = 8;
% params.simulation.t_steps = params.simulation.output_cadence*params.simulation.num_snapshots;
end

function data = loadData(ST)
data = struct;

if isfield(ST.params.fields,'dims')
    list = {'X','V','B'};
else
    list = {'X','V','B'};
end

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


list = {'eta','gamma','Prad','Pin','flag','mu'};
% list = {'eta','gamma','Prad','flag'};

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
        gammap = ST.data.(['sp' num2str(ss)]).gamma(pin,:);
        tmp = zeros(size(gammap));
        for ii=1:size(tmp,1)
            tmp(ii,:) = ...
                ( gammap(ii,:) - gammap(ii,1) )./gammap(ii,1);
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
        V = sqrt(1 - 1./gammap.^2); % in units of the speed of light
        Vpar = mean((V.*cos(eta)).^2,1);
        Vperp = mean((V.*sin(eta)).^2,1);
        
        figure(h1)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        try
            plot(ST.time,err(:,ss),'k-',ST.time,minerr(:,ss),'r:',ST.time,maxerr(:,ss),'r:')
%             hold on
%             plot(ST.time,tmp)
%             hold off
        catch
        end
        box on
        grid on
        xlabel('Time (s)','Interpreter','latex','FontSize',16)
        ylabel('$\Delta \mathcal{E}/\mathcal{E}_0$ (\%)','Interpreter','latex','FontSize',16)
        saveas(h1,[ST.path 'energy_conservation'],'fig')
        
        figure(h2)
        subplot(double(ST.params.simulation.num_species),1,double(ss))
        plot(ST.time,Vpar-Vperp,'k-')
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
x = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots+1);
fz = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots+1);
z = zeros(ST.params.simulation.num_species,nbins,ST.params.simulation.num_snapshots+1);

if ST.visible
    h1 = figure('Visible','on');
else
    h1 = figure('Visible','off');
end

set(h1,'name','Pitch angle PDF vs. time','numbertitle','off')
% h0 = figure('Visible',ST.visible);
% set(h0,'name','Energy PDF vs. time','numbertitle','off')
if ST.visible
    h = figure('Visible','on');
else
    h = figure('Visible','off');
end
set(h,'name','Pitch angle: variability','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    c = ST.params.scales.v;
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin & passing;
    eta = ST.data.(['sp' num2str(ss)]).eta(bool,:);
    Eo = ST.data.(['sp' num2str(ss)]).gamma(bool,:)*m*c^2/q;
    Eo = Eo/1E6;
    
%     try
%         tmp = eta(:,end);
%         save(['pitch_' num2str(ST.params.species.etao(ss)) '_E_' num2str(floor(Eo(1,1)))],'tmp')
%     catch
%     end
    
    if ~isempty(eta)
        mean_fx(ss,:) = mean(eta,1);
        std_fx(ss,:) = std(eta,0,1);
        skewness_fx(ss,:) = skewness(eta,0,1);
        kurtosis_fx(ss,:) = kurtosis(eta,1,1);
        
        mean_fz(ss,:) = mean(Eo,1);
        std_fz(ss,:) = std(Eo,0,1);
        skewness_fz(ss,:) = skewness(Eo,0,1);
        kurtosis_fz(ss,:) = kurtosis(Eo,1,1);
        
        for ii=1:ST.params.simulation.num_snapshots+1
            minVal = min(min( eta(:,ii) ));
            maxVal = max(max( eta(:,ii) ));
            x(ss,:,ii) = linspace(minVal,maxVal,nbins);
            Dx = mean(diff(x(ss,:,ii)));
            
            minVal = min(min( Eo(:,ii) ));
            maxVal = max(max( Eo(:,ii) ));
            z(ss,:,ii) = linspace(minVal,maxVal,nbins);
            Dz = mean(diff(z(ss,:,ii)));
            
            stats(1,ss,ii) = 100*mean(abs(eta(:,ii) - eta(:,1))./eta(:,1));
            stats(2,ss,ii) = 100*std(abs(eta(:,ii) - eta(:,1))./eta(:,1));
            [fx(ss,:,ii),~] = hist(eta(:,ii),x(ss,:,ii));
            fx(ss,:,ii) = fx(ss,:,ii)/(Dx*sum(fx(ss,:,ii)));
            [fz(ss,:,ii),~] = hist(Eo(:,ii),z(ss,:,ii));
            fz(ss,:,ii) = fz(ss,:,ii)/(Dz*sum(fz(ss,:,ii)));
        end
    end
    
    figure(h1)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    surf(ST.time,squeeze(x(ss,:,:)),log10(squeeze(fx(ss,:,:))),'LineStyle','none')
%     surf(ST.time,squeeze(x(ss,:)),squeeze(fx(ss,:,:)),'LineStyle','none')
%     axis([tmin tmax minVal maxVal])
    box on
    axis on
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
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
    
    
%     tmp = ...
%         ST.data.(['sp' num2str(ss)]).gamma(bool,:)*ST.params.species.m(ss)*ST.params.scales.v^2/abs(ST.params.species.q(ss));
%     tmp = tmp/1E6;
%     if ~isempty(tmp)
%         minVal = min(min( tmp ));
%         maxVal = max(max( tmp ));
%         y = linspace(minVal,maxVal,nbins);
%         fy = zeros(nbins,ST.params.simulation.num_snapshots+1);
%         for ii=1:ST.params.simulation.num_snapshots+1
%             [fy(:,ii),~] = hist(tmp(:,ii),y);
%         end
%     end
%     
%     figure(h0)
%     subplot(double(ST.params.simulation.num_species),1,double(ss))
%     if ~isempty(tmp)
%         surf(ST.time,y,log10(fy),'LineStyle','none')
%         axis([tmin tmax minVal maxVal])
%     end
%     box on
%     axis on
%     xlabel('Time $t$ (sec)','Interpreter','latex','FontSize',16)
%     ylabel('$\mathcal{E}(t)$ (MeV)','Interpreter','latex','FontSize',16)
%     colormap(jet(256))
end

if ST.visible
    h2 = figure('Visible','on');
else
    h2 = figure('Visible','off');
end
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

offset = floor(double(ST.params.simulation.num_snapshots+1)/N);

% z = linspace(-4,4,100);
% fz = exp( -0.5*z.^2 )/sqrt(2*pi);

if ST.visible
    h3 = figure('Visible','on');
else
    h3 = figure('Visible','off');
end
set(h3,'name','PDF pitch angle','numbertitle','off')
if ST.visible
    hh = figure('Visible','on');
else
    hh = figure('Visible','off');
end
set(hh,'name','PDF energy','NumberTitle','off')
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
%     plot(z,log10(fz),'k')
    
    for ss=1:ST.params.simulation.num_species    
        xAxis = ( x(ss,:,it) - mean(x(ss,:,it)) )/std(fx(ss,:,it));
        f = std(fx(ss,:,it))*squeeze(fx(ss,:,it))/trapz(x(ss,:,it),fx(ss,:,it));
        
        figure(h3)
%         plot(xAxis,f,'o:')
        plot(xAxis,log10(f),'o:')
                
        xAxis = z(ss,:,it);
        f = squeeze(fz(ss,:,it));
        
        figure(hh)
        plot(xAxis,f,'o:')
    end
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
% saveas(h0,[ST.path 'Energy_PDF_vs_time'],'fig')
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

if ST.visible
    h1 = figure('Visible','on');
else
    h1 = figure('Visible','off');
end
set(h1,'name','PDF of magnetic moment','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    tmp = ST.data.(['sp' num2str(ss)]).mu(pin,:);
    
    if ~isempty(tmp)
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
        
        for ii=1:ST.params.simulation.num_snapshots+1
            try
                [fx(ss,:,ii),~] = hist(tmp(:,ii),x(ss,:));
                fx(ss,:,ii) = fx(ss,:,ii)/( trapz(x(ss,:),fx(ss,:,ii)) );
%                 fx(ss,:,ii) = fx(ss,:,ii)/max(squeeze(fx(ss,:,ii)));
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

if ST.visible
    h2 = figure('Visible','on');
else
    h2 = figure('Visible','off');
end
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

MMD.fx = fx;
MMD.x = x;
MMD.stat1 = mean_f;
MMD.stat2 = std_f;
MMD.stat3 = skewness_f;
MMD.stat4 = kurtosis_f;
end

function poloidalPlaneDistributions(ST,nbins)
N = 3;
offset = floor(double(ST.params.simulation.num_snapshots+1)/N);

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
Bo = ST.params.fields.Bo;
Ro = ST.params.fields.Ro; % Major radius in meters.
lambda = ST.params.fields.lambda;
qo = ST.params.fields.qo;

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
    
    m = ST.params.species.m(ss);
    q = ST.params.species.q(ss);
    
    I = zeros(num_part,ST.params.simulation.num_snapshots+1);
    err = zeros(1,ST.params.simulation.num_snapshots+1);
    minerr = zeros(1,ST.params.simulation.num_snapshots+1);
    maxerr = zeros(1,ST.params.simulation.num_snapshots+1);
    
    for ii=1:ST.params.simulation.num_snapshots+1
        X = squeeze( ST.data.(['sp' num2str(ss)]).X(:,pin,ii) );
        V = squeeze( ST.data.(['sp' num2str(ss)]).V(:,pin,ii) );
        gammap = squeeze( ST.data.(['sp' num2str(ss)]).gamma(pin,ii) )';
        
        
        % Toroidal coordinates
        % r = radius, theta = poloidal angle, phi = toroidal angle
        r = sqrt( (sqrt(X(1,:).^2 + X(2,:).^2) - Ro).^2 + X(3,:).^2 );
        theta = atan2(X(3,:),sqrt(X(1,:).^2 + X(2,:).^2) - Ro);
        theta(theta<0) = theta(theta<0) + 2*pi;
        zeta = atan2(X(1,:),X(2,:));
        zeta(zeta<0) = zeta(zeta<0) + 2*pi;
        % Toroidal coordinates

        eta = r/Ro;
        wc = q*Bo./(m*gammap);

        dzeta_dt = (X(2,:).*V(1,:) - X(1,:).*V(2,:))./( X(1,:).^2 + X(2,:).^2 );
    
        T1 = dzeta_dt.*(1 + eta.*cos(theta)).^2;
        T2 = lambda^2*wc.*log(1 + (r/lambda).^2)/(2*qo*Ro^2);
        I(:,ii) = T1 - T2;
        
        tmp_vec = (I(:,ii) - I(:,1))./I(:,1);
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
    
    for pp=1:size(I,1)
        I(pp,:) = (I(pp,:) - I(pp,1))/I(pp,1);
    end
    
    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    try
        plot(ST.time,err,'k-',ST.time,minerr,'r:',ST.time,maxerr,'r:')
        hold on
        plot(ST.time,I)
        hold off
    catch
    end
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
    gammap = ST.data.(['sp' num2str(ss)]).gamma(pin,1)';
    wc = abs(q)*Bo./(gammap*m);
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

    S = 12*ones(size(gammap));
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

if isfield(ST.params.fields,'a')
    t = linspace(0,2*pi,200);
    Rs = ST.params.fields.Ro + ST.params.fields.a*cos(t);
    Zs = ST.params.fields.a*sin(t);
end

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
close(h0)

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
    Prad = squeeze(abs(ST.data.(['sp' num2str(ss)]).Prad(bool,end)));

    figure(h1)
    subplot(1,2,1)
    hold on
    plot(R,Z,'s','MarkerSize',4,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
%     plot(R,Z,'.','MarkerSize',10,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
%     plot3(R,Z,Prad,'s','MarkerSize',4,'MarkerFaceColor',colour(ss,:),'MarkerEdgeColor',colour(ss,:))
    hold off
    legends{ss} = ['$\eta_0 =$' num2str(ST.params.species.etao(ss)) '$^\circ$'];
end

figure(h1)
subplot(1,2,1)
legend(legends,'interpreter','latex','FontSize',12)
if isfield(ST.params.fields,'a')
    hold on
    plot(Rs,Zs,'r')
    hold off
end
box on
axis on
axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)

legends = cell(1,ST.params.simulation.num_species);
for ss=ST.params.simulation.num_species:-1:1
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
if isfield(ST.params.fields,'a')
    hold on
    plot(Rs,Zs,'r')
    box on
end
axis on
axis square
xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
saveas(h1,[ST.path 'IC'],'fig')
close(h1)
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
lamb = ST.params.fields.lambda;
try
    co = ST.params.fields.co;
catch
    qo = ST.params.fields.qa;
    co = a/lamb;
end
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

kB = 1.38E-23; % Boltzmann constant
Kc = 8.987E9; % Coulomb constant in N*m^2/C^2
mu0 = (4E-7)*pi; % Magnetic permeability
ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass
% % % % % % % % % % %

PR.models = zeros(3,ST.params.simulation.num_species);

h = figure('Visible',ST.visible);
set(h,'name','Model comparison: ratios','numbertitle','off')
h1 = figure('Visible',ST.visible);
set(h1,'name','Scatter plot: comparison','numbertitle','off')
h2 = figure('Visible',ST.visible);
set(h2,'name','Scatter plot: actual radiation','numbertitle','off')
h3 = figure('Visible',ST.visible);
set(h3,'name','Radiation bar plots','numbertitle','off')
for ss=1:ST.params.simulation.num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) );
    bool = pin & passing;
    aux = find(bool == 1);
    S = numel(aux);
    
    V = ST.data.(['sp' num2str(ss)]).V(:,bool,:);
    v = squeeze( sqrt( sum(V.^2,1) ) );
    gammap = ST.data.(['sp' num2str(ss)]).gamma(bool,:);
    eta = pi*ST.data.(['sp' num2str(ss)]).eta(bool,:)/180;
    
    gammapo = repmat(ST.params.species.gammao(ss),S,ST.params.simulation.num_snapshots+1);
    vo = ST.params.scales.v*sqrt(1 - 1./gammapo.^2);
    etao = pi*repmat(ST.params.species.etao(ss),S,ST.params.simulation.num_snapshots+1)/180;
    
%     gammapo = gammap;
%     etao = eta;
%     vo = v;
    
    X = ST.data.(['sp' num2str(ss)]).X(:,bool,:);
    try
        B = ST.data.(['sp' num2str(ss)]).B(:,bool,:);
    catch
        B = analyticalB(ST,X);
    end
    E = analyticalE(ST,X);
    
    vec_mag = zeros(size(gammap));
    for it=1:size(X,3)
        for ii=1:size(gammap,1)
            VxE = cross(squeeze(V(:,ii,it)),squeeze(E(:,ii,it)));
            VxB = cross(squeeze(V(:,ii,it)),squeeze(B(:,ii,it)));
            VxVxB = cross(squeeze(V(:,ii,it)),VxB);
            vec = VxE + VxVxB;
            vec_mag(ii,it) = sqrt( vec'*vec );
        end
    end

    % Approximation of <1/R^2> ~ sin^4(eta)/rg^2
    vperp = vo.*sin(etao);
    wc = q*ST.params.fields.Bo./(gammapo*m);
    rg = vperp./wc;
    kappa2 = sin(etao).^4./rg.^2;
 
    % Actual curvature
    kappa = q*vec_mag./(m*gammap.*v.^3);    

    % Landau-Lifshiftz radiation formula
    PR_LL = abs(ST.data.(['sp' num2str(ss)]).Prad(bool,:));
    
    % Larmor approximation using actual curvature
%     PR_L = 2*Kc*q^2*(gammap.*v).^4.*kappa.^2/(3*c^3);
    Tr = 6*pi*ep*(m*ST.params.scales.v)^3./(q^4*squeeze(sum(B.^2,1)));
    PR_L = gammap.*v.*(m*gammap.*v).*sin(eta).^2./Tr;
    
    % Larmor approximation using approximation for curvature
%     PR_app = 2*Kc*q^2*(gammapo.*v).^4.*kappa2/(3*c^3);
    Tr = 6*pi*ep*(m*ST.params.scales.v)^3/(q^4*ST.params.fields.Bo^2);
    PR_app = gammapo.*vo.*(m*gammapo.*vo).*sin(etao).^2/Tr;
    
    figure(h1)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    hold on
    plot(squeeze(180*eta(:,end)/pi),squeeze(PR_L(:,end)),'b.','MarkerSize',2)
    plot(squeeze(180*etao(:,end)/pi),squeeze(PR_app(:,end)),'r.','MarkerSize',2)
%     plot(squeeze(180*eta(:,end)/pi),squeeze(PR_app(:,end)),'r.','MarkerSize',2)
    hold off
    box on; grid on;
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    ylabel('$P_R$','Interpreter','latex','FontSize',16)
    
    figure(h2)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(squeeze(180*eta(:,end)/pi),squeeze(PR_LL(:,end)),'k.','MarkerSize',4)
    box on; grid on;
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    ylabel('$P_R$','Interpreter','latex','FontSize',16)
    
    try
        minVal = min( eta(:,end) );
        maxVal = max( eta(:,end) );
        x = linspace(minVal,maxVal,20);
        Dx = mean(diff(x));
        
        [fx,~] = hist(eta(:,end),x);
%         fx = fx/(Dx*sum(fx));
    catch
    end
    
    
    try
        minVal = min( PR_LL(:,end) );
        maxVal = max( PR_LL(:,end) );
        z = linspace(minVal,maxVal,20);
        Dz = mean(diff(z));
        
        [fz,~] = hist(PR_LL(:,end),z);
%         fz = fz/(Dz*sum(fz));
    catch
    end
    
    figure(h3)
    offset = 2*(double(ss) - 1);
    subplot(double(ST.params.simulation.num_species),2,offset + 1)
    bar(180*x/pi,fx)
    xlabel('$\theta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    subplot(double(ST.params.simulation.num_species),2,offset + 2)
    bar(z,fz)
    xlabel('$P_R$ (Watts/electron)','Interpreter','latex','FontSize',16)
    
%     PR_LL = mean(PR_LL,1);
    PR_LL = sum(PR_LL,1);
    PR_L = mean(PR_L,1);
    PR_app = mean(PR_app,1);
    
    RATIO1 = PR_L./PR_LL;
    RATIO2 = PR_app./PR_LL;   
    PR.models(:,ss) = [PR_LL(end), PR_L(end), PR_app(end)];
    
    figure(h)
    subplot(double(ST.params.simulation.num_species),1,double(ss))
    plot(ST.time,RATIO1,'k',ST.time,RATIO2,'r')
    box on; grid on;
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$P_R$','Interpreter','latex','FontSize',12)
end
saveas(h,[ST.path 'radiation_ratios'],'fig')
saveas(h1,[ST.path 'scatter_plot_comparison'],'fig')
saveas(h2,[ST.path 'scatter_plot_actual_radiation'],'fig')
saveas(h3,[ST.path 'radiation_pdfs'],'fig')
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
%     hold on
%     plot(eta_mean,Et_mean,'k',...
%         etao,Eo,'rs','MarkerFaceColor','r','MarkerSize',10)
%     hold off
    box on; grid on;
    xlabel('$\eta$ ($^\circ$)','Interpreter','latex','FontSize',16)
    ylabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)
end
saveas(h,[ST.path 'scatter_plot_E_vs_pitch'],'fig')

end

function P = synchrotronSpectrum(ST,filtered)
disp('Calculating spectrum of synchrotron radiation...')
P = struct;

numSnapshots = 3;
it1 = ST.params.simulation.num_snapshots + 1 - numSnapshots;
it2 = ST.params.simulation.num_snapshots + 1;

geometry = 'cylindrical';

upper_integration_limit = 200.0;

N = 500;
lambda_min = 450E-9;% in meters
lambda_max = 950E-9;% in meters
% lambda_min = 907E-9;% in meters
% lambda_max = 917E-9;% in meters
% lambda_min = 742E-9;% in meters
% lambda_max = 752E-9;% in meters

% lambda_min = 10E-9;% in meters
% lambda_max = 2E-5;% in meters

lambda_camera = linspace(lambda_min,lambda_max,N);
Dlambda_camera = mean(diff(lambda_camera));

rmin = 0;
try
    rmax = ST.params.fields.a;
    Nr = 25;
    Ntheta = 80;
    
    Rmin = 0.9;
    Rmax = 2.1;
    
    Zmin = -0.6;
    Zmax = 0.6;
    
    NR = 50;
    NZ = 50;
catch
    %     Ro = ST.params.fields.Ro;
    %     rmax = max([max(ST.params.fields.R) - Ro, Ro - min(ST.params.fields.R)]);
    rmax = 1.2;
    Nr = 40;
    Ntheta = 80;
    
    Rmin = 0.9;
    Rmax = 2.6;
    
    Zmin = -2.0;
    Zmax = 2.0;
    
    NR = 100;%100
    NZ = 235;%235
end

if strcmp(geometry,'poloidal')
    Psyn_total = zeros(Ntheta,Nr);
else
    Psyn_total = zeros(NZ,NR);
end


lambda_camera = 1E2*lambda_camera;

num_species = double(ST.params.simulation.num_species);

fh = figure;
spectrum_figure = figure;
numPanels = ceil(sqrt(num_species + 1));

L = cell(1,num_species);

% Poloidal distribution of the total radiated power
for ss=1:num_species
% for ss=6:6
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    Ro = ST.params.fields.Ro;
    
    if strcmp(geometry,'poloidal')
        Psyn = zeros(Ntheta,Nr);
    else
        Psyn = zeros(NZ,NR);
    end
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag(:,1:end),2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta(1:end) < 90,2) );
    bool = pin;% & passing;
    
    X = [];
    V = [];
    gammap = [];
    Prad = [];
    eta = [];
    for ii=it1:it2
        X = [X,ST.data.(['sp' num2str(ss)]).X(:,bool,ii)];
        V = [V,ST.data.(['sp' num2str(ss)]).V(:,bool,ii)];
        gammap = [gammap;ST.data.(['sp' num2str(ss)]).gamma(bool,ii)];
        eta = [eta;pi*ST.data.(['sp' num2str(ss)]).eta(bool,ii)/180];
        Prad = [Prad;abs(ST.data.(['sp' num2str(ss)]).Prad(bool,ii))];
    end
    v = squeeze( sqrt( sum(V.^2,1) ) )';
    
    [vp,~,camPos] = identifyVisibleParticles(X,V,gammap,false);
%     vp = true(size(gammap));
%     psi = ones(size(gammap));
    
    X(:,~vp) = [];
    V(:,~vp) = [];
    gammap(~vp) = [];
    v(~vp) = [];
    eta(~vp) = [];
    Prad(~vp) = [];
    
    numPart = size(camPos,1);
    
    if strcmp(geometry,'poloidal')
        % Toroidal coordinates
        % r = radius, theta = poloidal angle, zeta = toroidal angle
        r = squeeze(sqrt( (sqrt(sum(X(1:2,:).^2,1)) - Ro).^2 + X(3,:).^2));
        theta = atan2(squeeze(X(3,:)),squeeze(sqrt(sum(X(1:2,:).^2,1)) - Ro));
        theta(theta<0) = theta(theta<0) + 2*pi;

        Dr = (rmax - rmin)/Nr;
        r_grid = 0.5*Dr + (0:1:(Nr-1))*Dr;

        Dtheta = 2*pi/Ntheta;
        theta_grid = 0.5*Dtheta + (0:1:(Ntheta-1))*Dtheta;

        ir = floor(r/Dr) + 1;
        itheta = floor(theta/Dtheta) + 1;

        % % % Set-up of the grid of the poloidal plane
        x_grid = zeros(Ntheta,Nr);
        y_grid = zeros(Ntheta,Nr);
        for ii=1:Nr
            for jj=1:Ntheta
                x_grid(jj,ii) = Ro + r_grid(ii)*cos(theta_grid(jj));
                y_grid(jj,ii) = r_grid(ii)*sin(theta_grid(jj));
            end
        end
        % Toroidal coordinates
    else
        % Cylindrical coordinates
        R = sqrt(sum(X(1:2,:).^2,1));
        Z = X(3,:);

        DR = (Rmax - Rmin)/NR;
        DZ = (Zmax - Zmin)/NZ;

        R_grid = Rmin + 0.5*DR + (0:1:(NR-1))*DR;
        Z_grid = Zmin + 0.5*DZ + (0:1:(NZ-1))*DZ;

        iR = floor((R - Rmin)/DR) + 1;
        iZ = floor((Z - Zmin)/DZ) + 1;
        % Cylindrical coordinates
    end
    
    % % % Option for calculating the total Psyn WITHOUT wavelength filtering
    if (~filtered)
        for ii=1:numPart
            try
                if strcmp(geometry,'poloidal')
                    Psyn(itheta(ii),ir(ii)) = Psyn(itheta(ii),ir(ii)) + ...
                        Prad(ii);
                    
                    Psyn_total(itheta(ii),ir(ii)) = Psyn_total(itheta(ii),ir(ii)) + ...
                        Prad(ii);
                else
                    Psyn(iZ(ii),iR(ii)) = ...
                        Psyn(iZ(ii),iR(ii)) + Prad(ii);
                    
                    Psyn_total(iZ(ii),iR(ii)) = ...
                        Psyn_total(iZ(ii),iR(ii)) + Prad(ii);
                end
            catch
                disp('An issue calculating poloidal plane')
            end
        end
        
        figure(fh)
        subplot(numPanels,numPanels,ss)
        if strcmp(geometry,'poloidal')
            surf(x_grid,y_grid,Psyn,'LineStyle','none')
            axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid)])
        else
            contourf(R_grid,Z_grid,Psyn,'LineStyle','none')
            axis([min(R_grid) max(R_grid) min(Z_grid) max(Z_grid)])
        end
        colormap(jet(512))
        h = colorbar;
        ylabel(h,'$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
        view([0,90])
        axis equal; box on
        shading interp
        xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
        ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
        title(['$\theta_0 = $' num2str(ST.params.species.etao(ss)) '$^\circ$'],...
            'Interpreter','latex','FontSize',16)
        
    else % WITH filter
        
        try
            B = ST.data.(['sp' num2str(ss)]).B(:,bool,it);
            E = zeros(size(B));
        catch
            B = analyticalB(ST,X);
            E = analyticalE(ST,X);
        end

        vec_mag = zeros(numPart,1);
        Binormal = zeros(numPart,3);
        
        for ii=1:numPart
            VxE = cross(squeeze(V(:,ii)),squeeze(E(:,ii)));
            VxB = cross(squeeze(V(:,ii)),squeeze(B(:,ii)));
            VxVxB = cross(squeeze(V(:,ii)),VxB);
            vec = VxE + VxVxB;
            vec_mag(ii) = sqrt( vec'*vec );
            
            Binormal(ii,:) = vec/vec_mag(ii);
        end
        
        psi = abs(pi/2 - acos(dot(camPos',Binormal')));
        
        % Actual curvature
        k = q*vec_mag./(m*gammap.*v.^3);

        % % % % NO CUTOFF % % %
%         k = ...
%             q*ST.params.fields.Bo*sin(pi*ST.params.species.etao(ss)/180)...
%             /(ST.params.species.gammao(ss)*ST.params.species.m(ss)*v(1))*ones(size(v));
        
        % % % % Beyond this point all variables are in cgs units % % % %
        c = 1E2*ST.params.scales.v;
        qe = 3E9*q;
        m = 1E3*m;
        
        k = k/1E2;
        
        lambdac = (4/3)*pi*(1./gammap).^3./k;
        
%         lambda_camera = linspace(1E-7,5*lambdac(1),N); % % % % NO CUTOFF % % %
        
        I = find(lambdac > lambda_min);
        numEmittingPart = numel(I);
         
        Psyn_camera = zeros(numEmittingPart,N);
        disp('Decomposing radiation in wavelengths...')
        
        C0 = 4*pi*c*qe^2/sqrt(3);
        fun = @(x) besselk(5/3,x);
        y = (gammap.*psi').^2;
        for ii=1:numEmittingPart
            ind = I(ii);
            for jj=1:N
                lower_integration_limit = lambdac(ind)/lambda_camera(jj);
%                 if (lambda_camera(jj) < lambdac(ind)) && (lower_integration_limit < upper_integration_limit)             
%                 if (lower_integration_limit < upper_integration_limit)
%                     Q = integral(fun,lower_integration_limit,upper_integration_limit);
%                     C1 = 1/(gammap(ind)^2*lambda_camera(jj)^3);
%                     Psyn_camera(ii,jj) =  C0*C1*Q;
%                 end

                if (lambda_camera(jj) < lambdac(ind)) && isfinite(lambdac(ind))
                    zeta = 0.5*lower_integration_limit*(1 + y(ind))^(3/2);
                    D0 = 3*c*qe^2*k(ind)/(2*pi*lambda_camera(jj)^2);
                    
                    Psyn_camera(ii,jj) = ...
                        D0*lower_integration_limit^2*gammap(ind)^2*(1 + y(ind))^2*(besselk(2/3,zeta).^2 + ...
                        (y(ind)/(1 + y(ind)))*besselk(1/3,zeta).^2);
                end
            end
        end
        Psyn_camera = Psyn_camera/numEmittingPart;
                
        % With the data of Psyn_camera we make the poloidal plane
        % 'photograph'.
        disp('Calculating poloidal plane...')
        for ii=1:numEmittingPart
            ind = I(ii);
            try
                Psyn_integrated = trapz(lambda_camera,Psyn_camera(ii,:));
                if strcmp(geometry,'poloidal')
                    Psyn(itheta(ind),ir(ind)) = ...
                        Psyn(itheta(ind),ir(ind)) + Psyn_integrated;
                    
                    Psyn_total(itheta(ind),ir(ind)) = ...
                        Psyn_total(itheta(ind),ir(ind)) + Psyn_integrated;
                else
                    Psyn(iZ(ind),iR(ind)) = ...
                        Psyn(iZ(ind),iR(ind)) + Psyn_integrated;
                    
                    Psyn_total(iZ(ind),iR(ind)) = ...
                        Psyn_total(iZ(ind),iR(ind)) + Psyn_integrated;
                end
            catch
                disp('An issue calculating poloidal plane')
            end
        end
        
        lch = 1E7;
        Pch = 1E-7;
        
        lambdac = lch*lambdac;
        Psyn_camera = Pch*Psyn_camera;
        Psyn = Pch*Psyn;
        k = 1E2*k;
        
        % % % % Beyond this point all variables are in SI units % % % %      

        figure(fh)
        subplot(numPanels,numPanels,ss)
        if strcmp(geometry,'poloidal')
            surf(x_grid,y_grid,Psyn,'LineStyle','none')
%             axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid)])
        else
%             surfc(R_grid,Z_grid,Psyn,'LineStyle','none')
            contourf(R_grid,Z_grid,Psyn,'LineStyle','none')
%             axis([min(R_grid) max(R_grid) min(Z_grid) max(Z_grid)])
        end
        colormap(jet(512))
        h = colorbar;
        ylabel(h,'$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
        view([0,90])
        axis equal; box on
        shading interp
        xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
        ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
        title(['$\theta_0 = $' num2str(ST.params.species.etao(ss)) '$^\circ$'],...
            'Interpreter','latex','FontSize',16)
        
        figure(spectrum_figure)
        hold on
        plot(lch*lambda_camera,sum(Psyn_camera,1)/lch)
        hold off
        grid on;box on
        xlabel('$\lambda$ (nm)','Interpreter','latex','FontSize',16)
        ylabel('$P_{syn}$ (Watts/m)','Interpreter','latex','FontSize',16)
        L{ss} = ['$\theta_0 = $' num2str(ST.params.species.etao(ss)) '$^\circ$'];
    end
    
end



% % % % Final figures % % % %
if (~filtered)
    figure(fh)
    subplot(numPanels,numPanels,num_species+1)
    if strcmp(geometry,'poloidal')
        surf(x_grid,y_grid,Psyn_total_pol,'LineStyle','none')
%         axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid)])
    else
        surf(R_grid,Z_grid,Psyn_total,'LineStyle','none')
%         axis([min(R_grid) max(R_grid) min(Z_grid) max(Z_grid)])
    end
    colormap(jet(512))
    h = colorbar;
    ylabel(h,'$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
    view([0,90])
    axis equal; box on
    shading interp
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    title('Total $P_{syn}$','Interpreter','latex','FontSize',16)
else
    Psyn_total = 1E-7*Psyn_total;
    
    figure(fh)
    subplot(numPanels,numPanels,num_species+1)
    if strcmp(geometry,'poloidal')
        surf(x_grid,y_grid,Psyn_total_pol,'LineStyle','none')
        axis([min(x_grid) max(x_grid) min(y_grid) max(y_grid)])
    else
        surf(R_grid,Z_grid,Psyn_total,'LineStyle','none')
        axis([min(R_grid) max(R_grid) min(Z_grid) max(Z_grid)])
    end
    colormap(jet(512))
    h = colorbar;
    ylabel(h,'$P_{syn}$ (Watts)','Interpreter','latex','FontSize',16)
    view([0,90])
    axis equal; box on
    shading interp
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    title('Total $P_{syn}$','Interpreter','latex','FontSize',16)
    
    figure(spectrum_figure)
    legend(L,'Interpreter','latex','FontSize',12)
end

disp('Spectrum of synchrotron radiation: done!')
end

function VS = identifyVisibleParticles_prototype(ST)
VS = struct;

% Radial position of inner wall
Riw = 1; % in meters

% Radial and vertical position of the camera
Rc = 2.38; % in meters
Zc = 0.076; % in meters

it = ST.params.simulation.num_snapshots + 1;

num_species = double(ST.params.simulation.num_species);

h = figure;
numPanels = ceil(sqrt(num_species));

for ss=1:num_species
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag,2)); % confined particles
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta < 90,2) ); % passing particles
    % If bool = pin & passing, we consider confined passing particles
    % If bool = pin, we consider passing and trapped particles
    bool = pin;
    np = numel(find(bool == 1));
    
    X = squeeze(ST.data.(['sp' num2str(ss)]).X(:,bool,it));
    V = squeeze(ST.data.(['sp' num2str(ss)]).V(:,bool,it));
    gammap = ST.data.(['sp' num2str(ss)]).gamma(bool,it);
    
    % mea = maximum emission angle of synchrotron radiation
    mea = 1./gammap; % in rad. \psi ~ 1/\gamma
    
    xo = X(1,:);
    yo = X(2,:);
    zo = X(3,:);
    
    v = sqrt(sum(V.^2,1));
    vox = V(1,:)./v;
    voy = V(2,:)./v;
    voz = V(3,:)./v;
    
    % First we find the Z position where V hits the wall at Rc
    theta_f = zeros(1,np);
    Z_f = zeros(1,np);
    hitInnerWall = true(1,np);
    for ii=1:np
        p = zeros(1,3);
        
        % polinomial coefficients p(1)*x^2 + p(2)*x + p(3) = 0
        p(1) = vox(ii)^2 + voy(ii)^2;
        p(2) = 2*(xo(ii)*vox(ii) + yo(ii)*voy(ii));
        p(3) = xo(ii)^2 + yo(ii)^2 - Rc^2;
        
        r = roots(p);
        if all(r>0) || all(r<0)
            error(['Something wrong at ii=' num2str(ii)])
        end        
        t = max(r);
        Z_f(ii) = zo(ii) + voz(ii)*t;
        theta_f(ii) = atan2(yo(ii) + voy(ii)*t,xo(ii) + vox(ii)*t);
        if (theta_f(ii) < 0)
            theta_f(ii) = theta_f(ii) + 2*pi;
        end
        
        p(3) = xo(ii)^2 + yo(ii)^2 - Riw^2;
        r = roots(p);
        if isreal(r) && any(r>0)
            hitInnerWall(ii) = false;
        end
    end
    
    Ro = sqrt(xo(~hitInnerWall).^2 + yo(~hitInnerWall).^2);
    Zo = zo(~hitInnerWall);
    figure(h);
    subplot(numPanels,numPanels,ss)
    plot(Ro,Zo,'r.')
    axis equal
    
    xo(~hitInnerWall) = [];
    yo(~hitInnerWall) = [];
    zo(~hitInnerWall) = [];
    vox(~hitInnerWall) = [];
    voy(~hitInnerWall) = [];
    voz(~hitInnerWall) = [];
    theta_f(~hitInnerWall) = [];
    mea(~hitInnerWall) = [];
    
    nvp = numel(find(hitInnerWall == false));
    
    Ro = sqrt(xo.^2 + yo.^2);
    figure(h);
    subplot(numPanels,numPanels,ss)
    hold on
    plot(Ro,zo,'k.')
    hold off
    
    % Then, we calculate the angle between V and the position of the camera
    xc = Rc*cos(theta_f);
    yc = Rc*sin(theta_f);
    
    ax = xc - xo;
    ay = yc - yo;
    az = Zc - zo;
    
    a = sqrt(ax.^2 + ay.^2 + az.^2);
    
    ax = ax./a;
    ay = ay./a;
    az = az./a;
    
    angle = acos(ax.*vox + ay.*voy + az.*voz)';
    visible = angle < mea;
    
    Ro = sqrt(xo(visible).^2 + yo(visible).^2);
    figure(h);
    subplot(numPanels,numPanels,ss)
    hold on
    plot(Ro,zo(visible),'g.','MarkerSize',18)
    hold off
    
end


end

function [vp,psi,camPos] = identifyVisibleParticles(X,V,gammap,option)

% Radial position of inner wall
Riw = 1; % in meters

% Radial and vertical position of the camera
Rc = 2.38; % in meters
Zc = 0.076; % in meters
% Zc = 0;

np = numel(gammap);

% mea = maximum emission angle of synchrotron radiation
mea = 1./gammap; % in rad. \psi ~ 1/\gamma

xo = X(1,:);
yo = X(2,:);
zo = X(3,:);

Ro = sqrt(sum(X(1:2,:).^2,1));

v = sqrt(sum(V.^2,1));
vox = V(1,:)./v;
voy = V(2,:)./v;
voz = V(3,:)./v;

tmp_cam = zeros(np,3);

% First we find the Z position where V hits the wall at Rc
theta_f = zeros(1,np);
% Z_f = zeros(1,np);
hitInnerWall = false(1,np);
for ii=1:np
    if (Ro(ii) < Rc)
        p = zeros(1,3);
        
        % polinomial coefficients p(1)*x^2 + p(2)*x + p(3) = 0
        p(1) = vox(ii)^2 + voy(ii)^2;
        p(2) = 2*(xo(ii)*vox(ii) + yo(ii)*voy(ii));
        p(3) = xo(ii)^2 + yo(ii)^2 - Rc^2;
        
        r = roots(p);
        if all(r>0) || all(r<0)
            disp(['Something wrong at ii=' num2str(ii)])
            hitInnerWall(ii) = true;
        else
            t = max(r);
            
            % Cartesian coordinates where the tangent vector intersects the
            % outter wall of the axisymmetric device
            X_f = xo(ii) + vox(ii)*t;
            Y_f = yo(ii) + voy(ii)*t;
%             Z_f = zo(ii) + voz(ii)*t;
            
            theta_f(ii) = atan2(Y_f,X_f);
            if (theta_f(ii) < 0)
                theta_f(ii) = theta_f(ii) + 2*pi;
            end
            
            X_cam = Rc*cos(theta_f(ii));
            Y_cam = Rc*sin(theta_f(ii));
            
            tmp_cam(ii,:) = [X_cam - xo(ii), Y_cam - yo(ii), Zc - zo(ii)];
            tmp_cam(ii,:) = tmp_cam(ii,:)/sqrt(tmp_cam(ii,:)*tmp_cam(ii,:)');
            
            % % % Check if the unitary velocity vector hits the inner wall
            p(3) = xo(ii)^2 + yo(ii)^2 - Riw^2; 
            r = roots(p);
            if isreal(r) && any(r>0)
                hitInnerWall(ii) = true;
            end
        end
    else
        % The particle hits the outer wall
        hitInnerWall(ii) = true;
    end
end

I = find(hitInnerWall == false);

% Then, we calculate the angle between V and the position of the camera
xc = Rc*cos(theta_f(I));
yc = Rc*sin(theta_f(I));

ax = xc - xo(I);
ay = yc - yo(I);
az = Zc - zo(I);

a = sqrt(ax.^2 + ay.^2 + az.^2);

ax = ax./a;
ay = ay./a;
az = az./a;

psi = acos(ax.*vox(I) + ay.*voy(I) + az.*voz(I))';
visible = psi <= mea(I);

if ( option )
    Ro = sqrt(xo(I(visible)).^2 + yo(I(visible)).^2);
    figure
    plot(Ro,zo(I(visible)),'g.','MarkerSize',18)
end

vp = false(1,np);
vp(I(visible)) = true;

psi(~visible) = [];

camPos = tmp_cam(I(visible),:);

end

function [ip,X] = findVisibleParticles(X,V,angle,camera_params,option)

% Radial position of inner wall (central post)
Riw = 1; % in meters

% Radial and vertical position of the camera
Rc = camera_params.position(1); % in meters
Zc = camera_params.position(2); % in meters

np = numel(angle);

xo = X(1,:);
yo = X(2,:);
zo = X(3,:);

Ro = sqrt(sum(X(1:2,:).^2,1));
y_camera = -camera_params.focal_length*(zo - camera_params.position(2))./(camera_params.position(1) - Ro);

v = sqrt(sum(V.^2,1));
vox = V(1,:)./v;
voy = V(2,:)./v;
voz = V(3,:)./v;

% First we find the Z position where V hits the wall at Rc
tmp_cam = zeros(np,3);
theta_f = zeros(1,np);
discarded = false(1,np);
X_f = zeros(1,np);
Y_f = zeros(1,np);
Z_f = zeros(1,np);
for ii=1:np
    p = zeros(1,3);
    
    % polinomial coefficients p(1)*x^2 + p(2)*x + p(3) = 0
    p(1) = vox(ii)^2 + voy(ii)^2;
    p(2) = 2*(xo(ii)*vox(ii) + yo(ii)*voy(ii));
    p(3) = xo(ii)^2 + yo(ii)^2 - Rc^2;
    
    r = roots(p);
    if ~isreal(r) || (all(r>0) || all(r<0))
        disp(['Something wrong at ii=' num2str(ii)])
        discarded(ii) = true;
    else
        t = max(r);
        
        if (y_camera(ii) < camera_params.pixel_grid.ymin) || (y_camera(ii) > camera_params.pixel_grid.ymax)
            discarded(ii) = true;
        else
            % Cartesian coordinates where the tangent vector intersects the
            % outter wall of the axisymmetric device
            X_f(ii) = xo(ii) + vox(ii)*t;
            Y_f(ii) = yo(ii) + voy(ii)*t;
            Z_f(ii) = zo(ii) + voz(ii)*t;
            
            theta_f(ii) = atan2(Y_f(ii),X_f(ii));
            if (theta_f(ii) < 0)
                theta_f(ii) = theta_f(ii) + 2*pi;
            end
            
            X_cam = Rc*cos(theta_f(ii));
            Y_cam = Rc*sin(theta_f(ii));
            
            tmp_cam(ii,:) = [X_cam - xo(ii), Y_cam - yo(ii), Zc - zo(ii)];
            tmp_cam(ii,:) = tmp_cam(ii,:)/sqrt(tmp_cam(ii,:)*tmp_cam(ii,:)');
            
            % % % Check if the unitary velocity vector hits the inner wall
            p(3) = xo(ii)^2 + yo(ii)^2 - Riw^2;
            s = roots(p);
            if isreal(s) && any(s>0)
                discarded(ii) = true;
            end
        end
    end
end

I = find(discarded == false);

% Then, we calculate the angle between V and the position of the camera
xc = Rc*cos(theta_f(I));
yc = Rc*sin(theta_f(I));

ax = xc - xo(I);
ay = yc - yo(I);
az = Zc - zo(I);

a = sqrt(ax.^2 + ay.^2 + az.^2);

ax = ax./a;
ay = ay./a;
az = az./a;

psi = acos(ax.*vox(I) + ay.*voy(I) + az.*voz(I))';
visible = psi <= angle(I);

II = I(visible);

theta_f = theta_f(I(visible));
XX = [X_f;Y_f;Z_f];
for jj=1:numel(II)
    rotation_matrix = [cos(theta_f(jj)),sin(theta_f(jj));
    -sin(theta_f(jj)),cos(theta_f(jj))];
    
    X(1:2,II(jj)) = rotation_matrix*X(1:2,II(jj));
    
    XX(1:2,II(jj)) = rotation_matrix*XX(1:2,II(jj));
end

if ( option )
    Ro = sqrt(xo(II).^2 + yo(II).^2);
    figure
    subplot(2,2,1)
    plot(Ro,zo(II),'g.','MarkerSize',5)
    xlabel('$R$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    subplot(2,2,2)
    histogram(zo(II))
    xlabel('$Z$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$f(Z)$','Interpreter','latex','FontSize',16)
    subplot(2,2,3)
    plot3([X(1,II);XX(1,II)],[X(2,II);XX(2,II)],[X(3,II);XX(3,II)])
    xlabel('$X$ (m)','Interpreter','latex','FontSize',16)
    ylabel('$Y$ (m)','Interpreter','latex','FontSize',16)
    axis equal; box on
    subplot(2,2,4)
    plot(theta_f,'b','MarkerSize',5)
    xlabel('Particle number','Interpreter','latex','FontSize',16)
    ylabel('$\theta$ (rad)','Interpreter','latex','FontSize',16)
end

ip = false(1,np);
ip(II) = true;
end

function pixel_grid = setupCameraPixelGrid(camera_params,option)
% Here we set-up the pixel grid of the camera in a coordinate system such
% that the x-axis corresponds to the horizontal dimension, and the y-axis
% corresponds to the vertical dimension, with respect to the floor level.
xmin = -camera_params.size(1)/2;
xmax = camera_params.size(1)/2;
DX = camera_params.size(1)/camera_params.NX;
pixels_x = (xmin+0.5*DX):DX:(xmax-0.5*DX);

ymin = camera_params.position(2) - camera_params.size(2)/2;
ymax = camera_params.position(2) + camera_params.size(2)/2;
DY = camera_params.size(2)/camera_params.NY;
pixels_y = (ymin+0.5*DY):DY:(ymax-0.5*DY);

pixel_grid = struct;
[pixel_grid.X, pixel_grid.Y] = meshgrid(pixels_x,pixels_y);

pixel_grid.xmin = xmin;
pixel_grid.xmax = xmax;

pixel_grid.ymin = ymin;
pixel_grid.ymax = ymax;

pixel_grid.xnodes = pixels_x;
pixel_grid.ynodes = pixels_y;

pixel_grid.xedges = linspace(xmin,xmax,camera_params.NX+1);
pixel_grid.yedges = linspace(ymin,ymax,camera_params.NY+1);

if option
    [X,Y] = ...
        meshgrid(pixel_grid.xedges,pixel_grid.yedges);
    
    h = figure;
    subplot(1,2,1)
    plot(pixel_grid.X,pixel_grid.Y,'rs','MarkerSize',2)
    hold on
    plot(X,Y,'k',X',Y','k')
    hold off
    axis equal
    title('Camera pixel array','Interpreter','latex','FontSize',12)
    xlabel('$X$ (m)','Interpreter','latex','FontSize',14)
    ylabel('$Y$ (m)','Interpreter','latex','FontSize',14)
    
    
    % xy-plane with geometry of camera
    Rmax = ceil(camera_params.position(1));
    R = 0.5:0.5:Rmax;
    t = linspace(0,2*pi,50);
    for ii=1:numel(R)
        x = R(ii)*cos(t);
        y = R(ii)*sin(t);
        
        figure(h)
        subplot(1,2,2)
        hold on
        plot(x,y,'--','Color',[0.4,0.4,0.4])
        hold off
        axis equal
        box on
    end
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(camera_params.position(1),0,'ko','MarkerSize',5)
    hold off
    
    N = 10;
    ccd_x = linspace(xmin,xmax,N);
    ccd_y = -camera_params.focal_length*ones(1,N);
    
    a = camera_params.incline;
    
    ccd_x_r = ccd_x*cos(a) - ccd_y*sin(a);
    ccd_y_r = ccd_x*sin(a) + ccd_y*cos(a);
    
    ccd_x_r = ccd_x_r + camera_params.position(1);
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(ccd_x_r,ccd_y_r,'k','LineWidth',2)
    hold off
    
    % Main line of sight
    line_x = zeros(1,N);
    line_y = linspace(-camera_params.focal_length,Rmax,N);
    
    line_x_r = line_x*cos(a) - line_y*sin(a);
    line_y_r = line_x*sin(a) + line_y*cos(a);
    
    line_x_r = line_x_r + camera_params.position(1);
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(line_x_r,line_y_r,'r-.','LineWidth',1)
    hold off
    
    % Cone of sight min
    line_x = zeros(1,N);
    line_y = linspace(0,Rmax+camera_params.focal_length,N);
    
    b = camera_params.horizontal_angle_view;
    
    line_x_r = line_x*cos(b) + line_y*sin(b);
    line_y_r = -line_x*sin(b) + line_y*cos(b);
        
    line_x_r = line_x_r + xmin;
    line_y_r = line_y_r - camera_params.focal_length;
    
    line_x = line_x_r;
    line_y = line_y_r;
    line_x_r = line_x*cos(a) - line_y*sin(a);
    line_y_r = line_x*sin(a) + line_y*cos(a);
    
    line_x_r = line_x_r + camera_params.position(1);
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(line_x_r,line_y_r,'r','LineWidth',1)
    hold off
    
    % Cone of sight max
    line_x = zeros(1,N);
    line_y = linspace(0,Rmax+camera_params.focal_length,N);
    
    line_x_r = line_x*cos(b) - line_y*sin(b);
    line_y_r = line_x*sin(b) + line_y*cos(b);
        
    line_x_r = line_x_r + xmax;
    line_y_r = line_y_r - camera_params.focal_length;
    
    line_x = line_x_r;
    line_y = line_y_r;
    line_x_r = line_x*cos(a) - line_y*sin(a);
    line_y_r = line_x*sin(a) + line_y*cos(a);
    
    line_x_r = line_x_r + camera_params.position(1);
    
    figure(h)
    subplot(1,2,2)
    hold on
    plot(line_x_r,line_y_r,'r','LineWidth',1)
    hold off
end


end

function [rotation_angle,ip_in_pixel] = findRotationAngles(X,camera_params)
%%
% Here we compute the rotation angles in the toroidal direction based on
% the pixel array of the camera.

np = size(X,2);
Rcam = camera_params.position(1);
incline = camera_params.incline;

xnodes = camera_params.pixel_grid.xnodes;
xedges = camera_params.pixel_grid.xedges;

ynodes = camera_params.pixel_grid.ynodes;
yedges = camera_params.pixel_grid.yedges;

eta_range = zeros(2,camera_params.NX);
beta_range = zeros(2,camera_params.NX);
for ii=1:camera_params.NX
    eta_range(:,ii) = atan2(abs([xedges(ii);xedges(ii+1)]),camera_params.focal_length);
    if (xnodes(ii) < 0)
        beta_range(:,ii) = pi/2 - incline - eta_range(:,ii);
    else
        beta_range(:,ii) = pi/2 - eta_range(:,ii);
    end
end
cos_beta_range = cos(beta_range);

psi_range = zeros(2,camera_params.NY);
for ii=1:camera_params.NY
    psi_range(:,ii) = atan2([yedges(ii);yedges(ii+1)],camera_params.focal_length);
end

rotation_matrix = [cos(incline),sin(incline);
    -sin(incline),cos(incline)];

P = [X(1,:) - Rcam;X(2,:)];
tmp = sqrt(sum(P.^2,1));
P(1,:) = P(1,:)./tmp;
P(2,:) = P(2,:)./tmp;

PP = zeros(2,np);
for jj=1:np
    PP(:,jj) = rotation_matrix*P(:,jj);
end

ipos = PP(1,:) > 0;
ineg = PP(1,:) < 0;

cos_beta = zeros(1,np);
cos_beta(ineg) = -P(1,ineg);
cos_beta(ipos) = PP(1,ipos);

D = sqrt( (X(1,:) - Rcam).^2 + X(2,:).^2 );
psi = atan2(X(3,:),D);

ip_in_pixel_x = cell(1,camera_params.NX);
rotation_angle = cell(1,camera_params.NX);
for ii=1:camera_params.NX
    ip_in_pixel_y = cell(1,camera_params.NY);
    for jj=1:camera_params.NY
        I1 = find(psi >= psi_range(1,jj));
        I2 = find(psi < psi_range(2,jj));
        ip_in_pixel_y{jj} = intersect(I1,I2);
    end
    
    if (xnodes(ii) < 0)
        I1 = find(cos_beta <= cos_beta_range(1,ii));
        I2 = find(cos_beta > cos_beta_range(2,ii));
        ip_in_pixel_x{ii} = intersect(I1,I2);
        
        if ~isempty(ip_in_pixel_x{ii})
            P = X(1:2,ip_in_pixel_x{ii});
            
            Ro = -P;
            tmp = sqrt(sum(Ro.^2,1));
            Ro(1,:) = Ro(1,:)./tmp;
            Ro(2,:) = Ro(2,:)./tmp;
            
            Rc = [Rcam - P(1,:);-P(2,:)];
            tmp = sqrt(sum(Rc.^2,1));
            Rc(1,:) = Rc(1,:)./tmp;
            Rc(2,:) = Rc(2,:)./tmp;
            
            cos_theta = dot(Ro,Rc,1);
            rotation_angle{ii} = pi - acos(cos_beta(ip_in_pixel_x{ii})) - acos(cos_theta);
        end
    else
        I1 = find(cos_beta >= cos_beta_range(1,ii));
        I2 = find(cos_beta < cos_beta_range(2,ii));
        ip_in_pixel_x{ii} = intersect(I1,I2);
        
        if ~isempty(ip_in_pixel_x{ii})
            P = X(1:2,ip_in_pixel_x{ii});
            
            Ro = -P;
            tmp = sqrt(sum(Ro.^2,1));
            Ro(1,:) = Ro(1,:)./tmp;
            Ro(2,:) = Ro(2,:)./tmp;
            
            Rc = [Rcam - P(1,:);-P(2,:)];
            tmp = sqrt(sum(Rc.^2,1));
            Rc(1,:) = Rc(1,:)./tmp;
            Rc(2,:) = Rc(2,:)./tmp;
            
            cos_theta = dot(Ro,Rc,1);
            rotation_angle{ii} = acos(cos_beta(ip_in_pixel_x{ii})) + incline - acos(cos_theta);
        end
    end
end

ip_in_pixel = cell(camera_params.NX,camera_params.NY);
for ii=1:camera_params.NX
    for jj=1:camera_params.NY
        ip_in_pixel{ii,jj} = intersect(ip_in_pixel_x{ii},ip_in_pixel_y{jj});
    end
end

end

function SD = syntheticDiagnosticSynchrotron(ST)
% Synthetic diagnostic of camera for synchrotron radiation.
SD = struct;

% Here we define many snapshots will be used
numSnapshots = 10;
it1 = ST.params.simulation.num_snapshots + 1 - numSnapshots;
it2 = ST.params.simulation.num_snapshots + 1;

% Simulation parameters
num_species = double(ST.params.simulation.num_species);

% Resolution in wavelength
N = 200;
% Wavelength range (450 nm - 950 nm for visible light).
lambda_min = 450E-9;% in meters
lambda_max = 950E-9;% in meters
lambda = linspace(lambda_min,lambda_max,N);
lambda = 1E2*lambda; % in cm

% Camera parameters
camera_params = struct;
camera_params.NX = 30;
camera_params.NY = 25;
camera_params.size = [0.25,0.2]; % [horizontal size, vertical size] in meters
camera_params.focal_length = 0.4; % In meters
camera_params.position = [2.38,0.0]; % [R,Z] in meters
% The angle defined by the detector plane (pixel array) and the x-axis of a
% coordinate system where phi = 0, the toroidal angle, corresponds to the
% y-axis, and phi = 90 corresponds to the x-axis
camera_params.incline = 55; % in degrees
camera_params.incline = deg2rad(camera_params.incline);
camera_params.horizontal_angle_view = ...
    atan2(0.5*camera_params.size(1),camera_params.focal_length); % in radians
camera_params.vertical_angle_view = ...
    atan2(0.5*camera_params.size(2),camera_params.focal_length); % in radians
camera_params.pixel_grid = setupCameraPixelGrid(camera_params,true);

for ss=2:num_species
    q = abs(ST.params.species.q(ss));
    m = ST.params.species.m(ss);
    Ro = ST.params.fields.Ro;
    
    Psyn = zeros(camera_params.NX,camera_params.NY);
    
    pin = logical(all(ST.data.(['sp' num2str(ss)]).flag(:,it1:it2),2));
    passing = logical( all(ST.data.(['sp' num2str(ss)]).eta(it1:it2) < 90,2) );
    bool = pin;% & passing;
    
    X = [];
    V = [];
    gammap = [];
    Prad = [];
    eta = [];
    B = [];
    E = [];
    for ii=it1:it2
        X = [X,ST.data.(['sp' num2str(ss)]).X(:,bool,ii)];
        V = [V,ST.data.(['sp' num2str(ss)]).V(:,bool,ii)];
        gammap = [gammap;ST.data.(['sp' num2str(ss)]).gamma(bool,ii)];
        eta = [eta;pi*ST.data.(['sp' num2str(ss)]).eta(bool,ii)/180];
        Prad = [Prad;abs(ST.data.(['sp' num2str(ss)]).Prad(bool,ii))];
        B = [B,ST.data.(['sp' num2str(ss)]).B(:,bool,ii)];
        try
            E = [E,ST.data.(['sp' num2str(ss)]).B(:,bool,ii)];
        catch
            E = zeros(size(B));
        end
    end
    v = squeeze( sqrt( sum(V.^2,1) ) )';
    
    numPart = numel(eta);
    
    vec = zeros(numPart,1);
    for ii=1:numPart
        VxE = cross(squeeze(V(:,ii)),squeeze(E(:,ii)));
        VxB = cross(squeeze(V(:,ii)),squeeze(B(:,ii)));
        VxVxB = cross(squeeze(V(:,ii)),VxB);
        aux = VxE + VxVxB;
        vec(ii) = sqrt( aux'*aux );
    end
    clear VxE VxB VxVxB aux E B
    
    k = q*vec./(m*gammap.*v.^3); % Actual curvature
    
    % % % % Beyond this point all variables are in cgs units % % % %
    c = 1E2*ST.params.scales.v;
    qe = 3E9*q;
    m = 1E3*m;
    
    k = k/1E2;
    
    % This angle is used as the primer criterium for deciding whether the
    % particle is seen by the camera or not.
    angle = (3*k*lambda(end)/(4*pi)).^(1/3);
    
    [ip,X] = findVisibleParticles(X,V,angle,camera_params,false);
    
    % Here we drop the particles that are not seen by the camera
    X(:,~ip) = [];
    V(:,~ip) = [];
    gammap(~ip) = [];
    v(~ip) = [];
    eta(~ip) = [];
    Prad(~ip) = [];
    
    [rotation_angles,ip_in_pixel] = findRotationAngles(X,camera_params);
    

    
    numVisiblePart = numel(eta);    
end


end