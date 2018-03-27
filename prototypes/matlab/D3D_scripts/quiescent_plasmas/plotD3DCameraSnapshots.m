close all
clear variables

load('165826_camdata.mat')
% load('165139_camdata.mat')

pixelPlot=1;

%% Pre-amble

for k=1:length(flags.xLim)
    [~, movieLims(k)]=min(abs(fastcam.time/1e3-flags.xLim(k)));
end
movieBins=movieLims(1):1:movieLims(2);

time = fastcam.time(movieBins)/1E3-in.tShift/1E3;

% Shot 165826
% time1 = 2.0;
% time2 = 2.25;
% time3 = 2.5;

% Shot 165826
% time1 = 0.15;
% time2 = 2.0;
% time3 = 2.3;

time1 = 1.25;
time2 = 1.5;
time3 = 1.7;

it1 = find(time>time1,1);
it2 = find(time>time2,1);
it3 = find(time>time3,1);

flags.plotBin = [it1 it2 it3]; 
%% Make Simpler Profile Figure (WITH EFIT AND RAYMAP)


lin.style2={'-k','-r','-b','-g','-m','-y','-c'};
lin.style3={'--k','--r','--b','--g','--m','--y','--c'};
lin.colors={'k','r','b','g','m','y','c'};

% hf = figure('Position',[300 100 1300 600])
hf = figure;
hs = figure;

if flags.plotNegCounts
    colormap(darkb2r(in.zContourLim(1),in.zContourLim(2)));
end

Riw = 1.016;

% First Row
for k=1:length(in.ts)
    if pixelPlot==0
        
        I0 = find(raymap.majz(1,:)>0,1);
        R0 = raymap.majr(:,I0);
        Z0 = raymap.majz(:,I0);
        
        I1 = find(raymap.majz(1,:)>0.25,1);
        R1 = raymap.majr(:,I1);
        Z1 = raymap.majz(:,I1);
        
        I2 = find(raymap.majz(1,:)>-0.25,1);
        R2 = raymap.majr(:,I2);
        Z2 = raymap.majz(:,I2);
        
        I3 = find(raymap.majz(1,:)>0.4,1);
        R3 = raymap.majr(:,I3);
        Z3 = raymap.majz(:,I3);
        
        I4 = find(raymap.majz(1,:)>-0.4,1);
        R4 = raymap.majr(:,I4);
        Z4 = raymap.majz(:,I4);
        
        DZ = diff(raymap.majz,1,1);
        max_DZ = max(DZ,[],1);
        [~,I] = min(abs((max_DZ)));
        RI = raymap.majr(:,I);
        ZI = raymap.majz(:,I);
        
        [~,J] = min(abs(raymap.majr(:,I)-G.rmaxis));
        RJ = raymap.majr(J,:);
        ZJ = raymap.majz(J,:);
        
        [~,K] = min(abs(raymap.majr(:,I) - Riw));
        RK = raymap.majr(K,:);
        ZK = raymap.majz(K,:);
        
        A = squeeze(fastcam.(fastcam.moviToUse)(:,:,flags.plotBin(k)));
        B = reshape(A,[numel(A),1]);
        B(B<=0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,50);
        
        
        figure(hf)
        subplot(4,length(in.ts),k+length(in.ts)*[0 1 2])  
        title(strcat(num2str(in.shot),'.0',num2str(in.ts(k))),...
            'FontSize',flags.fontSize,'Color',lin.colors{k})
        hold on
        
        contourf(raymap.majr, raymap.majz,A,v,'LineStyle','none')
        ax = gca;caxis(ax.CLim);colormap(jet(1024));
        plot(R0,Z0,'-r',R1,Z1,'-g',R2,Z2,'-b',R3,Z3,'-k',R4,Z4,'-m','LineWidth',2)
        
        % Plot the box
%         for kk=1:length(in.box)
%             plot(raymap.majr(in.box(kk).xSlice+in.box(kk).size*[1 1 -1 -1 1],1),...
%                 raymap.majz(1,in.box(kk).ySlice+in.box(kk).size*[-1 1 1 -1 -1]),...
%                 lin.style2{kk},'LineWidth',2)
%         end
        axis equal;xlim([0.95 2.45]);ylim(in.mapYlim)
        
        
        xlabel('R (m)','FontSize',flags.fontSize)
        if k==1;  ylabel('Z (m)','FontSize',flags.fontSize); end
        set(get(gcf,'CurrentAxes'),'FontSize',flags.fontSize)
        grid on; box on
        
        C = max(A(:,I0));
        A0 = A(:,I0)/C;
        A1 = A(:,I1)/C;
        A2 = A(:,I2)/C;
        A3 = A(:,I3)/C;
        A4 = A(:,I4)/C;
        
        figure(hs);
        subplot(1,length(in.ts),k)
        plot(R0,A0,'-.r',R1,A1,'-.g',R2,A2,'-.b',R3,A3,'-.k',R4,A4,'-.m','LineWidth',1)
%         legend({['$Z=$' num2str(zo)],['$Z=$' num2str(Z1)],['$Z=$' num2str(Z2)],...
%             ['$Z=$' num2str(z3)],['$Z=$' num2str(z4)]},'Interpreter','latex')
        title(strcat(num2str(in.shot),'.0',num2str(in.ts(k))),...
            'FontSize',flags.fontSize,'Color',lin.colors{k})
        xlabel('$R$ (m)','FontSize',flags.fontSize,'Interpreter','latex')
        xlim([1 2.5])
        
                % if isfield(in,'EFIT') % EFIT OVERLAY
        
        % Define G
        if flags.evolveEFIT; G=Gfull.gdata(flags.GplotBin(k)); end
        
        figure(hf)
        plot(G.rmaxis,G.zmaxis,'sk','LineWidth',5)% Magnetic Axis
        % Overlay rational surfaces
        qToPlot=[1 1.5 2 3 4];
        G.psiLevels=linspace(G.psimag,G.psibry,length(G.qpsi));
        for jj=1:length(qToPlot)
            [~, qBins(jj)]=min(abs(G.qpsi-qToPlot(jj)));
        end
        %flags.cols2=colormap('copper');
        [G.Rgrid, G.Zgrid]=meshgrid(G.rg, G.zg);
        contour(G.Rgrid,G.Zgrid,G.psizr,...
            G.psiLevels(qBins),'--',...
            'LineWidth',1,'Color','k')
        
        figure(hf)
        plot(G.xlim, G.ylim,'-k','LineWidth',2) % LIMITER
        plot(G.rbbbs(G.rbbbs>0),G.zbbbs(G.rbbbs>0),'-m','LineWidth',2) % Boundary
        %   end
        
        hold off
    else % ----- PIXEL PLOT -----
%         A = transpose(squeeze(fastcam.(fastcam.moviToUse)(:,:,flags.plotBin(k))));
        
        A = transpose(squeeze(fastcam.(fastcam.moviToUse)(:,:,flags.plotBin(k)) - fastcam.moviBase));
        
        B = reshape(A,[numel(A),1]);
        B(B<=0) = [];
        minval = min(B);
        maxval = max(B);
        v = linspace(minval,maxval,30);
        
        NR = size(A,2);
        NZ = size(A,1);
        
        DZ = diff(raymap.majz,1,1);
        max_DZ = max(DZ,[],1);
        [~,I] = min(abs((max_DZ)));
        
        [~,J] = min(abs(raymap.majr(:,I)-G.rmaxis));
        
        [~,K] = min(abs(raymap.majr(:,I) - Riw));
        
        [~,L] = min(abs(raymap.majz(K,:)-G.zmaxis));
        
        [~,M] = min(abs(raymap.majz(K,:)));
              
        DR = (G.rmaxis-Riw)/(J-K+1);
        Ro = Riw-K*DR;
        Rf = Ro + NR*DR;
        rAxis = linspace(Ro,Rf,NR);       
        DZ = DR;

%         DZ = G.zmaxis/(L-M+1);
        Zo = -M*DZ;
        Zf = (NZ-M)*DZ;
        zAxis = linspace(Zo,Zf,NZ);

%         DR = DZ;
%         Ro = Riw-K*DR;
%         Rf = Ro + NR*DR;
%         rAxis = linspace(Ro,Rf,NR);
        
        I0 = find(zAxis>0,1);
        zo = zAxis(I0);
        
        I1 = find(zAxis>0.25,1);
        Z1 = zAxis(I1);
        
        I2 = find(zAxis>-0.25,1)-1;
        Z2 = zAxis(I2);
        
        I3 = find(zAxis>0.4,1);
        z3 = zAxis(I3);
        
        I4 = find(zAxis>-0.4,1)-1;
        z4 = zAxis(I4);

        figure(hf)
        subplot(4,length(in.ts),k+length(in.ts)*[0 1 2])
        title(strcat(num2str(in.shot),'.0',num2str(in.ts(k))),...
            'FontSize',flags.fontSize,'Color',lin.colors{k})
        hold on
        
        %         contourf(A,v,'LineStyle','none')
        contourf(rAxis,zAxis,A,v,'LineStyle','none')
        plot(rAxis,zo*ones(size(rAxis)),'-r',rAxis,Z1*ones(size(rAxis)),'-g',...
            rAxis,Z2*ones(size(rAxis)),'-b',rAxis,z3*ones(size(rAxis)),'-k',...
            rAxis,z4*ones(size(rAxis)),'-m','LineWidth',2)
        box on;grid minor;ax = gca;caxis(ax.CLim);cm = colormap(jet(1024));cm(1,:) = [1,1,1];colormap(cm);
        axis equal;xlim([0.95 2.45]);ylim(in.mapYlim)
        xlabel('$R$ (m)','FontSize',flags.fontSize,'Interpreter','latex')
        ylabel('$Z$ (m)','FontSize',flags.fontSize,'Interpreter','latex')
        
        C = max(A(I0,:));
        A0 = A(I0,:)/C;
        A1 = A(I1,:)/C;
        A2 = A(I2,:)/C;
        A3 = A(I3,:)/C;
        A4 = A(I4,:)/C;
        
        figure(hs);
        subplot(1,length(in.ts),k)
        plot(rAxis,A0,'r',rAxis,A1,'g',rAxis,A2,'b',...
            rAxis,A3,'k',rAxis,A4,'m','LineWidth',1)
        legend({['$Z=$' num2str(zo)],['$Z=$' num2str(Z1)],['$Z=$' num2str(Z2)],...
            ['$Z=$' num2str(z3)],['$Z=$' num2str(z4)]},'Interpreter','latex')
        title(strcat(num2str(in.shot),'.0',num2str(in.ts(k))),...
            'FontSize',flags.fontSize,'Color',lin.colors{k})
        xlabel('$R$ (m)','FontSize',flags.fontSize,'Interpreter','latex')
        axis([1 2.5 0 1])
        
        % if isfield(in,'EFIT') % EFIT OVERLAY
        
        % Define G
        if flags.evolveEFIT; G=Gfull.gdata(flags.GplotBin(k)); end
        
        figure(hf)
        plot(G.rmaxis,G.zmaxis,'sk','LineWidth',5)% Magnetic Axis
        % Overlay rational surfaces
        qToPlot=[1 1.5 2 3 4 5 6];
        G.psiLevels=linspace(G.psimag,G.psibry,length(G.qpsi));
        for jj=1:length(qToPlot)
            [~, qBins(jj)]=min(abs(G.qpsi-qToPlot(jj)));
        end
        %flags.cols2=colormap('copper');
        [G.Rgrid, G.Zgrid]=meshgrid(G.rg, G.zg);
        
        figure(hf)
        contour(G.Rgrid,G.Zgrid,G.psizr,...
            G.psiLevels(qBins),'--',...
            'LineWidth',1,'Color','k')
        
        figure(hf)
        plot(G.xlim, G.ylim,'-k','LineWidth',2) % LIMITER
        plot(G.rbbbs(G.rbbbs>0),G.zbbbs(G.rbbbs>0),'-m','LineWidth',2) % Boundary
        %   end
        
%         contour(transpose(squeeze(fastcam.(fastcam.moviToUse)(:,:,flags.plotBin(k)))),in.zContours)
        % Plot the box
%         for kk=1:length(in.box)
%             plot(in.box(kk).xSlice+in.box(kk).size*[1 1 -1 -1 1],...
%                 in.box(kk).ySlice+in.box(kk).size*[-1 1 1 -1 -1],...
%                 lin.style2{kk},'LineWidth',2)
%         end
%         xlim([0 in.picSize])
%         ylim([0 in.picSize])
%         
%         xlabel('X (pix)','FontSize',flags.fontSize)
%         if k==1;  ylabel('Y (pix)','FontSize',flags.fontSize); end
%         set(get(gcf,'CurrentAxes'),'FontSize',flags.fontSize)
%         grid on; box on

    hold off 
    end
end

% Integrate counts in box time trace
subplot(4,length(in.ts),length(in.ts)*3+[1:length(in.ts)])
hold on
for kk=1:length(in.box)
    plot(fastcam.time(movieBins)/1e3-in.tShift/1e3,...
        fastcam.int(kk).sig(movieBins)/max(fastcam.int(kk).sig(movieBins)),...
        lin.style2{kk},'LineWidth',2)
end
ydum=[0 1.2];
%ydum=get(gca,'Ylim');
for k=1:length(in.ts)
    plot(time(flags.plotBin(k))*ones(size(ydum)),ydum,lin.style3{k})
end
plot(in.tLM*[1 1]/1e3,ydum,'-r')
for k=1:2
    plot(in.baseline_ts(k)*[1 1]/1e3-in.tShift/1e3,ydum,'-m')
end

if flags.plotNegCounts==0;  ydum=[0 ydum(2)]; end
if flags.tShift; plot([0 0],ydum,'-k'); end
ylim(ydum);

xlabel(flags.xStr,'FontSize',flags.fontSize)
ylabel('box counts / max','FontSize',flags.fontSize)
set(get(gcf,'CurrentAxes'),'FontSize',flags.fontSize)
grid on; box on
xlim(flags.xLim-in.tShift/1e3)

%% Other code lines
% r = eq.oneD.rho(1,:);
% t = eq.oneD.t(:,1);

% figure;hold on;
% ne=mean(eq.oneD.raw.ne(53:end,:),1);ne=ne/max(ne(1));
% Te=mean(eq.oneD.raw.Te(53:end,:),1);Te=Te/max(Te(1));
% Zeff=mean(eq.oneD.final.Zeff(53:end,:),1);Zeff=Zeff/max(Zeff(1));
% plot(r,ne,'b',r,Te,'r',r,Zeff,'k');
% ylabel('$f(\rho)$','Interpreter','latex');
% xlabel('$\rho$','Interpreter','latex');


% yyaxis left;plot(r,mean(eq.oneD.raw.ne(53:end,:),1));xlabel('$\rho$','Interpreter','latex');ylabel('$n_e$','Interpreter','latex');yyaxis right;plot(r,mean(eq.oneD.raw.Te(53:end,:),1));;ylabel('$T_e$','Interpreter','latex');
% yyaxis left;plot(r,mean(eq.oneD.final.Zeff(53:end,:),1));xlabel('$\rho$','Interpreter','latex');ylabel('$Z_{eff}$','Interpreter','latex');yyaxis right;plot(r,mean(eq.oneD.final.E_over_ED(53:end-1,:),1));ylabel('$E/E_{CH}$','Interpreter','latex')
% 
% f=mean(eq.oneD.raw.ne(53:end,:),1);f=f/max(f);yyaxis left;plot(r,f);xlabel('$\rho$','Interpreter','latex');ylabel('$n_e$','Interpreter','latex');yyaxis right;f=mean(eq.oneD.raw.Te(53:end,:),1);f=f/max(f);plot(r,f);ylabel('$T_e$','Interpreter','latex');
% yyaxis left;plot(r,mean(eq.oneD.final.Zeff(53:end,:),1));xlabel('$\rho$','Interpreter','latex');ylabel('$Z_{eff}$','Interpreter','latex');yyaxis right;plot(r,mean(eq.oneD.final.E_over_ECH(53:end-1,:),1))
% 
% surf(eq.oneD.t,eq.oneD.rho,eq.oneD.final.E_over_ECH,'LineStyle','none');box on;colormap(jet);xlabel('$t$','Interpreter','latex');ylabel('$\rho$','Interpreter','latex')
% surf(eq.oneD.t,eq.oneD.rho,eq.oneD.raw.ne,'LineStyle','none');box on;colormap(jet);xlabel('$t$','Interpreter','latex');ylabel('$\rho$','Interpreter','latex')
% surf(eq.oneD.t,eq.oneD.rho,eq.oneD.raw.Te,'LineStyle','none');box on;colormap(jet);colormap(jet);xlabel('$t$','Interpreter','latex');ylabel('$\rho$','Interpreter','latex')
% surf(eq.oneD.t,eq.oneD.rho,log10(eq.oneD.calc.nu_ee),'LineStyle','none');box on;colormap(jet);

