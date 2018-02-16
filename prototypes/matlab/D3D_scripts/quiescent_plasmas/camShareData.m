close all
clear variables

%load('data/165826_camdata.mat')
load('165139_camdata.mat')

pixelPlot=0;

%% Pre-amble

for k=1:length(flags.xLim)
    [~, movieLims(k)]=min(abs(fastcam.time/1e3-flags.xLim(k)));
end
movieBins=movieLims(1):1:movieLims(2);


%% Make Simpler Profile Figure (WITH EFIT AND RAYMAP)


lin.style2={'-k','-r','-b','-g','-m','-y','-c'};
lin.style3={'--k','--r','--b','--g','--m','--y','--c'};
lin.colors={'k','r','b','g','m','y','c'};

figure('Position',[300 100 1300 600])

if flags.plotNegCounts
    colormap(darkb2r(in.zContourLim(1),in.zContourLim(2)));
end

% First Row
for k=1:length(in.ts)
    subplot(4,length(in.ts),k+length(in.ts)*[0 1 2])
    title(strcat(num2str(in.shot),'.0',num2str(in.ts(k))),...
        'FontSize',flags.fontSize,'Color',lin.colors{k})
    hold on
    if pixelPlot==0
        % if isfield(in,'EFIT') % EFIT OVERLAY
        
        % Define G
        if flags.evolveEFIT; G=Gfull.gdata(flags.GplotBin(k)); end
        
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
        
        plot(G.xlim, G.ylim,'-k','LineWidth',2) % LIMITER
        plot(G.rbbbs(G.rbbbs>0),G.zbbbs(G.rbbbs>0),'-k','LineWidth',5) % Boundary
        %   end
        
        contour(raymap.majr, raymap.majz,(squeeze(fastcam.(fastcam.moviToUse)(:,:,flags.plotBin(k)))),...
            in.zContours,'LineWidth',2)
        caxis([in.zContourLim])
        % Plot the box
        for kk=1:length(in.box)
            plot(raymap.majr(in.box(kk).xSlice+in.box(kk).size*[1 1 -1 -1 1],1),...
                raymap.majz(1,in.box(kk).ySlice+in.box(kk).size*[-1 1 1 -1 -1]),...
                lin.style2{kk},'LineWidth',2)
        end
        xlim([0.95 2.45])
        ylim(in.mapYlim)
        
        
        xlabel('R (m)','FontSize',flags.fontSize)
        if k==1;  ylabel('Z (m)','FontSize',flags.fontSize); end
        set(get(gcf,'CurrentAxes'),'FontSize',flags.fontSize)
        grid on; box on
        
    else % ----- PIXEL PLOT -----
        contour(transpose(squeeze(fastcam.(fastcam.moviToUse)(:,:,flags.plotBin(k)))),in.zContours)
        % Plot the box
        for kk=1:length(in.box)
            plot(in.box(kk).xSlice+in.box(kk).size*[1 1 -1 -1 1],...
                in.box(kk).ySlice+in.box(kk).size*[-1 1 1 -1 -1],...
                lin.style2{kk},'LineWidth',2)
        end
        xlim([0 in.picSize])
        ylim([0 in.picSize])
        
        xlabel('X (pix)','FontSize',flags.fontSize)
        if k==1;  ylabel('Y (pix)','FontSize',flags.fontSize); end
        set(get(gcf,'CurrentAxes'),'FontSize',flags.fontSize)
        grid on; box on
        
    end
    caxis([in.zContourLim])
    axis square
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
    plot(fastcam.time(flags.plotBin(k))*[1 1]/1e3-in.tShift/1e3,ydum,lin.style3{k})
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

