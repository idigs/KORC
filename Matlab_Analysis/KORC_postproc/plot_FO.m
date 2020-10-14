format longE

%% load data from file

%ST=diagnoseKORC('../','off',[0,101]);
load('ST.mat')

%ST=ST_EI_FOold;

%% extract data from ST structure

c=ST.params.scales.v;
me=ST.params.scales.m;
qe=ST.params.scales.q;

time=ST.time;

X=squeeze(ST.data.sp1.X(:,1,:));
Y=squeeze(ST.data.sp1.X(:,2,:));
Z=squeeze(ST.data.sp1.X(:,3,:));

[PHI,R,Z]=cart2pol(X,Y,Z);

BX=squeeze(ST.data.sp1.B(:,1,:));
BY=squeeze(ST.data.sp1.B(:,2,:));
BZ=squeeze(ST.data.sp1.B(:,3,:));

BMAG=sqrt(BX.^2+BY.^2+BZ.^2);

BR=BX.*cos(PHI)+BY.*sin(PHI);
BPHI=-BX.*sin(PHI)+BY.*cos(PHI) ; 

EX=squeeze(ST.data.sp1.E(:,1,:));
EY=squeeze(ST.data.sp1.E(:,2,:));

PSIp=squeeze(ST.data.sp1.PSIp(:,:));

EPHI=-EX.*sin(PHI)+EY.*cos(PHI)  ;

VX=squeeze(ST.data.sp1.V(:,1,:));
VY=squeeze(ST.data.sp1.V(:,2,:));
VZ=squeeze(ST.data.sp1.V(:,3,:));

VMAG=sqrt(VX.^2+VY.^2+VZ.^2);

gam=1./sqrt(1-(VMAG/c).^2);

ppll=gam*me./BMAG.*(VX.*BX+VY.*BY+VZ.*BZ);
mu=(gam.^2*me^2.*VMAG.^2-ppll.^2)./(2*me*BMAG);
  

pmag=sqrt(ppll.^2+mu.*BMAG*2*me);

%gam=ST.data.sp1.g(:,:);
eta=ST.data.sp1.eta(:,:);
xi=cos(deg2rad(eta));

vpll=ppll./(gam*me);

flag=ST.data.sp1.flag;

I=sum(qe*vpll.*flag.*BPHI./BMAG,1);

AveE=me*c^2/qe*sum(gam.*flag)./sum(flag);
AveE(isnan(AveE))=0;
AveKE=me*c^2/qe*sum((gam-1).*flag)./sum(flag);
AveKE(isnan(AveKE))=0;

%% construct background toroidal magnetic vector potential

PSIP=ST.params.fields.psi_p;
FLAG=ST.params.fields.Flag;

RF=ST.params.fields.R;
ZF=ST.params.fields.Z;

%% plot data
%close all

histpeta=1;
histreta=0;
histp=0;
histxi=0;
histRZ=1;
evoEave=1;
evoI=1;

timeslice=1;
timeslice=size(time,2);

%figure;
%hold on
%scatter(R(:,timeslice),Z(:,timeslice),'o')
%hold off
%daspect([1,1,1])


if histpeta==1
    fig=figure;
    hist=histogram2(eta(:,timeslice),pmag(:,timeslice)/(me*c),50,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    colormap jet
    colorbar()
    xlabel('Pitch Angle $\theta$','Interpreter','latex','FontSize',14)
    ylabel('$p$ ($m_ec$)','Interpreter','latex','FontSize',14)
%    axis([0,55,hist.YBinLimits(1),hist.YBinLimits(2)])
    hold off
    
    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'TEST13_V0','-dpdf')
end

if histreta==1
    fig=figure;
    hist=histogram2(R(:,timeslice),eta(:,timeslice),50,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    colormap jet
    colorbar()
    xlabel('$R (m)$','Interpreter','latex','FontSize',14)
    ylabel('Pitch Angle $\theta$','Interpreter','latex','FontSize',14)
%    axis([0,55,hist.YBinLimits(1),hist.YBinLimits(2)])
    hold off
    
    print(fig,'TEST13_V0','-dpdf')
end

if histp==1
    fig=figure;
    hist=histogram(pmag(:,timeslice)/(me*c),50,'Normalization','count');
    hold on
   
%    line([pmag(1,1)/(me*c),pmag(1,1)/(me*c)],[0,max(hist.Values)],...
%        'LineStyle','--','Linewidth',2.,'Color','red')

    xlabel('$p$ ($m_ec$)','Interpreter','latex','FontSize',14)
    ylabel('$f(p)$','Interpreter','latex','FontSize',14)
%    axis([0,21,0,max(hist.Values)])
    hold off
    
    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'TEST13_V0','-dpdf')
end

if histxi==1
    fig=figure;
    hist=histogram(xi(flag(:,timeslice)>0,timeslice),50,'Normalization','count');
    hold on
   
%    line([xi(1,1),xi(1,1)],[0,max(hist.Values)],...
%        'LineStyle','--','Linewidth',2.,'Color','red')

    xlabel('$\xi$','Interpreter','latex','FontSize',14)
    ylabel('$f(\xi)$','Interpreter','latex','FontSize',14)
%    axis([-1.1,1.1,0,inf])
    hold off
    
    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'TEST13_V0','-dpdf')
end

if histRZ==1
    fig=figure;
%    contour(RF,ZF,PSIP',50,'k','LineWidth',1,...
%        'LineStyle','--')
    hold on
%    contour(RF,ZF,FLAG',1,'k','LineWidth',2)
    hist=histogram2(R(:,timeslice),Z(:,timeslice),50,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    colormap jet
    colorbar()
    xlabel('${\rm R\,(m)}$','Interpreter','latex','FontSize',14)
    ylabel('${\rm Z\,(m)}$','Interpreter','latex','FontSize',14)
%    axis([1.4,1.85,-.175,.175])
    daspect([1,1,1])

    print(fig,'TEST13C_X0','-dpdf')
end

if evoEave==1
    fig=figure; 
    plot(time,AveE)
    xlabel('${\rm t\,(s)}$','Interpreter','latex','FontSize',14)
    ylabel('${\rm E_{\rm ave}\,(eV)}$','Interpreter','latex','FontSize',14)

    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'TEST13_evoEave','-dpdf')
    
if evoI==1
    fig=figure; 
    plot(time,I/I(1))
    xlabel('${\rm t\,(s)}$','Interpreter','latex','FontSize',14)
    ylabel('${\rm I/I(0)}$','Interpreter','latex','FontSize',14)

    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
         strcat('ne= ',num2str(ST.params.profiles.neo))};
    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    print(fig,'evoI','-dpdf')

end

end