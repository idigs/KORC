format longE
set(0,'defaultfigurecolor',[1 1 1])

%% load data from file

join=0;
if(join==1)
    
    load('../TEST36c/36c/ST.mat')
    ST36c=ST;
    load('../TEST36c/36ca/ST.mat')
    ST36ca=ST;
    
    ST.time=cat(2,ST36c.time,ST36ca.time);   
    ST.data.sp1.X=cat(3,ST36c.data.sp1.X,ST36ca.data.sp1.X);
    ST.data.sp1.V=cat(3,ST36c.data.sp1.V,ST36ca.data.sp1.V);
    ST.data.sp1.B=cat(3,ST36c.data.sp1.B,ST36ca.data.sp1.B);
    ST.data.sp1.E=cat(3,ST36c.data.sp1.E,ST36ca.data.sp1.E);
    ST.data.sp1.g=cat(2,ST36c.data.sp1.g,ST36ca.data.sp1.g);
    ST.data.sp1.eta=cat(2,ST36c.data.sp1.eta,ST36ca.data.sp1.eta);
    ST.data.sp1.flagCon=cat(2,ST36c.data.sp1.flagCon,ST36ca.data.sp1.flagCon);
    ST.data.sp1.flagCol=cat(2,ST36c.data.sp1.flagCol,ST36ca.data.sp1.flagCol);

else
    
    load('ST.mat')
    
end

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

PHI(PHI<0)=pi+(pi-abs(PHI(PHI<0))); 

R0F=2.604489052000000e+00; % magnetic axis as determined by EFIT
Z0F=-1.353408685000000e-01;

THETA=atan2(Z-Z0F,R-R0F);
THETA(THETA<0)=pi+(pi-abs(THETA(THETA<0))); 


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

flagCon=ST.data.sp1.flagCon;
flagCol=ST.data.sp1.flagCol;

active=flagCon.*flagCol;
lost=zeros(size(active));
lost(active<1)=1;

conflag=zeros(size(flagCon));
conflag(flagCon<1)=1;

conlossind=zeros(size(conflag));
for i=2:size(conlossind,2)
    for j=1:size(conlossind,1)
        if (conflag(j,i)-conflag(j,i-1)==1)
            conlossind(j,i)=1;
        end
    end
end

colflag=zeros(size(flagCol));
colflag(flagCol<1)=1;

Ipart=qe*vpll.*BPHI./(2*pi*R.*BMAG);
Ipart(active==0)=0;
I=sum(Ipart,1);


%gam(flag==0)=0;
E=me*c^2/qe*sum(gam.*active);
KE=me*c^2/qe*sum((gam-1).*active);

AveE=me*c^2/qe*sum(gam.*active)/size(R,1);
AveE(isnan(AveE))=0;
AveKE=me*c^2/qe*sum((gam-1).*active)/size(R,1);
AveKE(isnan(AveKE))=0;

%% construct background toroidal magnetic vector potential

PSIP=ST.params.fields.psi_p;
FLAG=ST.params.fields.Flag2D;

RF=ST.params.fields.R;
ZF=ST.params.fields.Z;

%% LOAD AORSA DATA

filename='/Users/21b/Desktop/AORSA/171089/AORSA_Data/bharvey-200MHz';

fileID = fopen(filename,'r');
formatSpec1 = '%u';
formatSpec2 = '%f';
formatSpec3 = '%f %f';

init = textscan(fileID,formatSpec1,2);
NRa=init{1}(1);
NZa=init{1}(2);


A = textscan(fileID,formatSpec2,NRa);
B = textscan(fileID,formatSpec2,NZa);
C = textscan(fileID,formatSpec3,NRa*NZa);
D = textscan(fileID,formatSpec3,NRa*NZa);
E = textscan(fileID,formatSpec3,NRa*NZa);
F = textscan(fileID,formatSpec3,NRa*NZa);
G = textscan(fileID,formatSpec3,NRa*NZa);
H = textscan(fileID,formatSpec3,NRa*NZa);
I = textscan(fileID,formatSpec3,NRa*NZa);
J = textscan(fileID,formatSpec3,NRa*NZa);
K = textscan(fileID,formatSpec3,NRa*NZa);
L = textscan(fileID,formatSpec2,1);
M = textscan(fileID,formatSpec2,NRa*NZa);

fclose(fileID);

%% MANIPULATE AORSA DATA

Ra=A{1};
Za=B{1};

[ZZa,RRa]=meshgrid(Za,Ra);

ReEX=reshape(F{1},[NZa,NRa]);
ImEX=reshape(F{2},[NZa,NRa]);
ReEY=reshape(G{1},[NZa,NRa]);
ImEY=reshape(G{2},[NZa,NRa]);
ReEZ=reshape(H{1},[NZa,NRa]);
ImEZ=reshape(H{2},[NZa,NRa]);

ReBX=reshape(I{1},[NZa,NRa]);
ImBX=reshape(I{2},[NZa,NRa]);
ReBY=reshape(J{1},[NZa,NRa]);
ImBY=reshape(J{2},[NZa,NRa]);
ReBZ=reshape(K{1},[NZa,NRa]);
ImBZ=reshape(K{2},[NZa,NRa]);

phi=0;
nn=35;

omega=2*pi*200*10^6;
atime=0;

BX1=real((ReBX+ImBX*1j)*exp(1j*(omega*atime+nn*phi)));
BY1=real((ReBY+ImBY*1j)*exp(1j*(omega*atime+nn*phi)));
BZ1=real((ReBZ+ImBZ*1j)*exp(1j*(omega*atime+nn*phi)));

EX1=real((ReEX+ImEX*1j)*exp(1j*(omega*atime+nn*phi)));
EY1=real((ReEY+ImEY*1j)*exp(1j*(omega*atime+nn*phi)));
EZ1=real((ReEZ+ImEZ*1j)*exp(1j*(omega*atime+nn*phi)));

%% LOAD EFIT DATA

filename='/Users/21b/Desktop/AORSA/171089/EquilibriumData/g171089.05450';

fileID = fopen(filename,'r');
formatSpec1 = '%f';
formatSpec4 = '%f %f %f %f';

init = textscan(fileID,'%s %s %s %s %d %d %d',1);

tmp=split(init{3},"#");
pulse=tmp{2};
tmp1=split(init{4},"m");
etime=str2double(tmp1{1})/1e3;

NRe=init{end-1};
NZe=init{end};

A = textscan(fileID,formatSpec1,149);

Rdim=A{1}(1); %Distance from inner edge to outer edge covered by EFIT, in [m]
Zdim=A{1}(2); %Distance from bottom to top (Z_axis) covered by EFIT, equally spaced around midplane, in [m]
R0=A{1}(3); % major radius [m]
R1=A{1}(4); % Inner edge [m]
Zmid=A{1}(5); % Midplane [m]
RmAxis=A{1}(6); % R of magnetic axis [m]
ZmAxis=A{1}(7); % Z of magnetic axis [m]
psiAxis=A{1}(8); % poloidal flux at magnetic axis [Wb]
psiSep=A{1}(9); % poloidal flux at separatrix [Wb]
Bt0=A{1}(10); % toroidal magnetic field at major radius in [T]
Ip=A{1}(11); % plasma current in [A]

A = textscan(fileID,formatSpec1,NRe);
tmp=A{1};

A = textscan(fileID,formatSpec1,NRe);
tmp1=A{1};

A = textscan(fileID,formatSpec1,NRe);
tmp2=A{1};

B = textscan(fileID,formatSpec1,NRe*NZe); % PSIP

A = textscan(fileID,formatSpec1,NRe);
tmp3=A{1};

A = textscan(fileID,formatSpec1,2);
Nlcfs=A{1}(1);
Nwall=A{1}(2);

C = textscan(fileID,formatSpec1,Nlcfs*2); % LCFS
D = textscan(fileID,formatSpec1,Nwall*2); % Limiter

fclose(fileID);

%% MANIPULATE EFIT DATA

RR=linspace(R1,R1+Rdim,NRe);
ZZ=linspace(Zmid-Zdim/2,Zmid+Zdim/2,NZe);

PSIP=reshape(B{1},[NZe,NRe]);
PSIP1=interp2(ZZ,RR,PSIP,ZZa,RRa);

tmp4=zeros(size(Nlcfs,1));
tmp5=zeros(size(Nlcfs,1));

for i=1:1:Nlcfs*2
    if mod(i,2)==1
        tmp4(int16(i/2))=C{1}(i);
    else
        tmp5(int16(i/2))=C{1}(i);
    end        
end

RLCFS=tmp4;
ZLCFS=tmp5;

tmp6=zeros(size(Nlcfs,1));
tmp7=zeros(size(Nlcfs,1));

for i=1:1:Nwall*2
    if mod(i,2)==1
        tmp6(int16(i/2))=D{1}(i);
    else
        tmp7(int16(i/2))=D{1}(i);
    end
end

Rlim=tmp6;
Zlim=tmp7;

[dPSIdZ,dPSIdR]=gradient(PSIP1,Za(2)-Za(1),Ra(2)-Ra(1));

BRe=-1./RRa.*dPSIdZ;
BZe=1./RRa.*dPSIdR;
BPHIe=zeros(size(BRe));
for ii=1:size(Ra,1)
    BPHIe(ii,:)=R0*Bt0/Ra(ii);
end

%% plot data
%close all

histpeta=0;
histreta=0;
histp=0;
histxi=0;
histRZ=0;
evoEave=0;
evoI=0;
RZmovie=0;
histRevo=0;
compconfine=0;
deconmovie=0;
decon2D=0;
deconwallpower=0;
histEZ=0;
deconloc=0;
losstrajectory=0;

%timeslice=1;
timeslice=size(time,2);

%figure;
%hold on
%scatter(R(:,timeslice),Z(:,timeslice),'o')
%hold off
%daspect([1,1,1])

if losstrajectory==1
    fig=figure;
    box on
    
    ipart=1;
    
    %plot(R(ipart,:),Z(ipart,:)) 
    plot(R,Z) 
    hold on
    contour(RF,ZF,FLAG',1,'k','LineWidth',2)
    plot(RLCFS,ZLCFS,'k--','linewidth',2)
    
    daspect([1,1,1])

    ylabel('$Z\,({\rm m})$','interpreter','latex','fontsize',16)
    xlabel('$R\,({\rm m})$','interpreter','latex','fontsize',16)
    
    set(gca,'Fontsize',16)
    
    print(fig,'trajectory_TEST36c','-depsc')
end

if compconfine==1
   
    plotdeltaB=1;
    plotdeltaB_IAEA=1;
    plotE_deltaB200=0;
    plotE_deltaB1000=0;
    ploteta=0;

    fig=figure;
    box on
    hold on
    
    if plotdeltaB_IAEA==1
        
        load('ST_36_c_analysis.mat')
        activec=ST_36_c_analysis.active;
        conflagc=ST_36_c_analysis.conflag;
        timec=ST_36_c_analysis.time;       

        load('ST_36_e_analysis.mat')
        activee=ST_36_e_analysis.active;
        conflage=ST_36_e_analysis.conflag;
        timee=ST_36_e_analysis.time;
        
        load('ST_36_e1_analysis.mat')
        activee1=ST_36_e1_analysis.active;
        conflage1=ST_36_e1_analysis.conflag;
        timee1=ST_36_e1_analysis.time;
        
        load('ST_36_e2_analysis.mat')
        activee2=ST_36_e2_analysis.active;
        conflage2=ST_36_e2_analysis.conflag;
        timee2=ST_36_e2_analysis.time;
        
        load('ST_36_e3_analysis.mat')
        activee3=ST_36_e3_analysis.active;
        conflage3=ST_36_e3_analysis.conflag;
        timee3=ST_36_e3_analysis.time;
        
        p2=plot(timec*10^6,sum(conflagc)/sum(activec(:,1))*100,'-','linewidth',3,'color',[0,0,1]);
        p4=plot(timee2*10^6,sum(conflage2)/sum(activee2(:,1))*100,'-','linewidth',3,'color','k');
        p5=plot(timee1*10^6,sum(conflage1)/sum(activee1(:,1))*100,'-','linewidth',3,'color',[.5,0,.5]);
        p6=plot(timee*10^6,sum(conflage)/sum(activee(:,1))*100,'-','linewidth',3,'color',[0,.5,0]);
        p7=plot(timee3*10^6,sum(conflage3)/sum(activee3(:,1))*100,'-','linewidth',3,'color',[1,0,0]);
        
        legend([p7,p2,p4,p5,p6],{'$5{\rm G}$','$10{\rm G}$',...
            '$20{\rm G}$','$50{\rm G}$','$100{\rm G}$'},...
            'Interpreter','latex','FontSize',20,'position',[.765,.24,.1,.3])
        
        title('$E_0=10{\rm MeV}, \eta_0=170^\circ$','Interpreter','latex')
        
    elseif plotdeltaB==1
        load('ST_36_b_analysis.mat')
        activeb=ST_36_b_analysis.active;
        conflagb=ST_36_b_analysis.conflag;
        timeb=ST_36_b_analysis.time;
        
        load('ST_36_c_analysis.mat')
        activec=ST_36_c_analysis.active;
        conflagc=ST_36_c_analysis.conflag;
        timec=ST_36_c_analysis.time;
        
        load('ST_36_d_analysis.mat')
        actived=ST_36_d_analysis.active;
        conflagd=ST_36_d_analysis.conflag;
        timed=ST_36_d_analysis.time;

        load('ST_36_e_analysis.mat')
        activee=ST_36_e_analysis.active;
        conflage=ST_36_e_analysis.conflag;
        timee=ST_36_e_analysis.time;
        
        load('ST_36_e1_analysis.mat')
        activee1=ST_36_e1_analysis.active;
        conflage1=ST_36_e1_analysis.conflag;
        timee1=ST_36_e1_analysis.time;
        
        load('ST_36_e2_analysis.mat')
        activee2=ST_36_e2_analysis.active;
        conflage2=ST_36_e2_analysis.conflag;
        timee2=ST_36_e2_analysis.time;
        
        load('ST_36_e3_analysis.mat')
        activee3=ST_36_e3_analysis.active;
        conflage3=ST_36_e3_analysis.conflag;
        timee3=ST_36_e3_analysis.time;
        
        p1=plot(timeb,sum(conflagb)/sum(activeb(:,1)),'-','linewidth',3,'color',[0,0,1]);
        p2=plot(timec,sum(conflagc)/sum(activec(:,1)),'-','linewidth',3,'color',[1,0,0]);
        p3=plot(timed,sum(conflagd)/sum(actived(:,1)),'-','linewidth',3,'color',[0,.5,0]);
        p4=plot(timee2,sum(conflage2)/sum(activee2(:,1)),'-','linewidth',3,'color','k');
        p5=plot(timee1,sum(conflage1)/sum(activee1(:,1)),'-','linewidth',3,'color',[.5,0,.5]);
        p6=plot(timee,sum(conflage)/sum(activee(:,1)),'-','linewidth',3,'color',[0,1,1]);
        p7=plot(timee3,sum(conflage3)/sum(activee3(:,1)),'-','linewidth',3,'color',[1,1,0]);
        
        legend([p3,p1,p7,p2,p4,p5,p6],{'$\delta B=0$','$\delta B=10{\rm G}$','$\delta B=50{\rm G}$','$\delta B=100{\rm G}$',...
            '$\delta B=200{\rm G}$','$\delta B=500{\rm G}$','$\delta B=1000{\rm G}$'},...
            'Interpreter','latex','FontSize',20,'location','southeast')
        
        title('$E_0=10{\rm MeV}, \eta_0=170^\circ$','Interpreter','latex')
    
    elseif plotE_deltaB200==1
        load('ST_36_e2_analysis.mat')
        activee2=ST_36_e2_analysis.active;
        conflage2=ST_36_e2_analysis.conflag;
        timee2=ST_36_e2_analysis.time;

        load('ST_36_f_analysis.mat')
        activef=ST_36_f_analysis.active;
        conflagf=ST_36_f_analysis.conflag;
        timef=ST_36_f_analysis.time;

        load('ST_36_f1_analysis.mat')
        activef1=ST_36_f1_analysis.active;
        conflagf1=ST_36_f1_analysis.conflag;
        timef1=ST_36_f1_analysis.time;
                
        load('ST_36_f2_analysis.mat')
        activef2=ST_36_f2_analysis.active;
        conflagf2=ST_36_f2_analysis.conflag;
        timef2=ST_36_f2_analysis.time;
        
        load('ST_36_f3_analysis.mat')
        activef3=ST_36_f3_analysis.active;
        conflagf3=ST_36_f3_analysis.conflag;
        timef3=ST_36_f3_analysis.time;        
        
        load('ST_36_f4_analysis.mat')
        activef4=ST_36_f4_analysis.active;
        conflagf4=ST_36_f4_analysis.conflag;
        timef4=ST_36_f4_analysis.time;
        
        p4=plot(timee2,sum(conflage2)/sum(activee2(:,1)),'-','linewidth',3,'color','k');
        p7=plot(timef,sum(conflagf)/sum(activef(:,1)),'-','linewidth',3,'color','[0,0,1]');
        p8=plot(timef1,sum(conflagf1)/sum(activef1(:,1)),'-','linewidth',3,'color',[1,0,0]);
        p9=plot(timef2,sum(conflagf2)/sum(activef2(:,1)),'-','linewidth',3,'color',[0,.5,0]);
        p9a=plot(timef3,sum(conflagf3)/sum(activef3(:,1)),'-','linewidth',3,'color',[.5,0,.5]);
        p9b=plot(timef4,sum(conflagf4)/sum(activef4(:,1)),'-','linewidth',3,'color',[0,1,1]);
        
        legend([p7,p4,p8,p9b,p9a,p9],{'$E_0=5{\rm MeV}$','$E_0=10{\rm MeV}$','$E_0=20{\rm MeV}$',...
        '$E_0=25{\rm MeV}$','$E_0=30{\rm MeV}$','$E_0=40{\rm MeV}$',},...
        'Interpreter','latex','FontSize',20,'location','southeast')
    
        title('$\delta B=200{\rm G}, \eta_0=170^\circ$','Interpreter','latex')
        
    elseif plotE_deltaB1000==1
        load('ST_36_e_analysis.mat')
        activee=ST_36_e_analysis.active;
        conflage=ST_36_e_analysis.conflag;
        timee=ST_36_e_analysis.time;
        
        load('ST_36_h_analysis.mat')
        activeh=ST_36_h_analysis.active;
        conflagh=ST_36_h_analysis.conflag;
        timeh=ST_36_h_analysis.time;
        
        load('ST_36_h1_analysis.mat')
        activeh1=ST_36_h1_analysis.active;
        conflagh1=ST_36_h1_analysis.conflag;
        timeh1=ST_36_h1_analysis.time;
        
        load('ST_36_h2_analysis.mat')
        activeh2=ST_36_h2_analysis.active;
        conflagh2=ST_36_h2_analysis.conflag;
        timeh2=ST_36_h2_analysis.time;
        
        load('ST_36_h3_analysis.mat')
        activeh3=ST_36_h3_analysis.active;
        conflagh3=ST_36_h3_analysis.conflag;
        timeh3=ST_36_h3_analysis.time;
        
        load('ST_36_h4_analysis.mat')
        activeh4=ST_36_h4_analysis.active;
        conflagh4=ST_36_h4_analysis.conflag;
        timeh4=ST_36_h4_analysis.time;
        
        load('ST_36_h5_analysis.mat')
        activeh5=ST_36_h5_analysis.active;
        conflagh5=ST_36_h5_analysis.conflag;
        timeh5=ST_36_h5_analysis.time;
        
        p6=plot(timee,sum(conflage)/sum(activee(:,1)),'-','linewidth',3,'color',[0,0,0]);
        p13=plot(timeh,sum(conflagh)/sum(activeh(:,1)),'-','linewidth',3,'color',[0,0,1]);
        p14=plot(timeh1,sum(conflagh1)/sum(activeh1(:,1)),'-','linewidth',3,'color',[1,0,0]);
        p15=plot(timeh2,sum(conflagh2)/sum(activeh2(:,1)),'-','linewidth',3,'color',[0,.5,0]);
        p16=plot(timeh3,sum(conflagh3)/sum(activeh3(:,1)),'-','linewidth',3,'color',[.5,0,.5]);
        p17=plot(timeh4,sum(conflagh4)/sum(activeh4(:,1)),'-','linewidth',3,'color',[0,1,1]);
        p18=plot(timeh5,sum(conflagh5)/sum(activeh5(:,1)),'-','linewidth',3,'color',[1,1,0]);
        
        legend([p13,p6,p14,p15,p16,p17,p18],{'$E_0=5{\rm MeV}$','$E_0=10{\rm MeV}$','$E_0=20{\rm MeV}$',...
            '$E_0=30{\rm MeV}$','$E_0=40{\rm MeV}$','$E_0=60{\rm MeV}$','$E_0=80{\rm MeV}$'},...
            'Interpreter','latex','FontSize',20,'location','southeast')
        
        title('$\delta B=1000{\rm G}, \eta_0=170^\circ$','Interpreter','latex')
        
    elseif ploteta==1
        load('ST_36_e2_analysis.mat')
        activee2=ST_36_e2_analysis.active;
        conflage2=ST_36_e2_analysis.conflag;
        timee2=ST_36_e2_analysis.time;
        
        p4=plot(timee2,sum(conflage2)/sum(activee2(:,1)),'-','linewidth',3,'color','k');
        p10=plot(timeg,sum(conflagg)/sum(activeg(:,1)),'-','linewidth',3,'color','[0,0,1]');
        p11=plot(timeg1,sum(conflagg1)/sum(activeg1(:,1)),'-','linewidth',3,'color',[1,0,0]);
        p12=plot(timeg2,sum(conflagg2)/sum(activeg2(:,1)),'-','linewidth',3,'color',[0,.5,0]);
        
        legend([p10,p4,p11,p12],{'$\eta_0=179^\circ$','$\eta_0=170^\circ$','$\eta_0=135^\circ$','$\eta_0=95^\circ$',},...
            'Interpreter','latex','FontSize',20,'location','east')
        
        title('$\delta B=100{\rm G}, E_0=10{\rm MeV}$','Interpreter','latex')
    end

     
    xlabel('${\rm t\,(\mu s)}$','Interpreter','latex','FontSize',20)  
    ylabel('RE Fraction Deconfined','FontSize',20,'Interpreter','latex')

    
    set(gca,'Fontsize',20)
    
    if plotdeltaB_IAEA==1
        print(fig,'compconfine_JET_95135_deltaB_IAEA','-dpdf')
    elseif plotdeltaB==1
        print(fig,'compconfine_JET_95135_deltaB','-depsc')
    elseif plotE_deltaB200==1
        print(fig,'compconfine_JET_95135_E_deltaB200','-depsc')
    elseif plotE_deltaB1000==1
        print(fig,'compconfine_JET_95135_E_deltaB1000','-depsc')
    elseif ploteta==1
        print(fig,'compconfine_JET_95135_eta','-depsc')
    end
    
end

if histRevo==1
    fig=figure;
    
    Redges=linspace(1.8,3.9,50);
    Rbins=zeros(size(Redges,2)-1,1);
    for ii=1:size(Rbins,1)
        Rbins(ii)=(Redges(ii)+Redges(ii+1))/2;
    end
    
    [hist1,~]=histcounts(R(active(:,1)>0,1),Redges,'Normalization','pdf'); 
    [hist2,~]=histcounts(R(active(:,end)>0,end),Redges,'Normalization','pdf'); 
    
    NRE0=sum(flagCon(:,1));
    NRE1=sum(flagCon(:,end));
    
    hold on
    box on
    
    plot(Rbins,hist1,'linewidth',3,'color',[0,0,1])
    plot(Rbins,hist2*NRE1/NRE0,'linewidth',3,'color',[1,0,0])
    plot(Rbins,hist2*NRE1/NRE0-hist1,'k--','linewidth',3)
    
    xlabel('$R({\rm m})$','Interpreter','latex')
    ylabel('$Rf_R({\rm m^{-1}})$','Interpreter','latex')
    title('$\eta_0=170^\circ$','Interpreter','latex')
    
    legtxt={strcat('$t=$',num2str(time(1))),strcat('$t=$',num2str(time(end))),...
        '$\Delta$'};
    
    xlim([min(Rbins),max(Rbins)])
    ylim([-.75,1.3])
    
    legend(legtxt,'Interpreter','latex','location','northeast','fontsize',16)
    
    set(gca,'Fontsize',16)
    
    print(fig,'histRevo_TEST36c','-depsc')
    
end

if RZmovie==1
    
    loops = size(time,2);
    F(loops) = struct('cdata',[],'colormap',[]);
    
    figure('position',[10,10,425,500]);
    vidfile = VideoWriter('RZmovie_TEST36c.mp4','MPEG-4');
    open(vidfile);

    Rbins=linspace(1.8,3.9,50);
    Zbins=linspace(-1.6,2,50);
    
    maxhist=0;
    for ii=1:size(time,2)
        [hist,~]=histcounts2(R(active(:,ii)>0,ii),...
            Z(active(:,ii)>0,ii),Rbins,Zbins,'Normalization','count');        
        
        maxhisttmp=max(hist(:));
        maxhist=max(maxhist,maxhisttmp);        
    end
    
    colormap(plasma)
    
    for ii=1:size(time,2)
        hist=histogram2(R(active(:,ii)>0,ii),...
            Z(active(:,ii)>0,ii),Rbins,Zbins,'Normalization','count',...
            'DisplayStyle','tile','LineStyle','none',...
            'Normalization','count');               
        hold on
        contour(RF,ZF,FLAG',1,'k','LineWidth',2)
        hold off
    
        format shortG
        title(strcat('t=',num2str(time(ii),'%1.2e'),'s'))
        ylabel('$Z\,({\rm m})$','interpreter','latex','fontsize',16)
        xlabel('$R\,({\rm m})$','interpreter','latex','fontsize',16)

        daspect([1,1,1])
        set(gca,'Fontsize',16)
        

        cb=colorbar('position',[.85,.11,.03,.815]);
        caxis([0,maxhist]);
        cb.Label.String='$N_{\rm RE}$';
        cb.Label.Interpreter='latex';
        cb.Label.FontSize=16;
        
        colormap(plasma)
        set(gcf,'color','w');
        
        drawnow
        F(ii)=getframe(gcf);
        writeVideo(vidfile,F(ii));
    end
    
    close(vidfile)
end

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
    
    print(fig,'histpeta_TEST36c','-dpdf')
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
    
    print(fig,'histreta_TEST36c','-dpdf')
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
    
    print(fig,'histp_TEST36c','-dpdf')
end

if histxi==1
    fig=figure;
    hist=histogram(xi(flagCon(:,timeslice)>0,timeslice),50,'Normalization','count');
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
    
    print(fig,'histxi_TEST36c','-dpdf')
end

if histRZ==1
    fig=figure('position',[10,10,425,500]);
    box on
    hold on
    
    Rbins=linspace(.95,2.35,50);
    Zbins=linspace(-1.,1.25,50);
    
    maxc=35;
    maxe=33;
    maxe1=32;
    maxe2=34;
    
    maxh=28;
    maxh1=23;
    maxh2=6;
    maxh3=1;
    maxh4=4;
    maxh5=1;
    
    maxf=35;
    maxf1=30;
    maxf2=35;
    maxf3=35;
    maxf4=38;
    
    %contour(RF,ZF,FLAG',1,'k','LineWidth',2)
    hist=histogram2(R(active(:,timeslice)>0,timeslice),Z(active(:,timeslice)>0,timeslice),Rbins,Zbins,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    plot(RLCFS,ZLCFS,'k--','linewidth',2)
    plot(Rlim,Zlim,'k','linewidth',2)
    colormap plasma
    
    cb=colorbar('position',[.8,.11,.04,.815]);
    %caxis([0,maxf4]);
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=16;
        
    xlim([.9,2.45])
    
    xlabel('${\rm R\,(m)}$','Interpreter','latex','FontSize',14)
    ylabel('${\rm Z\,(m)}$','Interpreter','latex','FontSize',14)
%    axis([1.4,1.85,-.175,.175])
    daspect([1,1,1])

    %title('$\delta B=200{\rm G}$, $E_0=5{\rm MeV}$, $\eta_0=170^\circ$','interpreter','latex')
    
    set(gca,'Fontsize',16)
    
    print(fig,'histRZ_TEST36c','-depsc')
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
    
    print(fig,'evoEave_TEST36c','-dpdf')
    
end
    
if evoI==1
    
    fig=figure;
    box on
    hold on
    
    p1=plot(time,I/I(1),'linewidth',3,'color','k');
    p2=plot(time,E/E(1),'linewidth',3,'color',[0,0,1]);
    
    xlabel('${\rm t\,(s)}$','Interpreter','latex','FontSize',20)
    ylabel('Normalized Current and Energy','FontSize',20,'Interpreter','latex')
    
    yyaxis right
    HFS=R<2.5;
    LFS=R>2.5;
    
    set(gca,'Fontsize',20)
    
    startfrom=1;
    
    p4=plot(time(startfrom:end),(sum(HFS(:,startfrom:end).*conflag(:,startfrom:end))-sum(HFS(:,startfrom).*conflag(:,startfrom)))/...
        sum(active(:,startfrom)),'-','linewidth',3,'color',[1,0,0]);
    p5=plot(time(startfrom:end),(sum(LFS(:,startfrom:end).*conflag(:,startfrom:end))-sum(LFS(:,startfrom).*conflag(:,startfrom)))/...
        sum(active(:,startfrom)),'-','linewidth',3,'color',[.5,0,.5]);
    
    ylabel('RE Fraction','FontSize',20,'Interpreter','latex')
    set(gca,'ycolor','k')
    yyaxis left

    legend([p1,p2],{'KORC $I_{\rm RE}$','KORC $\mathcal{E}_{\rm RE}$'},...
        'Interpreter','latex','FontSize',16,'location','west')
    
    ah1=axes('position',get(gca,'position'),'visible','off');
    legend(ah1,[p4,p5],{'Lost to HFS','Lost to LFS'},...
        'Location','east','Interpreter','latex','FontSize',16)
    
%    txt={strcat('E:',num2str(ST.params.fields.Eo)),...
%         strcat('Impurity:',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('Z0= ',num2str(ST.params.collisions_ms.Zj)),...
%         strcat('dt= ',num2str(ST.params.simulation.dt)),...
%         strcat('ne= ',num2str(ST.params.profiles.neo))};
%    text(.25,.9,txt,'FontSize',14,'Units','normalized');
    
    set(gca,'Fontsize',20)
    
    print(fig,'evol_combo_TEST36c','-depsc')

end

if histEZ==1
    fig=figure;
    
    Ebinedges=linspace(0,max(gam(:)),55);
    Zbinedges=linspace(min(Z(:)),max(Z(:)),50);
    
    hist=histogram2(gam(conflag(:,timeslice)>0,timeslice),...
        Z(conflag(:,timeslice)>0,timeslice),Ebinedges,Zbinedges,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    
    set(gca,'Fontsize',12)
    
    title('JET 95135: Deconfined E VS Z','fontsize',16)
    ylabel('$Z\,({\rm m})$','interpreter','latex','fontsize',16)
    xlabel('$E/(m_ec^2)$','interpreter','latex','fontsize',16)
    
    cb=colorbar();
    cb.Label.String='$N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=16;
    
    print(fig,'JET_95135_REenergyZ_TEST36c','-dpdf')
    
    
    Zvals=Zbinedges(1:end-1)+hist.BinWidth(2)/2;
    Evals=Ebinedges(1:end-1)+hist.BinWidth(1)/2;
    
    intE=zeros(hist.NumBins(2),1);
    for ii=1:hist.NumBins(2)
        for jj=1:hist.NumBins(1)
            intE(ii)=intE(ii)+hist.Values(jj,ii)*Evals(jj);
        end
    end
    
    flagRZ=contourc(RF,ZF,double(FLAG)');
    Rin=flagRZ(1,end-1);
    
    intE=intE/(2*pi*Rin);
    
    fig=figure;
    
    semilogx(intE,Zvals,'linewidth',2)
    
    set(gca,'Fontsize',12)
    
    title('JET 95135: RE Energy per unit toroidal length deposited along Z','fontsize',16)
    ylabel('$Z\,({\rm m})$','interpreter','latex','fontsize',16)
    xlabel('$E/(2\pi R_{\rm lim}m_ec^2)\,[{\rm m^{-1}}]$','interpreter','latex','fontsize',16)
    
    print(fig,'JET_95135_REintenergyZ_TEST36c','-dpdf')
    
end

if deconwallpower==1
    
    ExpREs=2.813e5*2*pi*1.682/(c*qe);
    
    ID3D0=2.813e5;
    IKORC0=-I(1);
    
    RErat=ExpREs/size(R,1);
    RErat1=ID3D0/IKORC0;
    
    Eloss=conlossind.*gam*me*c^2*RErat1;

    maxZ=0.3;
    minZ=-0.6;
    
    bins=[5,10,15,20,25,30,32,35,37,40,42,45,47,50,52,55,57,60,62,65,67,70,75,80,85,90,95,100,...
        110,120,130,140,150,160];
    ZEpowerfluxmax=zeros(size(bins));
    
    for kk=1:size(bins,2)
        Zbinedges=linspace(minZ,maxZ,bins(kk)+1);
        dZ=Zbinedges(2)-Zbinedges(1);
        Zvals=zeros(1,size(Zbinedges,2)-1);
        for i=1:size(Zvals,2)
            Zvals(i)=(Zbinedges(i)+Zbinedges(i+1))/2;
        end

        ZEloss=zeros(size(Zbinedges,2)-1,size(Eloss,2));

        for i=1:size(Eloss,2)
            for j=1:size(Eloss,1)
                if Eloss(j,i)>0 && Z(j,i)<maxZ && Z(j,i)>minZ
                    Zind=int8(floor((Z(j,i)-Zbinedges(1))/dZ))+int8(1);
                    ZEloss(Zind,i)=ZEloss(Zind,i)+Eloss(j,i);
                end
            end       
        end

        ZEpower=ZEloss/(time(2)-time(1));

        flagRZ=contourc(RF,ZF,double(FLAG)');
        Rin=flagRZ(1,end-1);

        ZEpowerflux=ZEpower/(2*pi*Rin*dZ);

        ZEpowerfluxlevels=linspace(min(ZEpowerflux(:)),max(ZEpowerflux(:)),25)/10^6;
        
        ZEpowerfluxmax(kk)=max(ZEpowerflux(:));
    end
        

    Zbinedges=linspace(minZ,maxZ,70+1);
    dZ=Zbinedges(2)-Zbinedges(1);
    Zvals=zeros(1,size(Zbinedges,2)-1);
    for i=1:size(Zvals,2)
        Zvals(i)=(Zbinedges(i)+Zbinedges(i+1))/2;
    end

    ZEloss=zeros(size(Zbinedges,2)-1,size(Eloss,2));

    for i=1:size(Eloss,2)
        for j=1:size(Eloss,1)
            if Eloss(j,i)>0 && Z(j,i)<maxZ && Z(j,i)>minZ
                Zind=int8(floor((Z(j,i)-Zbinedges(1))/dZ))+int8(1);
                ZEloss(Zind,i)=ZEloss(Zind,i)+Eloss(j,i);
            end
        end       
    end

    %ZEpower=ZEloss/(time(2)-time(1))*sqrt(2.699999999999991e-02);
    ZEpower=ZEloss/(time(2)-time(1));

    Epower=sum(ZEpower,1);
    
    flagRZ=contourc(RF,ZF,double(FLAG)');
    Rin=flagRZ(1,end-1);

    ZEpowerflux=ZEpower/(2*pi*Rin*dZ);

    ZEpowerfluxlevels=linspace(0,max(ZEpowerflux(:)),50)/10^6;

    Zcon=Z;
    Zcon(conflag<1)=NaN;
    

    
    fig=figure;
    
    box on
        
    newmap=plasma;
    newmap(1,:)=[1,1,1];
    
    colormap(newmap)
    
    contourf(time,Zvals,ZEpowerflux/10^6,ZEpowerfluxlevels,'linecolor','none')
    
    set(gca,'Fontsize',20)

    ylabel('$Z\,({\rm m})$','interpreter','latex','fontsize',22)
    xlabel('$t\,(\rm s)$','interpreter','latex','fontsize',22)
    
    cb=colorbar();
    cb.Label.String='RE Power Flux $[{\rm MJ/m^2}]$';
    cb.Label.Interpreter='latex';
    cb.FontSize=20;
    cb.Label.FontSize=22;
    colormap(newmap)
    
    
    print(fig,'D3D177301_REpowerflash_TEST36c','-djpeg')
    
    outputdata=0;    
    if outputdata==1
        filename='fig14.h5';
        h5create(filename,'/time',size(time))
        h5write(filename,'/time',time)

        h5create(filename,'/Zvals',size(Zvals))
        h5write(filename,'/Zvals',Zvals)    

        h5create(filename,'/ZEpowerflux',size(ZEpowerflux))
        h5write(filename,'/ZEpowerflux',ZEpowerflux)
    end
    
end

if decon2D==1
    fig=figure;
    
    time2D=ones(size(Z));
    for ii=1:size(time,2)
        time2D(:,ii)=time2D(:,ii)*time(ii);
    end
    
    Tbinedges=linspace(time(1),time(end),size(time,2)+1);
    Zbinedges=linspace(min(Z(:)),max(Z(:)),50);
    
    hist=histogram2(time2D(conflag>0),...
        Z(conflag>0),Tbinedges,Zbinedges,...
        'DisplayStyle','tile','LineStyle','none',...
        'Normalization','count');
    
    set(gca,'Fontsize',12)
    
    title('DIII-D 177301: Cumulative RE deposition along inner wall','fontsize',16)
    ylabel('$Z\,({\rm m})$','interpreter','latex','fontsize',16)
    xlabel('$t\,(\rm s)$','interpreter','latex','fontsize',16)
    
    cb=colorbar();
    cb.Label.String='Cumulative $N_{\rm RE}$';
    cb.Label.Interpreter='latex';
    cb.Label.FontSize=16;
    
    print(fig,'D3D177301_REdecon_TEST36c','-dpdf')
    
end

if deconmovie==1
    
    loops = size(time,2);
    F(loops) = struct('cdata',[],'colormap',[]);
    
    figure;
    vidfile = VideoWriter('decon_TEST36c.mp4','MPEG-4');
    open(vidfile);
       
    PHI2pi=mod(PHI,2*pi);
    
    flagRZ=contourc(RF,ZF,double(FLAG)');
    
    PHIbinedges=linspace(0,2*pi,50);
    %Zbinedges=linspace(flagRZ(2,end-49),flagRZ(2,end-1),50);
    Zbinedges=linspace(min(Z(:)),max(Z(:)),50);
    
    maxhist=0;
    for ii=1:loops  
        [hist,~,~]=histcounts2(PHI2pi(conflag(:,ii)>0,ii),...
            Z(conflag(:,ii)>0,ii),PHIbinedges,Zbinedges,'Normalization','count');        
        
        maxhisttmp=max(hist(:));
        maxhist=max(maxhist,maxhisttmp);        
    end
    
    
    for ii=1:loops                   
        hist=histogram2(PHI2pi(conflag(:,ii)>0,ii),...
            Z(conflag(:,ii)>0,ii),PHIbinedges,Zbinedges,...
            'DisplayStyle','tile','LineStyle','none',...
            'Normalization','count');
        
        format shortG
        title(strcat('t=',num2str(time(ii),'%1.4f'),'s'))
        ylabel('$Z\,({\rm m})$','interpreter','latex','fontsize',14)
        xlabel('$\phi\,(\rm rad)$','interpreter','latex','fontsize',14)
        
        cb=colorbar();
        cb.Label.String='Cumulative $N_{\rm RE}$';
        cb.Label.Interpreter='latex';
        cb.Label.FontSize=14;
        caxis([0,maxhist]);
        
        drawnow
        F(ii)=getframe(gcf);
        writeVideo(vidfile,F(ii));
    end
    
    close(vidfile)
end

if deconloc==1
   fig=figure;
   
   maxtheta=max(rad2deg(THETA(conflag(:,end)>0,end)));
   mintheta=min(rad2deg(THETA(conflag(:,end)>0,end)));
   
   
   scatter(rad2deg(PHI(conflag(:,end)>0,end)),rad2deg(THETA(conflag(:,end)>0,end)),'+')
   
   axis([0,360,0,360])
   
   xlabel('$\phi(^\circ)$','interpreter','latex')
   ylabel('$\theta_{\rm geo}(^\circ)$','interpreter','latex')
   title('$\delta B=10{\rm G}$, $E_0=10{\rm MeV}$, $\eta_0=170^\circ$','interpreter','latex')
   
   
   txt=strcat('$\%$',num2str(round(100*(maxtheta-mintheta)/360,1)),' wetted $\theta_{\rm geo}$');
   text(.7,.075,txt,'FontSize',16,'Units','normalized','interpreter','latex');
   
   set(gca,'Fontsize',18)
   box on
   
   print(fig,'JET_95135_deconloc_TEST36c','-dpdf')
end
