set(0,'defaultfigurecolor',[1 1 1])
format longE

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
time=0;

BX1=real((ReBX+ImBX*1j)*exp(1j*(omega*time+nn*phi)));
BY1=real((ReBY+ImBY*1j)*exp(1j*(omega*time+nn*phi)));
BZ1=real((ReBZ+ImBZ*1j)*exp(1j*(omega*time+nn*phi)));

EX1=real((ReEX+ImEX*1j)*exp(1j*(omega*time+nn*phi)));
EY1=real((ReEY+ImEY*1j)*exp(1j*(omega*time+nn*phi)));
EZ1=real((ReEZ+ImEZ*1j)*exp(1j*(omega*time+nn*phi)));

Rq=2.25;
Zq=0.0;

ReBX1q=interp2(ZZa,RRa,ReBX,Zq,Rq,'spline');
ReBY1q=interp2(ZZa,RRa,ReBY,Zq,Rq,'spline');
ReBZ1q=interp2(ZZa,RRa,ReBZ,Zq,Rq,'spline');

ImBX1q=interp2(ZZa,RRa,ImBX,Zq,Rq,'spline');
ImBY1q=interp2(ZZa,RRa,ImBY,Zq,Rq,'spline');
ImBZ1q=interp2(ZZa,RRa,ImBZ,Zq,Rq,'spline');

BX1q=interp2(ZZa,RRa,BX1,Zq,Rq,'spline');
BY1q=interp2(ZZa,RRa,BY1,Zq,Rq,'spline');
BZ1q=interp2(ZZa,RRa,BZ1,Zq,Rq,'spline');

%% LOAD EFIT DATA

filename='/Users/21b/Desktop/AORSA/171089/EquilibriumData/g171089.05450';

fileID = fopen(filename,'r');
formatSpec1 = '%f';
formatSpec4 = '%f %f %f %f';

init = textscan(fileID,'%s %s %s %s %d %d %d',1);

tmp=split(init{3},"#");
pulse=tmp{2};
tmp1=split(init{4},"m");
time=str2double(tmp1{1})/1e3;

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

Re=linspace(R1,R1+Rdim,NRe);
Ze=linspace(Zmid-Zdim/2,Zmid+Zdim/2,NZe);

[ZZe,RRe]=meshgrid(Ze,Re);

PSIP=reshape(B{1},[NZe,NRe]);

PSIP1=interp2(ZZe,RRe,PSIP,ZZa,RRa);

%PSIP2=permute(PSIP1,[2,1]);

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

FLAG=inpolygon(RRa,ZZa,Rlim,Zlim);

%% analysis

[dPSIdZ,dPSIdR]=gradient(PSIP1,Za(2)-Za(1),Ra(2)-Ra(1));

BR=-1./RRa.*dPSIdZ;
BZ=1./RRa.*dPSIdR;
BPHI=zeros(size(BR));
for ii=1:size(Ra,1)
    BPHI(ii,:)=R0*Bt0/Ra(ii);
end

Rq=2.25;
Zq=0.0;

PSIq=interp2(ZZa,RRa,PSIP1,Zq,Rq,'spline');
dPSIdRq=interp2(ZZa,RRa,dPSIdR,Zq,Rq,'spline');
dPSIdZq=interp2(ZZa,RRa,dPSIdZ,Zq,Rq,'spline');
BRq=interp2(ZZa,RRa,BR,Zq,Rq,'spline');
BPHIq=interp2(ZZa,RRa,BPHI,Zq,Rq,'spline');
BZq=interp2(ZZa,RRa,BZ,Zq,Rq,'spline');

%% calculating safety factor

flagPSI=interp2(ZZa,RRa,PSIP1,ZLCFS,RLCFS,'spline');

PSIPlim=min(flagPSI(:));

PSIP0=min(PSIP1(:));

PSIParray=linspace(PSIP0,PSIPlim,50);
PSINarray=(PSIParray-PSIP0)/(PSIPlim-PSIP0);

PSIT=zeros(1,50);
for jj=1:50
    BPHItmp=BPHI;
    BPHItmp(PSIP1>PSIParray(jj))=0;
    PSIT(jj)=trapz(Ra,trapz(Za,BPHItmp,1));
end

sfac=gradient(smooth(PSIT,15),2*pi*(PSIParray(2)-PSIParray(1)));

%% define toroidal fields

%!!! B0 should be positive, KORC field computations put in the (-) sign!!!%

Z0=0.;

B0=abs(Bt0);
E0=0.;

%% Save to KORC input file

sav=1;

if sav==1
    
    delete AORSA_D3D_171089.h5
    
    fields2hdf(Ra,[],Za,ReBX,ReBY,[],[],ReBZ,ImBX,ImBY,[],[],ImBZ,[],...
        ReEX,ReEY,ReEZ,ImEX,ImEY,ImEZ,[],[],[],PSIP1,...
        FLAG,[],[],[],'AORSA_D3D_171089.h5',B0,E0,RmAxis,ZmAxis,psiAxis,psiSep)
end

%% Plotting

figure;
plot(PSINarray,abs(sfac),'linewidth',3) 
hold on
plot(linspace(0,1,size(tmp3,1)),tmp3,'linewidth',3)

title('Safety factor','interpreter','latex')
xlabel('$\psi_N$','interpreter','latex')
ylabel('$q$','interpreter','latex')
legend({'Calculated','EFIT'},'location','northwest')

set(gca,'fontsize',16)
   

fig=figure;

contourf(Ra,Za,PSIP1',20,'linecolor','none')
hold on
plot(Rlim,Zlim,'k','linewidth',2)
plot(RLCFS,ZLCFS,'k--','linewidth',1)
scatter(RmAxis,ZmAxis,100,'+','k')
colorbar()
daspect([1,1,1])
title('$-\psi_p/2\pi$','interpreter','latex')
xlabel('$R({\rm m})$','interpreter','latex')
ylabel('$Z({\rm m})$','interpreter','latex')
set(gca,'fontsize',16)

fig=figure('position',[100,100,900,350]);

subplot(1,3,1)
contourf(Ra,Za,BR',20,'linecolor','none')
hold on
plot(Rlim,Zlim,'k','linewidth',2)
plot(RLCFS,ZLCFS,'k--','linewidth',1)
scatter(RmAxis,ZmAxis,100,'+','k')
colorbar()
daspect([1,1,1])
title('$B_{R,0}({\rm T})$','interpreter','latex')
xlabel('$R({\rm m})$','interpreter','latex')
ylabel('$Z({\rm m})$','interpreter','latex')
set(gca,'fontsize',16)

subplot(1,3,2)
contourf(Ra,Za,BPHI',20,'linecolor','none')
hold on
plot(Rlim,Zlim,'k','linewidth',2)
plot(RLCFS,ZLCFS,'k--','linewidth',1)
scatter(RmAxis,ZmAxis,100,'+','k')
colorbar()
daspect([1,1,1])
title('$B_{\phi,0}({\rm T})$','interpreter','latex')
xlabel('$R({\rm m})$','interpreter','latex')
%ylabel('$Z({\rm m})$','interpreter','latex')
yticklabels([])
set(gca,'fontsize',16)

subplot(1,3,3)
contourf(Ra,Za,BZ',20,'linecolor','none')
hold on
plot(Rlim,Zlim,'k','linewidth',2)
plot(RLCFS,ZLCFS,'k--','linewidth',1)
scatter(RmAxis,ZmAxis,100,'+','k')
colorbar()
daspect([1,1,1])
title('$B_{Z,0}({\rm T})$','interpreter','latex')
xlabel('$R({\rm m})$','interpreter','latex')
%ylabel('$Z({\rm m})$','interpreter','latex')
yticklabels([])
set(gca,'fontsize',16)

fig=figure('Position', [10 10 900 350]);

colormap redblue

subplot(1,3,1)

maxBX1=max(BX1(:));
minBX1=min(BX1(:));
maxabsBX1=max(maxBX1,abs(minBX1));
BX1cont=linspace(-maxabsBX1,maxabsBX1,100);

contourf(Ra,Za,BX1',BX1cont,'linecolor','none')
hold on
plot(Rlim,Zlim,'k','linewidth',2)
plot(RLCFS,ZLCFS,'k--','linewidth',1)
scatter(RmAxis,ZmAxis,100,'+','k')
caxis([-maxabsBX1,maxabsBX1])
colorbar()
daspect([1 1 1])
title('$B_{X,1}({\rm T})$','interpreter','latex')
xlabel('$R({\rm m})$','interpreter','latex')
ylabel('$Z({\rm m})$','interpreter','latex')
set(gca,'fontsize',16)

subplot(1,3,2)

maxBY1=max(BY1(:));
minBY1=min(BY1(:));
maxabsBY1=max(maxBY1,abs(minBY1));
BY1cont=linspace(-maxabsBY1,maxabsBY1,100);

contourf(Ra,Za,BY1',BY1cont,'linecolor','none')
hold on
plot(Rlim,Zlim,'k','linewidth',2)
plot(RLCFS,ZLCFS,'k--','linewidth',1)
scatter(RmAxis,ZmAxis,100,'+','k')
caxis([-maxabsBY1,maxabsBY1])
colorbar()
daspect([1 1 1])
title('$B_{Y,1}({\rm T})$','interpreter','latex')
xlabel('$R({\rm m})$','interpreter','latex')
%ylabel('$Z({\rm m})$','interpreter','latex')
yticklabels([])
set(gca,'fontsize',16)

subplot(1,3,3)

maxBZ1=max(BZ1(:));
minBZ1=min(BZ1(:));
maxabsBZ1=max(maxBZ1,abs(minBZ1));
BZ1cont=linspace(-maxabsBZ1,maxabsBZ1,100);

contourf(Ra,Za,BZ1',BZ1cont,'linecolor','none')
hold on
plot(Rlim,Zlim,'k','linewidth',2)
plot(RLCFS,ZLCFS,'k--','linewidth',1)
scatter(RmAxis,ZmAxis,100,'+','k')
caxis([-maxabsBZ1,maxabsBZ1])
colorbar()
daspect([1 1 1])
title('$B_{Z,1}({\rm T})$','interpreter','latex')
xlabel('$R({\rm m})$','interpreter','latex')
%ylabel('$Z({\rm m})$','interpreter','latex')
yticklabels([])
set(gca,'fontsize',16)
