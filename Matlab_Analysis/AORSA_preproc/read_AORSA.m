%%read_AORSA.m
% presently, plotting uses information from read_EFIT.m

format longE
set(0,'defaultfigurecolor',[1 1 1])


%% LOAD DATA

filename='/Users/21b/Desktop/AORSA/171089/AORSA_Data/bharvey-200MHz';

fileID = fopen(filename,'r');
formatSpec1 = '%u';
formatSpec2 = '%f';
formatSpec3 = '%f %f';

init = textscan(fileID,formatSpec1,2);
NR=init{1}(1);
NZ=init{1}(2);


A = textscan(fileID,formatSpec2,NR);
B = textscan(fileID,formatSpec2,NZ);
C = textscan(fileID,formatSpec3,NR*NZ);
D = textscan(fileID,formatSpec3,NR*NZ);
E = textscan(fileID,formatSpec3,NR*NZ);
F = textscan(fileID,formatSpec3,NR*NZ);
G = textscan(fileID,formatSpec3,NR*NZ);
H = textscan(fileID,formatSpec3,NR*NZ);
I = textscan(fileID,formatSpec3,NR*NZ);
J = textscan(fileID,formatSpec3,NR*NZ);
K = textscan(fileID,formatSpec3,NR*NZ);
L = textscan(fileID,formatSpec2,1);
M = textscan(fileID,formatSpec2,NR*NZ);

fclose(fileID);

%% MANIPULATE DATA

R=A{1};
Z=B{1};

[ZZ,RR]=meshgrid(Z,R);

ReEX=reshape(F{1},[NZ,NR]);
ImEX=reshape(F{2},[NZ,NR]);
ReEY=reshape(G{1},[NZ,NR]);
ImEY=reshape(G{2},[NZ,NR]);
ReEZ=reshape(H{1},[NZ,NR]);
ImEZ=reshape(H{2},[NZ,NR]);

ReBX=reshape(I{1},[NZ,NR]);
ImBX=reshape(I{2},[NZ,NR]);
ReBY=reshape(J{1},[NZ,NR]);
ImBY=reshape(J{2},[NZ,NR]);
ReBZ=reshape(K{1},[NZ,NR]);
ImBZ=reshape(K{2},[NZ,NR]);

rho=reshape(M{1},[NZ,NR]);

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

ReBX1q=interp2(ZZ,RR,ReBX,Zq,Rq,'spline');
ReBY1q=interp2(ZZ,RR,ReBY,Zq,Rq,'spline');
ReBZ1q=interp2(ZZ,RR,ReBZ,Zq,Rq,'spline');

ImBX1q=interp2(ZZ,RR,ImBX,Zq,Rq,'spline');
ImBY1q=interp2(ZZ,RR,ImBY,Zq,Rq,'spline');
ImBZ1q=interp2(ZZ,RR,ImBZ,Zq,Rq,'spline');

BX1q=interp2(ZZ,RR,BX1,Zq,Rq,'spline');
BY1q=interp2(ZZ,RR,BY1,Zq,Rq,'spline');
BZ1q=interp2(ZZ,RR,BZ1,Zq,Rq,'spline');


%% PLOT DATA

fig=figure('Position', [10 10 900 350]);

colormap redblue

subplot(1,3,1)

maxBX1=max(BX1(:));
minBX1=min(BX1(:));
maxabsBX1=max(maxBX1,abs(minBX1));
BX1cont=linspace(-maxabsBX1,maxabsBX1,100);

contourf(R,Z,BX1',BX1cont,'linecolor','none')
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

contourf(R,Z,BY1',BY1cont,'linecolor','none')
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

contourf(R,Z,BZ1',BZ1cont,'linecolor','none')
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

print(fig,'D3D_171089_AORSA_B1','-depsc')