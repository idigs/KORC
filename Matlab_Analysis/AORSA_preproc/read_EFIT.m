%%read_EFIT.m

format longE
set(0,'defaultfigurecolor',[1 1 1])

%% LOAD DATA

filename='/Users/21b/Desktop/AORSA/171089/EquilibriumData/g171089.05450';

fileID = fopen(filename,'r');
formatSpec1 = '%f';
formatSpec4 = '%f %f %f %f';

init = textscan(fileID,'%s %s %s %s %d %d %d',1);

tmp=split(init{3},"#");
pulse=tmp{2};
tmp1=split(init{4},"m");
time=str2double(tmp1{1})/1e3;

NR=init{end-1};
NZ=init{end};

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

A = textscan(fileID,formatSpec1,NR);
tmp=A{1};

A = textscan(fileID,formatSpec1,NR);
tmp1=A{1};

A = textscan(fileID,formatSpec1,NR);
tmp2=A{1};

B = textscan(fileID,formatSpec1,NR*NZ); % PSIP

A = textscan(fileID,formatSpec1,NR);
tmp3=A{1};

A = textscan(fileID,formatSpec1,2);
Nlcfs=A{1}(1);
Nwall=A{1}(2);

C = textscan(fileID,formatSpec1,Nlcfs*2); % LCFS
D = textscan(fileID,formatSpec1,Nwall*2); % Limiter

fclose(fileID);

%% MANIPULATE DATA

R=linspace(R1,R1+Rdim,NR);
Z=linspace(Zmid-Zdim/2,Zmid+Zdim/2,NZ);

PSI=reshape(B{1},[NZ,NR]);

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

[dPSIdZ,dPSIdR]=gradient(PSI,Z(2)-Z(1),R(2)-R(1));
[ZZ,RR]=meshgrid(Z,R);

BR=-1./RR.*dPSIdZ;
BZ=1./RR.*dPSIdR;
BPHI=zeros(size(BR));
for ii=1:size(R,2)
    BPHI(ii,:)=R0*Bt0/R(ii);
end

Rq=2.25;
Zq=0.0;

PSIq=interp2(ZZ,RR,PSI,Zq,Rq);
dPSIdRq=interp2(ZZ,RR,dPSIdR,Zq,Rq);
dPSIdZq=interp2(ZZ,RR,dPSIdZ,Zq,Rq);
BRq=interp2(ZZ,RR,BR,Zq,Rq);
BPHIq=interp2(ZZ,RR,BPHI,Zq,Rq);
BZq=interp2(ZZ,RR,BZ,Zq,Rq);

%% calculating safety factor

flagPSI=interp2(ZZ,RR,PSI,Zlim,Rlim,'spline');

PSIPlim=min(flagPSI(:));

PSIP0=min(PSI(:));

PSIParray=linspace(PSIP0,PSIPlim,50);
PSINarray=(PSIParray-PSIP0)/(PSIPlim-PSIP0);

PSIT=zeros(1,50);
for jj=1:50
    BPHItmp=BPHI;
    BPHItmp(PSI>PSIParray(jj))=0;
    PSIT(jj)=trapz(R,trapz(Z,BPHItmp,2));
end

sfac=gradient(smooth(PSIT,15),2*pi*(PSIParray(2)-PSIParray(1)));


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

contourf(R,Z,PSI',20,'linecolor','none')
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
contourf(R,Z,BR',20,'linecolor','none')
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
contourf(R,Z,BPHI',20,'linecolor','none')
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
contourf(R,Z,BZ',20,'linecolor','none')
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

print(fig,'JET_95135_EFIT_B0','-depsc')
