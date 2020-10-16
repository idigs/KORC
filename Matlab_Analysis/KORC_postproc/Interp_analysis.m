format longE

%% Load KORC output file

load('ST.mat')

%ST=ST;

c=ST.params.scales.v;
me=ST.params.scales.m;
qe=ST.params.scales.q;

%% Analyze field data

%PSIP=ST.params.fields.psi_p2D;
BR=ST.params.fields.BR;
BPHI=ST.params.fields.BPHI;
BZ=ST.params.fields.BZ;
ER=ST.params.fields.ER;
EPHI=ST.params.fields.EPHI;
EZ=ST.params.fields.EZ;
if isfield(ST.params.fields,'Flag3D')
    FLAG=ST.params.fields.Flag3D;
else isfield(ST.params.fields,'Flag2D');
    FLAG=ST.params.fields.Flag2D;
end
    
RF=ST.params.fields.R;
ZF=ST.params.fields.Z;

if isfield(ST.params.fields,'Flag3D')
    PHIF=ST.params.fields.PHI;
    [PP,RR,ZZ]=meshgrid(PHIF,RF,ZF);
    
    
    BX=BR.*cos(PP)-BPHI.*sin(PP);
    BY=BR.*sin(PP)+BPHI.*cos(PP);
end


calcpsip=0;
if calcpsip==1

    [ZZ,RR]=meshgrid(ZF,RF);

    [dBRdZ,dBRdR]=gradient(BR,ZF(2)-ZF(1),RF(2)-RF(1));
    [dBZdZ,dBZdR]=gradient(BZ,ZF(2)-ZF(1),RF(2)-RF(1));

    divb=dBRdR+dBZdZ;

    B=sqrt(BR.^2+BPHI.^2+BZ.^2);

    NR=size(RF,1);
    NZ=size(ZF,1);

    B0=ST.params.fields.Bo;
    R0=ST.params.fields.Ro;

    psi=zeros(size(BR));

    for j=2:size(psi,2)
       psi(1,j)=psi(1,j-1)-(ZF(2)-ZF(1))*BR(1,j); 
    end
    for i=2:size(psi,1)
       psi(i,:)=psi(i-1,:)+(RF(2)-RF(1))*BZ(i,:); 
    end

    [dpsidZ,dpsidR]=gradient(psi,ZF(2)-ZF(1),RF(2)-RF(1));


    BR1=-1./RR.*dpsidZ;
    BZ1=1./RR.*dpsidR;
end


%% Analyze particle data

%R=squeeze(ST.data.sp1.Y(:,1,:));
%Z=squeeze(ST.data.sp1.Y(:,3,:));

X=squeeze(ST.data.sp1.X(:,1,:));
Y=squeeze(ST.data.sp1.X(:,2,:));
Z=squeeze(ST.data.sp1.X(:,3,:));

[PHI,R,Z]=cart2pol(X,Y,Z);

%PSIp=ST.data.sp1.PSIp;

BXp=squeeze(ST.data.sp1.B(:,1,:));
BYp=squeeze(ST.data.sp1.B(:,2,:));
BZp=squeeze(ST.data.sp1.B(:,3,:));

BRp=BXp.*cos(PHI)+BYp.*sin(PHI);
BPHIp=-BXp.*sin(PHI)+BYp.*cos(PHI);

Bp=sqrt(BRp.^2+BPHIp.^2+BZp.^2);

EXp=squeeze(ST.data.sp1.E(:,1,:));
EYp=squeeze(ST.data.sp1.E(:,2,:));
EZp=squeeze(ST.data.sp1.E(:,3,:));

ERp=EXp.*cos(PHI)+EYp.*sin(PHI);
EPHIp=-EXp.*sin(PHI)+EYp.*cos(PHI);

VXp=squeeze(ST.data.sp1.V(:,1,:));
VYp=squeeze(ST.data.sp1.V(:,2,:));
VZp=squeeze(ST.data.sp1.V(:,3,:));
gam=ST.data.sp1.g;
eta=ST.data.sp1.eta;
Vp=sqrt(VXp.^2+VYp.^2+VZp.^2);

gyrow=qe*Bp./(gam'*me);
gyroiw=1./gyrow;
larmorp=Vp.*sin(deg2rad(eta')).*gyroiw;

%% Plotting data

close all

Bfield=1;
Bfieldpart=0;
psip=0;
dpsip=0;
psippart=0;
Efield=0;
Ephipart=0;


if Ephipart==1
   figure;   
   s=scatter(R(:,tsliceP),Z(:,tsliceP),[],EPHIp(:,tsliceP),'filled');
   s.MarkerEdgeColor = 'k';
   hold on
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('EPHI_{part}')
   xlabel('R(m)')
   ylabel('Z(m)')
end

if Bfieldpart==1
   figure;

   subplot(1,3,1)
   scatter(R(:,tsliceP),Z(:,tsliceP),[],BRp(:,tsliceP))
   hold on
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('BR_{part}')
   xlabel('R(m)')
   ylabel('Z(m)')

   subplot(1,3,2)
   scatter(R(:,tsliceP),Z(:,tsliceP),[],BPHIp(:,tsliceP))
   hold on
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('BPHI_{part}')
   xlabel('R(m)')
   ylabel('Z(m)')
   
   
   subplot(1,3,3)
   scatter(R(:,tsliceP),Z(:,tsliceP),[],BZp(:,tsliceP))
   hold on
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('BZ_{part}')
   xlabel('R(m)')
   ylabel('Z(m)')
   
end

if Efield==1
    fig=figure('Position', [10 10 1000 350]); 

    subplot(1,3,1)
    hold on
    if ndims(ER)==2
        contourf(RF,ZF,squeeze(ER(:,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,:))',1,'k')
    else
        contourf(RF,ZF,squeeze(ER(:,1,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,1,:))',1,'k') 
    end
    s=scatter(R(:),Z(:),[],ERp(:),'filled');
    s.MarkerEdgeColor = 'k';    
    daspect([1,1,1])
    colormap redblue
    colorbar()
    maxB=max(ER(:));
    minB=min(ER(:));
    maxB=max(maxB,abs(minB));
    caxis([-maxB,maxB])  
    title('ER')
    xlabel('R(m)')
    ylabel('Z(m)')
    hold off    
    
    subplot(1,3,2)
    hold on
    if ndims(ER)==2
        contourf(RF,ZF,squeeze(EPHI(:,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,:))',1,'k')
    else
        contourf(RF,ZF,squeeze(EPHI(:,1,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,1,:))',1,'k') 
    end
    s=scatter(R(:),Z(:),[],EPHIp(:),'filled');
    s.MarkerEdgeColor = 'k';    
    daspect([1,1,1])
    colormap redblue
    colorbar()
    maxB=max(EPHI(:));
    minB=min(EPHI(:));
    maxB=max(maxB,abs(minB));
    caxis([-maxB,maxB])  
    title('EPHI')
    xlabel('R(m)')
    ylabel('Z(m)')
    hold off
    
    subplot(1,3,3)
    hold on
    if ndims(ER)==2
        contourf(RF,ZF,squeeze(EZ(:,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,:))',1,'k')
    else
        contourf(RF,ZF,squeeze(EZ(:,1,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,1,:))',1,'k') 
    end
    s=scatter(R(:),Z(:),[],EZp(:),'filled');
    s.MarkerEdgeColor = 'k';    
    daspect([1,1,1])
    colormap redblue
    colorbar()
    maxB=max(EZ(:));
    minB=min(EZ(:));
    maxB=max(maxB,abs(minB));
    caxis([-maxB,maxB])   
    title('EZ')
    xlabel('R(m)')
    ylabel('Z(m)')
    hold off    
end

if Bfield==1
    fig=figure('Position', [10 10 1000 350]);

    subplot(1,3,1)
    hold on
    if ndims(BR)==2
        contourf(RF,ZF,squeeze(BR(:,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,:))',1,'k')
    else
        contourf(RF,ZF,squeeze(BR(:,1,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,1,:))',1,'k')
    end
    s=scatter(R(:,:),Z(:,:),[],BRp(:,:),'filled');
    s.MarkerEdgeColor = 'k';
    daspect([1,1,1])
    colormap plasma
    colorbar()
    maxBR=max(BR(:));
    maxBR=max(maxBR,max(BRp(:)));
    minBR=min(BR(:));
    minBR=min(minBR,min(BRp(:)));
    caxis([minBR,maxBR])   
    title('BR')
    xlabel('R(m)')
    ylabel('Z(m)')
    hold off
   
    subplot(1,3,2)
    hold on
    if ndims(BPHI)==2
        contourf(RF,ZF,squeeze(BPHI(:,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,:))',1,'k')
    else
        contourf(RF,ZF,squeeze(BPHI(:,1,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,1,:))',1,'k')
    end
    s=scatter(R(:,:),Z(:,:),[],BPHIp(:,:),'filled');
    s.MarkerEdgeColor = 'k';
    daspect([1,1,1])
    colormap plasma
    colorbar()
    maxBPHI=max(BPHI(:));
    maxBPHI=max(maxBPHI,max(BPHIp(:)));
    minBPHI=min(BPHI(:));
    minBPHI=min(minBPHI,min(BPHIp(:)));
    caxis([minBPHI,maxBPHI])
    title('BPHI')
    xlabel('R(m)')
    ylabel('Z(m)')
    hold off

    subplot(1,3,3)
    hold on
    if ndims(BZ)==2
        contourf(RF,ZF,squeeze(BZ(:,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,:))',1,'k')
    else
        contourf(RF,ZF,squeeze(BZ(:,1,:))',200,'LineStyle','none')
        contour(RF,ZF,squeeze(FLAG(:,1,:))',1,'k')
    end
    s=scatter(R(:,:),Z(:,:),[],BZp(:,:),'filled');
    s.MarkerEdgeColor = 'k';
    daspect([1,1,1])
    colormap plasma
    colorbar()
    maxBZ=max(BZ(:));
    maxBZ=max(maxBZ,max(BZp(:)));
    minBZ=min(BZ(:));
    minBZ=min(minBZ,min(BZp(:)));
    caxis([minBZ,maxBZ])
    title('BZ')
    xlabel('R(m)')
    ylabel('Z(m)')
    hold off
end
   
if dpsip==1
   figure;
   
   subplot(1,3,1)
   contourf(RF,ZF,squeeze(dPSIPdR(:,tsliceF,:))',200,'LineStyle','none')
   hold on
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('dPSIPdR')
   xlabel('R(m)')
   ylabel('Z(m)')
   hold off
   
   subplot(1,3,2)
   contourf(RF,ZF,squeeze(dPSIPdT(:,tsliceF,:))',200,'LineStyle','none')
   hold on
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('dPSIPdT')
   xlabel('R(m)')
   ylabel('Z(m)')
   hold off
   
   subplot(1,3,3)
   contourf(RF,ZF,squeeze(dPSIPdZ(:,tsliceF,:))',200,'LineStyle','none')
   hold on
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('dPSIPdZ')
   xlabel('R(m)')
   ylabel('Z(m)')
   hold off
end

if psippart==1
   figure;
   scatter(R(:,tsliceP),Z(:,tsliceP),[],PSIp(:,tsliceP))
   hold on
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('PSIp_{part}')
   xlabel('R(m)')
   ylabel('Z(m)')
   
end

if psip==1
   figure;
   contourf(RF,ZF,squeeze(PSIP(:,tsliceF,:))',50)
   hold on
   scatter(R(:,tsliceP),Z(:,tsliceP),[],PSIp(:,tsliceP))
   contour(RF,ZF,squeeze(FLAG(:,tsliceF,:))',1,'k')
   daspect([1,1,1])
   colormap jet
   colorbar()
   title('PSIP')
   xlabel('R(m)')
   ylabel('Z(m)')
   
end