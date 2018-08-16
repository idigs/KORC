% plot out data for ORNL data file
close all
clear all

path = '/Volumes/LeopoJD/DIII-D/Data/D3D_Eric/';

myfile = sprintf([path 'ORNL_RE_data_2017.mat']);
load(myfile)
 
% data loaded is: 
% tVjfitms - time vector for jfit [ms]
% psiMjfit - matrices of poloidal flux [T-m^2] from jfit on time grid tVjfitms
% rhoMjfit = matrices of normalized minor radius (r/a) from jfit
% rg - radial grid for jfit [m]
% zg - vertical grid from jfit [m]
% B0 - toroidal magnetic field on axis [T]
% Ip - toroidal plasma current [kA]
% R0 - radius of magnetic axis [m]. Can reconstruct B vector in normal fashion: BT = B0*R0/R, BZ = (dpsi/dR)/R, BR = -(dpsi/dZ)/R
% Vloop - toroidal loop voltage [V]
% rhoV - normalized radius vector for densities
% nArIV - density of neutral argon vs r/a [1e13/cm^3]
% nArIIV - density of Ar+ vs r/a [1e13/cm^3]
% nArIIIV - density of Ar2+ vs r/a [1e13/cm^3]
% nDIV - density of D vs r/a [1e13/cm^3]
% nDIIV - density of D+ vs r/a [1e13/cm^3]. Electron density ne can be obtained from charge conservation
% tplotV - time vector for IR images [ms]
% BIRrA - matrices of IR brightness [1e13 photons/cm^2/s/ster/nm], averaged over wavelength range 3 - 5 um, at times tplotV
% RIRV - radius vector for IR images at tangency plane of image [m]
% ZIRV - vertical vector for IR images at tangecy plane of image [m]
% EfitVc - kinetic energy grid for estimated RE distribution function [MeV]
% ffitV1 - RE distribution function [/cm^3/MeV] on energy grid EfitVc
% thetafitV1 - mean pitch angle [radians] of REs on energy grid EfitVc
% RvertIR, ZvertIR - coordinates of IR camera exit pupil [m]
% limdatax - outline of DIII-D walls for plotting on images [m]
% BvisA - matrix of visible brightness [1e13 photons/cm^2/s/ster/nm] averaged over wavelength range 700 - 800 nm at time tvisimg
% RvisV - radius vector for vis image at tangency plane of image [m]
% ZvisV = vertical vector for vis image at tangency plane of image [m]
% rotangle = rotation angle of vis image [degrees] 
% RtanV - radius vector for vis spectroscopy (600 nm) at tangecy plane in midplane [m]
% BforV - brightness in forward direction along midplane [1e13 photons/cm^2/s/ster/nm]
% BrevV - brightness in reverse direction along midplane [1e13 photons/cm^2/s/ster/nm]

pscaleval = 2; % scale factor for plotting of vis image
 
% Generate FLAG
[RG,ZG]=meshgrid(rg,zg);
load([path 'wall.mat']);

nr = numel(rg);
nz = numel(zg);

% fig = figure;
% 
% FLAG = zeros(nr,nz);
% for rr=1:nr
%     for zz=1:nz
%         a = atan2(zg(zz),rg(rr)-Ro);
%         if (a<0); a=a+2*pi; end;
%         r = sqrt( (rg(rr)-Ro).^2 + zg(zz).^2 );
%         
%         [~,ia] = min(abs(a-angle));
%         
%         if ( (ia>1) && (ia~=numel(angle)) )
%             d = a - angle(ia);
%             dp = a - angle(ia+1);
%             dm = a - angle(ia-1);
%         elseif (ia == numel(angle))
%             d = sign(a - angle(ia));
%             dp = sign(d);
%             dm = -sign(d);
%         elseif (ia == 1)
%             d = sign(a - angle(ia));
%             if (sign(d) ~= 0)
%                 dp = -sign(d);
%                 dm = sign(d);
%             else
%                 dp = 1;
%             end
%         end
%         
%         if ( sign(d) == sign(dp) )
%             m = (radius(ia) - radius(ia-1))/(angle(ia) - angle(ia-1));
%             y = m*(a - angle(ia-1)) + radius(ia-1);
%         else
%             m = (radius(ia+1) - radius(ia))/(angle(ia+1) - angle(ia));
%             y = m*(a - angle(ia)) + radius(ia);
%         end
%         
%         if (r > y)
%             FLAG(rr,zz) = 0;
%             figure(fig)
%             hold on;plot(rg(rr),zg(zz),'r.');hold off
%         else
%             FLAG(rr,zz) = 1;
%             figure(fig)
%             hold on;plot(rg(rr),zg(zz),'k.');hold off
%         end
%     end
% end


fig1 = figure;
ne = nDIIV + nArIIV + 2*nArIIIV;
Zeff = (nDIIV + nArIIV + 4*nArIIIV)./ne;

ep = 8.854E-12;% Electric permitivity
c = 2.9979E8; % Speed of light
me = 9.109382E-31; % Electron mass
qe = 1.602176E-19; % Electron charge

To = 2.0;
Te = -rhoV + To;
Clog = 32.2 - 1.15*log10(ne*1E19) + 2.3*log10(Te);
ECH = qe^3*(ne*1E19).*Clog/(4*pi*ep^2*me*c^2);
Ebar = (Vloop/(2*pi*R0))./ECH;
C = 2*Ebar./(Zeff + 1);

figure(fig1)

subplot(5,1,1)
plot(rhoV,ne*1E19)
xlabel('$\rho/a$','interpreter','latex')
ylabel('$n_e$ (m$^{-3}$)','interpreter','latex')

subplot(5,1,2)
plot(rhoV,Zeff)
xlabel('$\rho/a$','interpreter','latex')
ylabel('$Z_{eff}$','interpreter','latex')

subplot(5,1,3)
plot(rhoV,Te)
xlabel('$\rho/a$','interpreter','latex')
ylabel('$T_e$ (eV)','interpreter','latex')

subplot(5,1,4)
plot(rhoV,Ebar)
xlabel('$\rho/a$','interpreter','latex')
ylabel('$\bar{E}$','interpreter','latex')

subplot(5,1,5)
plot(rhoV,C)
xlabel('$\rho/a$','interpreter','latex')
ylabel('$\mathcal{C}$','interpreter','latex')

% **** plot density data *****
densfigure = figure;
subplot(2,1,1)
plot(rhoV,nArIV,'k')
hold on
plot(rhoV,nArIIV,'r')
plot(rhoV,nArIIIV,'b')
legend({'Ar','Ar$^{+}$','Ar$^{+2}$'},'interpreter','latex')
xlabel('rho')
ylabel('[1e13/cm^3]')
title('Ar densities vs radius')

subplot(2,1,2)
plot(rhoV,nDIV,'k')
hold on
plot(rhoV,nDIIV,'r')
legend({'D','D$^{+}$'},'interpreter','latex')
xlabel('rho')
ylabel('[1e13/cm^3]')
title('D densities vs radius')
 
% plot initial distribution function estimate
fEfigure = figure;
subplot(2,2,1)
loglog(EfitVc,ffitV1,'k.-') % starting distribution function
xlabel('electron energy [MeV]')
ylabel('f_E [/cm^3/MeV]')
axis([1e-3 100 1e5 1e15])
title('distribution function')
subplot(2,2,2)
f = sin(thetafitV1);
fy = rad2deg(thetafitV1);
semilogx(EfitVc,fy,'k.-');
xlabel('electron energy [MeV]')
ylabel('sin(\theta)')
axis([1e-3 100 0 max(fy)])
title('pitch angle')

subplot(2,2,4)
plot3(EfitVc,rad2deg(thetafitV1),ffitV1,'k.-')
grid minor;box on;
xlabel('electron energy [MeV]')
ylabel('\theta')
zlabel('f_E [/cm^3/MeV]')

c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass

g = 1 + 1E6*EfitVc*qe/(me*c^2);
v = c*sqrt(1 - 1./g.^2);
vpar = v.*cos(thetafitV1);
jpar = vpar.*ffitV1;

figure(fEfigure)
subplot(2,2,3)
semilogx(EfitVc,jpar,'k.-');
xlabel('electron energy [MeV]')
ylabel('$j_{\parallel}$','Interpreter','latex')


lambda = [3E-6 5E-6]; % 3 um
fRE_E = ffitV1/trapz(EfitVc,ffitV1);
eta = thetafitV1;
for ll=1:numel(lambda)
    P = zeros(size(g));
    
    for gg=1:numel(g)
        P(gg) =  fRE_E(gg)*singleParticleSpectrum(B0,lambda(ll),g(gg),eta(gg));
    end
    subplot(2,2,4)
    hold on;
    plot3(EfitVc,rad2deg(thetafitV1),P,'.-')
    hold off
    grid minor;box on;
end


% plot out IR camera images
fnorm = 0.2; % normalization constant for images
IRcamfigure = figure;
fIR1 = BIRrA(1,:,:);
fIR1 = squeeze(fIR1);
fIR1 = flipud(fIR1);
fIR1 = fliplr(fIR1); % flip left right so CP is on left side like for efits and jfits
fIR1 = fnorm*fIR1; % normalize so values go from 0-100 for plot
subplot(2,2,1)
image(RIRV,ZIRV,fIR1)
set(gca,'Ydir','Normal')
hold on
plot(limdatax(2,:),limdatax(1,:),'c')
iit1jfit = max(find(tVjfitms<=tplotV(1)));
rhoMuse = squeeze(rhoMjfit(iit1jfit,:,:));
psiMuse = squeeze(psiMjfit(iit1jfit,:,:));
ctr = contour(rg,zg,rhoMuse,[0.95 0.95],'g'); % contour used for separatrix
% contour(rg,zg,rhoMuse,7,'b'); % contour used for separatrix
contour(rg,zg,psiMuse,7,'b'); % contour used for separatrix
myString = ['IR cam, t = ' num2str(tplotV(1),'%0.5g')];
title(myString)
axis equal
colormap('jet');%colormap('hot');
xlabel('R [m]')
ylabel('Z [m]')

[V,I] = min(abs(rhoMuse));
[~,II] = min(V);
figure(IRcamfigure)
hold on
plot(rg(II),zg(I(II)),'g+','LineWidth',2)
hold off

ctr(:,1:2) = [];
Ro = rg(II);
Zo = zg(I(II));

% rspx = sqrt( (ctr(1,:)-Ro).^2 + (ctr(2,:)-Zo).^2 );
% aspx = atan2(ctr(2,:)-Zo,ctr(1,:)-Ro);
% aspx(aspx<0) = aspx(aspx<0) + 2*pi;
% aspx = flip(aspx);
% rspx = flip(rspx);
% [~,I] = min(aspx);
% aspx = [aspx(I:end) aspx(1:I-1) 2*pi];
% rspx = [rspx(I:end) rspx(1:I-1) rspx(I)];
% 
% figure(fig)
% hold on;
% plot(ctr(1,:),ctr(2,:),'g','LineWidth',2);
% plot(limdatax(2,:),limdatax(1,:),'g','LineWidth',2)
% hold off
% 
% FLUX = -psiMuse';
% 
% ZEFF = zeros(nr,nz);
% NE = zeros(nr,nz);
% TE = zeros(nr,nz);
% for rr=1:nr
%     for zz=1:nz
%         a = atan2(zg(zz)-Zo,rg(rr)-Ro);
%         if (a<0); a=a+2*pi; end;
%         r = sqrt( (rg(rr)-Ro).^2 + (zg(zz)-Zo).^2 );
%         
%         rspx_i = interp1(aspx,rspx,a,'pchip');
%         
%         if (r < rspx_i)
%             ZEFF(rr,zz) = interp1(rhoV,Zeff,r/rspx_i,'pchip');
%             NE(rr,zz) = interp1(rhoV,ne,r/rspx_i,'pchip');
%             TE(rr,zz) = interp1(rhoV,Te,r/rspx_i,'pchip');
%             figure(fig)
%             hold on;plot(rg(rr),zg(zz),'ms');hold off
%         else
%             ZEFF(rr,zz) = Zeff(end);
%             NE(rr,zz) = ne(end);
%             TE(rr,zz) = Te(end);
%             figure(fig)
%             hold on;plot(rg(rr),zg(zz),'mx');hold off
%         end        
%     end
% end
% 
% fields2hdf(rg,[],zg,[],[],[],[],[],[],FLUX,FLAG,NE,TE,ZEFF,['JFIT_D3D_157576_t' num2str(tplotV(1)) '.h5'],-B0,Vloop/(2*pi*Ro),Ro,Zo)


fIR1 = BIRrA(2,:,:);
fIR1 = squeeze(fIR1);
fIR1 = flipud(fIR1);
fIR1 = fliplr(fIR1); % flip left right so CP is on left side like for efits and jfits
fIR1 = fnorm*fIR1; % normalize so values go from 0-100 for plot
figure(IRcamfigure)
subplot(2,2,2)
image(RIRV,ZIRV,fIR1)
set(gca,'Ydir','Normal')
hold on
plot(limdatax(2,:),limdatax(1,:),'c')
iit1jfit = max(find(tVjfitms<=tplotV(2)));
rhoMuse = squeeze(rhoMjfit(iit1jfit,:,:));
psiMuse = squeeze(psiMjfit(iit1jfit,:,:));
ctr = contour(rg,zg,rhoMuse,[0.95 0.95],'g'); % contour used for separatrix
contour(rg,zg,rhoMuse,7,'b'); % contour used for separatrix
myString = ['IR cam, t = ' num2str(tplotV(2),'%0.5g')];
title(myString)
axis equal
colormap('jet')%colormap('hot')
xlabel('R [m]')
ylabel('Z [m]')
 
[V,I] = min(abs(rhoMuse));
[~,II] = min(V);
figure(IRcamfigure)
hold on
plot(rg(II),zg(I(II)),'g+','LineWidth',2)
hold off

ctr(:,1:2) = [];
Ro = rg(II);
Zo = zg(I(II));

% rspx = sqrt( (ctr(1,:)-Ro).^2 + (ctr(2,:)-Zo).^2 );
% aspx = atan2(ctr(2,:)-Zo,ctr(1,:)-Ro);
% aspx(aspx<0) = aspx(aspx<0) + 2*pi;
% aspx = flip(aspx);
% rspx = flip(rspx);
% [~,I] = min(aspx);
% aspx = [aspx(I:end) aspx(1:I-1) 2*pi];
% rspx = [rspx(I:end) rspx(1:I-1) rspx(I)];
% 
% fig2 = figure;
% hold on;
% plot(ctr(1,:),ctr(2,:),'g','LineWidth',2);
% plot(limdatax(2,:),limdatax(1,:),'g','LineWidth',2)
% hold off
% 
% FLUX = -psiMuse';
% 
% ZEFF = zeros(nr,nz);
% NE = zeros(nr,nz);
% TE = zeros(nr,nz);
% for rr=1:nr
%     for zz=1:nz
%         a = atan2(zg(zz)-Zo,rg(rr)-Ro);
%         if (a<0); a=a+2*pi; end
%         r = sqrt( (rg(rr)-Ro).^2 + (zg(zz)-Zo).^2 );
%         
%         rspx_i = interp1(aspx,rspx,a,'pchip');
%         
%         if (r < rspx_i)
%             ZEFF(rr,zz) = interp1(rhoV,Zeff,r/rspx_i,'pchip');
%             NE(rr,zz) = interp1(rhoV,ne,r/rspx_i,'pchip');
%             TE(rr,zz) = interp1(rhoV,Te,r/rspx_i,'pchip');
%             figure(fig2)
%             hold on;plot(rg(rr),zg(zz),'ms');hold off
%         else
%             ZEFF(rr,zz) = Zeff(end);
%             NE(rr,zz) = ne(end);
%             TE(rr,zz) = Te(end);
%             figure(fig2)
%             hold on;plot(rg(rr),zg(zz),'mx');hold off
%         end        
%     end
% end
% 
% fields2hdf(rg,[],zg,[],[],[],[],[],[],FLUX,FLAG,NE,TE,ZEFF,['JFIT_D3D_157576_t' num2str(tplotV(2)) '.h5'],-B0,Vloop/(2*pi*Ro),Ro,Zo)

fIR1 = BIRrA(3,:,:);
fIR1 = squeeze(fIR1);
fIR1 = flipud(fIR1);
fIR1 = fliplr(fIR1); % flip left right so CP is on left side like for efits and jfits
fIR1 = fnorm*fIR1; % normalize so values go from 0-100 for plot
figure(IRcamfigure)
subplot(2,2,3)
image(RIRV,ZIRV,fIR1)
set(gca,'Ydir','Normal')
hold on
plot(limdatax(2,:),limdatax(1,:),'c')
iit1jfit = max(find(tVjfitms<=tplotV(3)));
rhoMuse = squeeze(rhoMjfit(iit1jfit,:,:));
psiMuse = squeeze(psiMjfit(iit1jfit,:,:));
ctr = contour(rg,zg,rhoMuse,[0.95 0.95],'g'); % contour used for separatrix
contour(rg,zg,rhoMuse,7,'b'); % contour used for separatrix
myString = ['IR cam, t = ' num2str(tplotV(3),'%0.5g')];
title(myString)
axis equal
colormap('jet')%colormap('hot')
xlabel('R [m]')
ylabel('Z [m]')
 
[V,I] = min(abs(rhoMuse));
[~,II] = min(V);
figure(IRcamfigure)
hold on
plot(rg(II),zg(I(II)),'g+','LineWidth',2)
hold off

ctr(:,1:2) = [];
Ro = rg(II);
Zo = zg(I(II));

% rspx = sqrt( (ctr(1,:)-Ro).^2 + (ctr(2,:)-Zo).^2 );
% aspx = atan2(ctr(2,:)-Zo,ctr(1,:)-Ro);
% aspx(aspx<0) = aspx(aspx<0) + 2*pi;
% aspx = flip(aspx);
% rspx = flip(rspx);
% [~,I] = min(aspx);
% aspx = [aspx(I:end) aspx(1:I-1) 2*pi];
% rspx = [rspx(I:end) rspx(1:I-1) rspx(I)];
% 
% fig3 = figure;
% hold on;
% plot(ctr(1,:),ctr(2,:),'g','LineWidth',2);
% plot(limdatax(2,:),limdatax(1,:),'g','LineWidth',2)
% hold off
% 
% FLUX = -psiMuse';
% 
% ZEFF = zeros(nr,nz);
% NE = zeros(nr,nz);
% TE = zeros(nr,nz);
% for rr=1:nr
%     for zz=1:nz
%         a = atan2(zg(zz)-Zo,rg(rr)-Ro);
%         if (a<0); a=a+2*pi; end;
%         r = sqrt( (rg(rr)-Ro).^2 + (zg(zz)-Zo).^2 );
%         
%         rspx_i = interp1(aspx,rspx,a,'pchip');
%         
%         if (r < rspx_i)
%             ZEFF(rr,zz) = interp1(rhoV,Zeff,r/rspx_i,'pchip');
%             NE(rr,zz) = interp1(rhoV,ne,r/rspx_i,'pchip');
%             TE(rr,zz) = interp1(rhoV,Te,r/rspx_i,'pchip');
%             figure(fig3)
%             hold on;plot(rg(rr),zg(zz),'ms');hold off
%         else
%             ZEFF(rr,zz) = Zeff(end);
%             NE(rr,zz) = ne(end);
%             TE(rr,zz) = Te(end);
%             figure(fig3)
%             hold on;plot(rg(rr),zg(zz),'mx');hold off
%         end        
%     end
% end
% 
% fields2hdf(rg,[],zg,[],[],[],[],[],[],FLUX,FLAG,NE,TE,ZEFF,['JFIT_D3D_157576_t' num2str(tplotV(3)) '.h5'],-B0,Vloop/(2*pi*Ro),Ro,Zo)
 
fIR1 = BIRrA(4,:,:);
fIR1 = squeeze(fIR1);
fIR1 = flipud(fIR1);
fIR1 = fliplr(fIR1); % flip left right so CP is on left side like for efits and jfits
BIRV = fIR1(410,:); % horizontal slice across downward shifted image
fIR1 = fnorm*fIR1; % normalize so values go from 0-100 for plot
figure(IRcamfigure)
subplot(2,2,4)
image(RIRV,ZIRV,fIR1)
set(gca,'Ydir','Normal')
hold on
plot(limdatax(2,:),limdatax(1,:),'c')
iit1jfit = max(find(tVjfitms<=tplotV(3)));
rhoMuse = squeeze(rhoMjfit(iit1jfit,:,:));
ctr = contour(rg,zg,rhoMuse,[0.95 0.95],'g'); % contour used for separatrix
contour(rg,zg,rhoMuse,7,'b'); % contour used for separatrix
myString = ['IR cam, t = ' num2str(tplotV(4),'%0.5g')];
title(myString)
axis equal
colormap('jet')%colormap('hot')
xlabel('R [m]')
ylabel('Z [m]')

[V,I] = min(abs(rhoMuse));
[~,II] = min(V);
figure(IRcamfigure)
hold on
plot(rg(II),zg(I(II)),'g+','LineWidth',2)
hold off

ctr(:,1:2) = [];
Ro = rg(II);
Zo = zg(I(II));

% rspx = sqrt( (ctr(1,:)-Ro).^2 + (ctr(2,:)-Zo).^2 );
% aspx = atan2(ctr(2,:)-Zo,ctr(1,:)-Ro);
% aspx(aspx<0) = aspx(aspx<0) + 2*pi;
% aspx = flip(aspx);
% rspx = flip(rspx);
% [~,I] = min(aspx);
% aspx = [aspx(I:end) aspx(1:I-1) 2*pi];
% rspx = [rspx(I:end) rspx(1:I-1) rspx(I)];
% 
% fig4 = figure;
% hold on;
% plot(ctr(1,:),ctr(2,:),'g','LineWidth',2);
% plot(limdatax(2,:),limdatax(1,:),'g','LineWidth',2)
% hold off
% 
% FLUX = -psiMuse';
% 
% ZEFF = zeros(nr,nz);
% NE = zeros(nr,nz);
% TE = zeros(nr,nz);
% for rr=1:nr
%     for zz=1:nz
%         a = atan2(zg(zz)-Zo,rg(rr)-Ro);
%         if (a<0); a=a+2*pi; end
%         r = sqrt( (rg(rr)-Ro).^2 + (zg(zz)-Zo).^2 );
%         
%         rspx_i = interp1(aspx,rspx,a,'pchip');
%         
%         if (r < rspx_i)
%             ZEFF(rr,zz) = interp1(rhoV,Zeff,r/rspx_i,'pchip');
%             NE(rr,zz) = interp1(rhoV,ne,r/rspx_i,'pchip');
%             TE(rr,zz) = interp1(rhoV,Te,r/rspx_i,'pchip');
%             figure(fig4)
%             hold on;plot(rg(rr),zg(zz),'ms');hold off
%         else
%             ZEFF(rr,zz) = Zeff(end);
%             NE(rr,zz) = ne(end);
%             TE(rr,zz) = Te(end);
%             figure(fig4)
%             hold on;plot(rg(rr),zg(zz),'mx');hold off
%         end        
%     end
% end
% 
% fields2hdf(rg,[],zg,[],[],[],[],[],[],FLUX,FLAG,NE,TE,ZEFF,['JFIT_D3D_157576_t' num2str(tplotV(4)) '.h5'],-B0,Vloop/(2*pi*Ro),Ro,Zo)
 
% *** plot VIS synchrotron image ****
myA = flipud(BvisA);
pAuse = pscaleval*myA; % normalize matrix so values go from 0-100 for plot 
% pAuse = imrotate(pAuse,rotangle);
BvisV = 0*pAuse(200,:); 
for iiZ = 1:10
  BvisV = BvisV + pAuse(190 + iiZ,:);
end
BvisV = BvisV/(10*pscaleval);% slice across middle of vis image [1e13 photons/cm^2/s/ster]
visfigure = figure;
image(RvisV,ZvisV,pAuse)  % show camera image
hold on
iit1jfit = max(find(tVjfitms<=tplotV(1)));
rhoMuse = squeeze(rhoMjfit(iit1jfit,:,:));
psiMuse = squeeze(psiMjfit(iit1jfit,:,:));
contour(rg,zg,rhoMuse,[0.95 0.95],'g') % contour used for separatrix for centered synch image
myString = ['VIS cam, t = ' num2str(tvisimg,'%0.5g')];
title(myString)
set(gca,'Ydir','Normal')
axis equal
xlabel('R [m]')
ylabel('Z [m]')
 
% plot CER forward and reverse synchrotron data (from t33 - t40 and t1 - t3 at 610 nm)
% tangency radius define as 1.7 m for getting tangency plane. Compare with IR and VIS images horizontal slice
cerfigure = figure;
plot(RtanV,BcforV,'rd')
hold on
plot(RtanV,10*BcrevV,'ks')
% plot(RvisV,BvisV,'g')
plot(RIRV,BIRV/10,'b')
xlabel('R [m]')
ylabel('Bright [1e13 photons/cm^2/s/ster]')
legend('forward [610 nm]','backward [610 nm] (x10)','forward [750 nm]','forward [4 um] (/10)')
title('Midplane synchrotron brightness')
xlim([1 2.5])



