% plot out data for ORNL data file
close all
clear all

path = '/media/l8c/FantomHD/DIII-D/Data/D3D_Eric/';

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

% **** plot density data *****
ne = nDIIV + nArIIV + 2*nArIIIV;

densfigure = figure;
plot(rhoV,ne,'k-',rhoV,nArIV,'r--',rhoV,nArIIV,'b--',rhoV,nArIIIV,'m--',rhoV,nDIV,'c.-',rhoV,nDIIV,'.-')
legend({'$n_e$','Ar','Ar$^{+}$','Ar$^{+2}$','D','D$^{+}$'},'interpreter','latex')
xlabel('$\rho$','interpreter','latex')
ylabel('[$1\times 10^{13}$/cm$^3$]','interpreter','latex')
title('Electron density vs radius')
 
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
ctr = contour(rg,zg,rhoMuse,[0.95 0.95],'g'); % contour used for separatrix
contour(rg,zg,rhoMuse,7,'b'); % contour used for separatrix
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
hold on
plot(rg(II),zg(I(II)),'g+','LineWidth',2)
hold off

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
hold on
plot(rg(II),zg(I(II)),'g+','LineWidth',2)
hold off
 
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
hold on
plot(rg(II),zg(I(II)),'g+','LineWidth',2)
hold off
 
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



