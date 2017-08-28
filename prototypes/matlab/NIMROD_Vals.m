% Script to read Val's magnetic fields

close all
clear all

[r,z,varr]=readcon(1100,1); % load data

[R,Z,Br]=onegrid(r,z,varr,14); % Obtain magnetic field components

[~,~,Bz]=onegrid(r,z,varr,15); % Obtain magnetic field components

[~,~,Bphi]=onegrid(r,z,varr,16); % Obtain magnetic field components

[~,~,n]=onegrid(r,z,varr,36); % Obtain magnetic field components

figure;
subplot(1,3,1)
surf(R,Z,squeeze(Br(:,:,1)),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90]);axis equal

subplot(1,3,2)
surf(R,Z,squeeze(Bphi(:,:,1)),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90]);axis equal

subplot(1,3,3)
surf(R,Z,squeeze(Bz(:,:,1)),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90]);axis equal


[X3,Y3,Z3,Br3]=make_field_real(R,Z,Br,32);

[~,~,~,Bz3]=make_field_real(R,Z,Bz,32);

[~,~,~,Bphi3]=make_field_real(R,Z,Bphi,32);

[~,~,~,n3]=make_field_real(R,Z,n,32);

B = sqrt(Br3.^2 + Bz3.^2 + Bphi3.^2);

% figure;
% ii = 1;
% subplot(1,3,1)
% surf(R,Z,B(:,:,ii),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90]);axis equal
% 
% ii = 10;
% subplot(1,3,2)
% surf(R,Z,B(:,:,ii),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90]);axis equal
% 
% ii = 20;
% subplot(1,3,3)
% surf(R,Z,B(:,:,ii),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90]);axis equal

Rmin = min(min(R));
Rmax = max(max(R));
Zmin = min(min(Z));
Zmax = max(max(Z));

% N = 10;
% hold on;
% plot(linspace(Rmin,Rmax,N),Zmax*ones(1,N),'k')
% plot(linspace(Rmin,Rmax,N),Zmin*ones(1,N),'k')
% plot(Rmax*ones(1,N),linspace(Zmin,Zmax,N),'k')
% plot(Rmin*ones(1,N),linspace(Zmin,Zmax,N),'k')
% hold off
%%
Bslice = squeeze(Br3(:,:,1));
BI = scatteredInterpolant(reshape(R,[numel(R) 1]),reshape(Z,[numel(Z) 1]),reshape(Bslice,[numel(Bslice) 1]));

NR=150;
NZ=250;
A = zeros(NR,NZ);
rAxis=linspace(Rmin,Rmax,NR);
zAxis=linspace(Zmin,Zmax,NZ);
for rr=1:NR
    for zz=1:NZ
        A(rr,zz) = BI(rAxis(rr),zAxis(zz));
    end
end

figure;surf(rAxis,zAxis,A','linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90]);axis equal
hold on;plot3(R(end,:),Z(end,:),max(max(A))*ones(1,numel(R(end,:))),'k');hold off