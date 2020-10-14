% Script to generate the electromagnetic fields of Val's simulations of
% thermal quench.
% t = [     25         0.0041
%           50	     0.0196
%          100	     0.0696
%          300	     0.2696
%          500	     0.4696
%          700	     0.6696
%          800	     0.7051
%          900	     0.7345
%         1000	     0.7674
%         1050	     0.7799
%         1100	     0.8008
%         1200	     0.8785
%         1300	     0.9645
%         1350	     1.0145
%         1500	     1.1645
%         1700	     1.3645
%         1900	     1.5645];% Diverted plasma simulation

function S = getNIMRODFields(id,NR,NPHI,NZ)
% [BR,BPHI,BZ,R,Z,Ro,Zo,Bo]
S = struct;
F = struct;
close all

[r,z,varr]=readcon_NIMUW(id,0); % load data

[R_NIMROD,Z_NIMROD,Br]=onegrid(r,z,varr,14); % Obtain magnetic field components
[~,~,Bz]=onegrid(r,z,varr,15); % Obtain magnetic field components
[~,~,Bphi]=onegrid(r,z,varr,16); % Obtain magnetic field components

B = sqrt(Br.^2 + Bphi.^2 + Bz.^2);

[~,~,Jr]=onegrid(r,z,varr,20); % Obtain current density components
[~,~,Jz]=onegrid(r,z,varr,21); % Obtain current density components
[~,~,Jphi]=onegrid(r,z,varr,22); % Obtain current density components

[~,~,Vr]=onegrid(r,z,varr,26); % Obtain velocity components
[~,~,Vz]=onegrid(r,z,varr,27); % Obtain velocity components
[~,~,Vphi]=onegrid(r,z,varr,28); % Obtain velocity components

[~,~,diffprof]=onegrid(r,z,varr,13); % Obtain diffusion profile

%[~,~,n]=onegrid(r,z,varr,36); % Obtain magnetic field components

%[~,~,n1]=onegrid(r,z,varr,57);

[S.Ro,S.Zo,S.Bo,S.bi] = getMagneticAxisParams(id);% bound(R,Z)

% [X3D,Y3D,Z3D,F.Br3D]=make_field_real(R_NIMROD,Z_NIMROD,Br,NPHI); % Br3D(Z,R,PHI)
[~,~,~,F.Br3D]=make_field_real(R_NIMROD,Z_NIMROD,Br,NPHI); % Br3D(Z,R,PHI)
% X3D(:,:,1) = [];Y3D(:,:,1) = [];Z3D(:,:,1) = [];
% X3D = flip(X3D,3);Y3D = flip(Y3D,3);Z3D = flip(Z3D,3);
F.Br3D(:,:,1) = [];
% F.Br3D = flip(F.Br3D,3);


[~,~,~,F.Bz3D]=make_field_real(R_NIMROD,Z_NIMROD,Bz,NPHI);
F.Bz3D(:,:,1) = [];
% F.Bz3D = flip(F.Bz3D,3);

[~,~,~,F.Bphi3D]=make_field_real(R_NIMROD,Z_NIMROD,Bphi,NPHI);
F.Bphi3D = -F.Bphi3D;
F.Bphi3D(:,:,1) = [];
% F.Bphi3D = flip(F.Bphi3D,3);

% [~,~,~,n3D]=make_field_real(R_NIMROD,Z_NIMROD,n,NPHI);
% n3D(:,:,1) = [];
% % n3D = flip(n3D,3);
% 
% [~,~,~,n13D]=make_field_real(R_NIMROD,Z_NIMROD,n,NPHI);
% n13D(:,:,1) = [];
% % n13D = flip(n13D,3);

[~,~,~,F.Jr3D]=make_field_real(R_NIMROD,Z_NIMROD,Jr,NPHI);
F.Jr3D(:,:,1) = [];
[~,~,~,F.Jz3D]=make_field_real(R_NIMROD,Z_NIMROD,Jz,NPHI);
F.Jz3D(:,:,1) = [];
[~,~,~,F.Jphi3D]=make_field_real(R_NIMROD,Z_NIMROD,Jphi,NPHI);
F.Jphi3D = -F.Jphi3D;
F.Jphi3D(:,:,1) = [];

[~,~,~,F.Vr3D]=make_field_real(R_NIMROD,Z_NIMROD,Vr,NPHI);
F.Vr3D(:,:,1) = [];
[~,~,~,F.Vz3D]=make_field_real(R_NIMROD,Z_NIMROD,Vz,NPHI);
F.Vz3D(:,:,1) = [];
[~,~,~,F.Vphi3D]=make_field_real(R_NIMROD,Z_NIMROD,Vphi,NPHI);
F.Vphi3D = -F.Vphi3D;
F.Vphi3D(:,:,1) = [];

[~,~,~,F.diff3D]=make_field_real(R_NIMROD,Z_NIMROD,diffprof,NPHI);
F.diff3D(:,:,1) = [];

% Fields to be used in HDF5 file
S.BR = zeros(NR,NPHI,NZ);
S.BPHI = zeros(NR,NPHI,NZ);
S.BZ = zeros(NR,NPHI,NZ);
S.FLAG = zeros(NR,NPHI,NZ);
S.JR = zeros(NR,NPHI,NZ);
S.JPHI = zeros(NR,NPHI,NZ);
S.JZ = zeros(NR,NPHI,NZ);
S.VR = zeros(NR,NPHI,NZ);
S.VPHI = zeros(NR,NPHI,NZ);
S.VZ = zeros(NR,NPHI,NZ);
S.DIFFPROF = zeros(NR,NPHI,NZ);
% Fields to be used in HDF5 file
% 
% RAxis = squeeze(X3D(:,:,1));
% ZAxis = squeeze(Z3D(:,:,1));

Rmin = min(min(R_NIMROD));Rmax = max(max(R_NIMROD));
S.R = linspace(0.9*Rmin,1.1*Rmax,NR);

Zmin = min(min(Z_NIMROD));Zmax = max(max(Z_NIMROD));
S.Z = linspace(1.1*Zmin,1.1*Zmax,NZ);

Dphi = 2*pi/NPHI;
S.PHI = zeros(1,NPHI);
for ii=1:NPHI
    S.PHI(ii) = (ii - 1)*Dphi;
end

list_NIMROD = {'Bz3D','Bphi3D','Jr3D','Jphi3D','Jz3D','Vr3D','Vphi3D','Vz3D','diff3D'};
list_outputs = {'BZ','BPHI','JR','JPHI','JZ','VR','VPHI','VZ','DIFFPROF'};

h = figure;
subplot(2,3,1)
surf(R_NIMROD,Z_NIMROD,F.Br3D(:,:,1),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90])
xlabel('$R$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
title('$B_R$','Interpreter','latex')
box on;axis equal;axis([Rmin Rmax Zmin Zmax]);colorbar('location','northoutside')
subplot(2,3,2)
surf(R_NIMROD,Z_NIMROD,F.Bphi3D(:,:,1),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90])
xlabel('$R$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
title('$B_\phi$','Interpreter','latex')
box on;axis equal;axis([Rmin Rmax Zmin Zmax]);colorbar('location','northoutside')
subplot(2,3,3)
surf(R_NIMROD,Z_NIMROD,F.Bz3D(:,:,1),'linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90])
xlabel('$R$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
title('$B_Z$','Interpreter','latex')
box on;axis equal;axis([Rmin Rmax Zmin Zmax]);colorbar('location','northoutside')

for pp=1:NPHI
    A = squeeze(F.Br3D(:,:,pp));
    SI = scatteredInterpolant(reshape(Z_NIMROD,[numel(Z_NIMROD) 1]),reshape(R_NIMROD,[numel(R_NIMROD) 1]),reshape(A,[numel(A) 1]));
    disp(['phi: ' num2str(pp)])
    for rr=1:NR
        for zz=1:NZ
            if checkIfInDomain(S.R(rr),S.Z(zz),S)
                S.BR(rr,pp,zz) = SI(S.Z(zz),S.R(rr));
                S.FLAG(rr,pp,zz) = 1;
            else
                S.FLAG(rr,pp,zz) = 0;
            end
        end
    end    
end
figure(h)
subplot(2,3,4)
D = squeeze(S.BR(:,1,:));
surf(S.R,S.Z,D','linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90])
box on;axis equal;axis([min(S.R) max(S.R) min(S.Z) max(S.Z)]);colorbar('location','northoutside')
xlabel('$R$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
title('Interpolated $B_R$','Interpreter','latex')


for ll=1:numel(list_NIMROD)
    for pp=1:NPHI
        A = squeeze(F.(list_NIMROD{ll})(:,:,pp));
        SI = scatteredInterpolant(reshape(Z_NIMROD,[numel(Z_NIMROD) 1]),reshape(R_NIMROD,[numel(R_NIMROD) 1]),reshape(A,[numel(A) 1]));
        disp([list_NIMROD{ll} ' phi: ' num2str(pp)])
        for rr=1:NR
            for zz=1:NZ
                if S.FLAG(rr,pp,zz) == 1
                    S.(list_outputs{ll})(rr,pp,zz) = SI(S.Z(zz),S.R(rr));
                end
            end
        end
    end
end

figure(h)
subplot(2,3,5)
D = squeeze(S.BPHI(:,1,:));
surf(S.R,S.Z,D','linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90])
box on;axis equal;axis([min(S.R) max(S.R) min(S.Z) max(S.Z)]);colorbar('location','northoutside')
xlabel('$R$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
title('Interpolated $B_\phi$','Interpreter','latex')

figure(h)
subplot(2,3,6)
D = squeeze(S.BZ(:,1,:));
surf(S.R,S.Z,D','linestyle','none');colormap(jet);xlabel('R');ylabel('Z');view([0,90])
box on;axis equal;axis([min(S.R) max(S.R) min(S.Z) max(S.Z)]);colorbar('location','northoutside')
xlabel('$R$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
title('Interpolated $B_Z$','Interpreter','latex')

save(strcat('NIMROD_MST',num2str(id,'%05i')))

%load handel
%sound(y,Fs)
end


function [Ro,Zo,Bo,bi] = getMagneticAxisParams(id)

[r,z,varr]=readcon(id,0); % load data

[R_NIMROD,Z_NIMROD,Br]=onegrid(r,z,varr,14); % Obtain magnetic field components

[~,~,Bz]=onegrid(r,z,varr,15); % Obtain magnetic field components

[~,~,Bphi]=onegrid(r,z,varr,16); % Obtain magnetic field components

B = sqrt(Br.^2 + Bphi.^2 + Bz.^2);

[~,~,n]=onegrid(r,z,varr,36); % Obtain magnetic field components

% Finding the magnetic axis
[A,~]= max(n(:,:,1));
[~,J] = max(A);
[~,I] = max(n(:,J,1));

Ro = R_NIMROD(I,J);
Zo = Z_NIMROD(I,J);
Bo = B(I,J,1);


rb = sqrt((R_NIMROD(end,:)-Ro).^2 + (Z_NIMROD(end,:)-Zo).^2);
tb = atan2(Z_NIMROD(end,:)-Zo,R_NIMROD(end,:)-Ro);
tb(tb<0) = tb(tb<0) + 2*pi;

[~,I] = min(tb);

tb = [tb(I:end) tb(2:I-1)];
rb = [rb(I:end) rb(2:I-1)];

bi = spline(tb,rb);

end

function l = checkIfInDomain(R,Z,S)

r = sqrt((R-S.Ro).^2 + (Z-S.Zo).^2);
t = atan2(Z-S.Zo,R-S.Ro);
if (t<0);t = t + 2*pi;end

if r > ppval(S.bi,t)
    l = 0;
else
    l = 1;
end

end