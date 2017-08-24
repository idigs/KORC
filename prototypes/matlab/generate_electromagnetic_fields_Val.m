% Script to generate the electromagnetic fields of Val's simulations of
% thermal quench.

close all
clear all

Npplanes = 36;
NR = 150;
NZ = 150;
Nphi = Npplanes;


t = [     25         0.0041 
          50	     0.0196
         100	     0.0696
         300	     0.2696
         500	     0.4696
         700	     0.6696
         800	     0.7051
         900	     0.7345
        1000	     0.7674
        1050	     0.7799
        1100	     0.8008
        1200	     0.8785
        1300	     0.9645
        1350	     1.0145
        1500	     1.1645
        1700	     1.3645
        1900	     1.5645];% Diverted plasma simulation

[r,z,varr]=readcon(25,0); % load data

[R,Z,Br]=onegrid(r,z,varr,14); % Obtain magnetic field components

[~,~,Bz]=onegrid(r,z,varr,15); % Obtain magnetic field components

[~,~,Bphi]=onegrid(r,z,varr,16); % Obtain magnetic field components

[~,~,n]=onegrid(r,z,varr,36); % Obtain magnetic field components

[~,~,n1]=onegrid(r,z,varr,57); 


[X3D,Y3D,Z3D,Br3D]=make_field_real(R,Z,Br,Npplanes);
X3D(:,:,end) = [];Y3D(:,:,end) = [];Z3D(:,:,end) = [];
Br3D(:,:,end) = [];
X3D = flip(X3D,3);Y3D = flip(Y3D,3);Z3D = flip(Z3D,3);
Br3D = flip(Br3D,3);

[~,~,~,Bz3D]=make_field_real(R,Z,Bz,Npplanes);
Bz3D(:,:,end) = [];
Bz3D = flip(Bz3D,3);

[~,~,~,Bphi3D]=make_field_real(R,Z,Bphi,Npplanes);
Bphi3D(:,:,end) = [];
Bphi3D = flip(Bphi3D,3);

B = sqrt(Br3D.^2 + Bphi3D.^2 + Bz3D.^2);

[~,~,~,n3D]=make_field_real(R,Z,n,Npplanes);
n3D(:,:,end) = [];
n3D = flip(n3D,3);

[~,~,~,n13D]=make_field_real(R,Z,n,Npplanes);
n13D(:,:,end) = [];
n13D = flip(n13D,3);



% SI = scatteredInterpolant(reshape(X3D,[numel(X3D) 1]),reshape(Z,[numel(Z) 1]),reshape(Bslice,[numel(Bslice) 1]));

