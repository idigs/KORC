% Generate data structures for fields2hdf.m from Carlos' data.
close all
R = B.R(1,:);
Z = B.Z(:,1);
BR = B.Br;
BZ = B.Bz;
BPHI = B.Bt;

subplot(1,3,1)
surf(R,Z,BR,'LineStyle','none')
colormap(jet);colorbar
view([0 90])
xlabel('R')
ylabel('Z')

subplot(1,3,2)
surf(R,Z,BPHI,'LineStyle','none')
colormap(jet);colorbar
view([0 90])
xlabel('R')
ylabel('Z')

subplot(1,3,3)
surf(R,Z,BZ,'LineStyle','none')
colormap(jet);colorbar
view([0 90])
xlabel('R')
ylabel('Z')

[A,I] = min(B.psirz);
[~,J] = min(A);
b = sqrt(BR.^2+BPHI.^2+BZ.^2);

FLAG = zeros(size(BR));
FLAG(2:end-1,2:end-1) = 1;

fields2hdf(R,[],Z,BR',BPHI',BZ',[],FLAG','D3D_165139.h5',b(J,I(J)),R(I(J)),Z(J))