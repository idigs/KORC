% Generate data structures for fields2hdf.m from Carlos' data.
R = B.R(1,:);
Z = B.Z(:,1);
BR = B.Br;
BZ = B.Bz;
BPHI = B.Bt;

surf(R,Z,BPHI)

[A,I] = min(B.psirz);
[~,J] = min(A);
b = sqrt(BR.^2+BPHI.^2+BZ.^2);

FLAG = zeros(size(BR));
FLAG(2:end-1,2:end-1) = 1;

fields2hdf(R,[],Z,BR,BPHI,BZ,[],FLAG,'D3D_165139.h5',b(J,I(J)),R(I(J)),Z(J))