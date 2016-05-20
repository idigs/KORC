% name = 'fields/xpand_iter3D_sc4_bet015_I87_hi_acc.dat';
% ST = pOrbs(name,'XPANDER','3D',[150,100,150],[1E6,1E-2,10],[-1,1],[7,0,0],[0.985384856627944,0]);

filename = 'ITER_3D.h5'

R = squeeze(B.R(1,:,1));
Z = squeeze(B.Z(:,1,1));
phi=squeeze(B.phi(1,1,:));
h5create(filename, '/NR', [1])
h5write(filename, '/NR',B.NR)
h5create(filename, '/NZ', [1])
h5write(filename, '/NZ',B.NZ)
h5create(filename, '/NPHI', [1])
h5write(filename, '/NPHI',B.Nphi)
h5create(filename, '/R', size(R))
h5write(filename, '/R',R)
h5create(filename, '/Z', size(Z'))
h5write(filename, '/Z',Z')
h5create(filename, '/PHI', size(phi'))
h5write(filename, '/PHI',phi')

dims = [B.NR,B.Nphi,B.NZ];
BR = zeros(dims);
BPHI = zeros(dims);
BZ = zeros(dims);

FR = zeros(fliplr(dims));
FPHI = zeros(fliplr(dims));
FZ = zeros(fliplr(dims));

for ir=1:B.NR
    for iphi=1:B.Nphi
        for iz=1:B.NZ
            BR(ir,iphi,iz) = B.BR(ir,iz,iphi);
            BPHI(ir,iphi,iz) = B.Bphi(ir,iz,iphi);
            BZ(ir,iphi,iz) = B.BZ(ir,iz,iphi);
        end
    end
end

for ir=1:B.NR
    for iphi=1:B.Nphi
        for iz=1:B.NZ
            FR(iz,iphi,ir) = BR(ir,iphi,iz);
            FPHI(iz,iphi,ir) = BPHI(ir,iphi,iz);
            FZ(iz,iphi,ir) = BZ(ir,iphi,iz);
        end
    end
end


fileID = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
dataspaceID = H5S.create_simple(3, dims, dims);
dsetname = '/BR';
datasetID = H5D.create(fileID,dsetname,'H5T_NATIVE_DOUBLE',dataspaceID,'H5P_DEFAULT');
H5D.write(datasetID,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',FR);
H5D.close(datasetID);
H5S.close(dataspaceID);

dataspaceID = H5S.create_simple(3, dims, dims);
dsetname = '/BPHI';
datasetID = H5D.create(fileID,dsetname,'H5T_NATIVE_DOUBLE',dataspaceID,'H5P_DEFAULT');
H5D.write(datasetID,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',FPHI);
H5D.close(datasetID);
H5S.close(dataspaceID);

dataspaceID = H5S.create_simple(3, dims, dims);
dsetname = '/BZ';
datasetID = H5D.create(fileID,dsetname,'H5T_NATIVE_DOUBLE',dataspaceID,'H5P_DEFAULT');
H5D.write(datasetID,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',FZ);
H5D.close(datasetID);
H5S.close(dataspaceID);
H5F.close(fileID);
%%
dims = [10,5,2];

data = zeros(dims);
for ii=1:dims(1)
    for jj=1:dims(2)
        for zz=1:dims(3)
            data(ii,jj,zz) = ii + (jj-1)*dims(1) + (zz-1)*dims(1)*dims(2);
        end
    end
end
data

h5data = zeros(fliplr(dims));
for ii=1:dims(1)
    for jj=1:dims(2)
        for zz=1:dims(3)
            h5data(zz,jj,ii) = data(ii,jj,zz);
        end
    end
end


filename = fullfile('my_file.h5');
fileID = H5F.create(filename,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
dataspaceID = H5S.create_simple(3, dims, dims);
dsetname = '/F';
datasetID = H5D.create(fileID,dsetname,'H5T_NATIVE_DOUBLE',dataspaceID,'H5P_DEFAULT');
H5D.write(datasetID,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',h5data);
H5D.close(datasetID);
H5S.close(dataspaceID);
H5F.close(fileID);
%%
dims = [10,5];
data = zeros(dims);
for ii=1:dims(1)
    for jj=1:dims(2)
        data(ii,jj) = ii + (jj-1)*dims(1);
    end
end
data
h5dims = fliplr(dims)
filename = fullfile('my_file.h5');
fileID = H5F.create(filename,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
dataspaceID = H5S.create_simple(2, dims, dims);
dsetname = '/F';
datasetID = H5D.create(fileID,dsetname,'H5T_NATIVE_DOUBLE',dataspaceID,'H5P_DEFAULT');
H5D.write(datasetID,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',data);
H5D.close(datasetID);
H5S.close(dataspaceID);
H5F.close(fileID);
