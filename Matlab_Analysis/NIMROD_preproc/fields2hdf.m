function fields2hdf(R,PHI,Z,BR,BPHI,BZ,ER,EPHI,EZ,FLUX,FLAG,NE,TE,ZEFF,outputfile,Bo,Eo,Ro,Zo)
% size(A) = [numel(R),numel(PHI),numel(Z)], where A can be any of the field
% size(A) = [numel(R),numel(Z)], where A can be any of the field
% components or magnetic flux.
% Example for JFIT D3D fields:
% fields2hdf(R,[],Z,[],[],[],D,[],outputfile,2.19,1.695,0.0)
% Example for EFIT axisymnumelmetric D3D magnetic fields
% fields2hdf(R,[],Z,BR,BPHI,BZ,[],FLAG,'D3D_165826_1.h5',1.4096,1.7033,0.025)
% Example for ITER fields using XPANDER fields
% fields2hdf(R,PHI,Z,BR,BPHI,BZ,[],[],'ITER.h5')
% Example for DIIID fields using VMEC fields
% fields2hdf(R,[],Z,BR,BPHI,BZ,[],FLAG,'D3D_VMEC.h5',2.4538,1.6282,0.0282)
% Example for NIMROD fields
% fields2hdf(R,PHI,Z,BR,BPHI,BZ,[],[],'NIMROD_DIVERTED.h5',2.1170,1.7272,0.0142)
% fields2hdf(S.R,S.PHI,S.Z,S.BR,S.BPHI,S.BZ,[],S.FLAG,'NIMROD_DIVERTED_1100D.h5',S.Bo,S.Ro,S.Zo)

narginchk(9,19)

NR = numel(R);
NPHI = numel(PHI);
NZ = numel(Z);


if ~isempty(PHI)
    

    
    dims = [NR,NPHI,NZ];
    
    dsetname = '/NPHI';
    h5create(outputfile,dsetname, [1])
    h5write(outputfile,dsetname,dims(2))
    
    dsetname = '/PHI';
    h5create(outputfile,dsetname, dims(2))
    h5write(outputfile,dsetname,PHI)
    
    dsetname = '/NZ';
    h5create(outputfile,dsetname, [1])
    h5write(outputfile,dsetname,dims(3))
    
    dsetname = '/Z';
    h5create(outputfile,dsetname,dims(3))
    h5write(outputfile,dsetname',Z)
else
    dims = [NR,NZ];
    
    dsetname = '/NZ';
    h5create(outputfile,dsetname, [1])
    h5write(outputfile,dsetname,dims(2))
    
    dsetname = '/Z';
    h5create(outputfile,dsetname,dims(2))
    h5write(outputfile,dsetname',Z)
end

dsetname = '/NR';
h5create(outputfile,dsetname, [1])
h5write(outputfile,dsetname,dims(1))

dsetname = '/R';
h5create(outputfile,dsetname, dims(1))
h5write(outputfile,dsetname,R)

if ~isempty(FLUX)
    dsetname = '/PSIp';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,FLUX)
end

if ~isempty(NE)
    dsetname = '/ne';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,NE)
end

if ~isempty(TE)
    dsetname = '/Te';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,TE)
end

if ~isempty(ZEFF)
    dsetname = '/Zeff';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,ZEFF)
end

if ~isempty(BR)
    dsetname = '/BR';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,BR)
end

if ~isempty(BPHI)
    dsetname = '/BPHI';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,BPHI)
end

if ~isempty(BZ)
    dsetname = '/BZ';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,BZ)
end

if ~isempty(ER)
    dsetname = '/ER';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,ER)
end

if ~isempty(EPHI)
    dsetname = '/EPHI';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,EPHI)
end

if ~isempty(EZ)
    dsetname = '/EZ';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,EZ)
end

if ~isempty(FLAG)
    dsetname = '/FLAG';
    h5create(outputfile,dsetname,dims,'DataType','int8')
    h5write(outputfile,dsetname,int8(FLAG))
end

if nargin > 8
    dsetname = '/Bo';
    h5create(outputfile,dsetname,[1])
    h5write(outputfile,dsetname,Bo)
    
    dsetname = '/Eo';
    h5create(outputfile,dsetname,[1])
    h5write(outputfile,dsetname,Eo)
    
    dsetname = '/Ro';
    h5create(outputfile,dsetname, [1])
    h5write(outputfile,dsetname,Ro)
    
    dsetname = '/Zo';
    h5create(outputfile,dsetname, [1])
    h5write(outputfile,dsetname,Zo)
end

end