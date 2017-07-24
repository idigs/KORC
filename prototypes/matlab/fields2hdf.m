function fields2hdf(R,PHI,Z,BR,BPHI,BZ,F,outputfile,Bo,Ro)
% size(A) = [numel(R),numel(PHI),numel(Z)], where A can be any of the field
% components or magnetic flux.
% Example for JFIT D3D fields:
% fields2hdf(R,PHI,Z,[],[],[],D,outputfile,2.19,1.695)
% Example for EFIT axisymmetric D3D magnetic fields
% fields2hdf(R,[],Z,BR,BPHI,BZ,[],'D3D_Q.h5')
% Example for ITER fields using XPANDER fields
% fields2hdf(R,PHI,Z,BR,BPHI,BZ,[],'ITER.h5')

narginchk(8,10)

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

if ~isempty(F)
    dsetname = '/PSIp';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,F)
end

if ~isempty(BR)
    dsetname = '/BR';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,BR)
end

if ~isempty(BR)
    dsetname = '/BPHI';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,BPHI)
end

if ~isempty(BR)
    dsetname = '/BZ';
    h5create(outputfile,dsetname,dims)
    h5write(outputfile,dsetname,BZ)
end

if nargin > 8
    dsetname = '/Bo';
    h5create(outputfile,dsetname,[1])
    h5write(outputfile,dsetname,Bo)
    
    dsetname = '/Ro';
    h5create(outputfile,dsetname, [1])
    h5write(outputfile,dsetname,Ro)
end
end