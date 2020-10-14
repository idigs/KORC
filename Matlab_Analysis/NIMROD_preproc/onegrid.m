function [R,Z,fld]=onegrid(r,z,varr,ivar)
%function [R,Z,fld]=onegrid(r,z,varr,ivar)
%puts the NIMROD data stored in sereval blocks on to one large grid
%inputs r,z,varr are outputs from readcon.m, ivar is the number of the
%variable(s) to be output as fld (see readcon.m). The R,Z outputs are the 
%grid for fld.  

if ivar<32 
iv1=ivar+3;
elseif ivar<56
iv1=ivar+1;
else
iv1=ivar+(size(varr,5)-55)/2;
end

mxbl=size(varr,1)-1;
mybl=size(varr,2)-1;

for ir=1:size(r,3)
    if (r(end,1,1)==r(1,1,ir))
        nybl=ir-1;
        nxbl=size(r,3)/nybl;
        break
    end
end

for ibl=1:nxbl
for jbl=1:nybl
R((ibl-1)*mxbl+1:(ibl-1)*mxbl+(mxbl+1),(jbl-1)*mybl+1:(jbl-1)*mybl+(mybl+1))=r(:,:,nybl*(ibl-1)+jbl);
Z((ibl-1)*mxbl+1:(ibl-1)*mxbl+(mxbl+1),(jbl-1)*mybl+1:(jbl-1)*mybl+(mybl+1))=z(:,:,nybl*(ibl-1)+jbl);
fld((ibl-1)*mxbl+1:(ibl-1)*mxbl+(mxbl+1),(jbl-1)*mybl+1:(jbl-1)*mybl+(mybl+1),:)=...
      varr(:,:,:,nybl*(ibl-1)+jbl,ivar)+i*varr(:,:,:,nybl*(ibl-1)+jbl,iv1);
end
end
