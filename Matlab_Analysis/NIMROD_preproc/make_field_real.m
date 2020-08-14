function [X3,Y3,Z3,rfield]=make_field_real(R,Z,field,nplanes);
nfour=size(field,3);

x = input('prompt');

rfield=zeros(size(R,1),size(R,2),nplanes+1);

for jp=1:nplanes+1
  phi=(jp-1)*2*pi/nplanes;
  X3(:,:,jp)=R(:,:).*cos(-phi);
  Y3(:,:,jp)=R(:,:).*sin(-phi);
  Z3(:,:,jp)=Z(:,:);
  for jm=1:nfour
    rfield(:,:,jp)=rfield(:,:,jp)+2*real(field(:,:,jm)*exp(i*phi*(jm-1)));
  end
    rfield(:,:,jp)=rfield(:,:,jp)-field(:,:,1);
end
    

