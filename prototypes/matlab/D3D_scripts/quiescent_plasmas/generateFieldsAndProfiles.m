clear all
close all
clc

% load('sh165139_06012_Bs.mat');
load('sh165826_05635_Bs.mat');

timeB = B.time;

R = G.rg;
Z = G.zg;

NR = numel(R);
NZ = numel(Z);

BR = B.Br; % NZ x NR
BPHI = B.Bt;  % NZ x NR
BZ = B.Bz; % NZ x NR
FLUX = G.psirz; % NZ x NR

Z0 = G.zmaxis;
R0 = G.rmaxis;

E0 = 0.6/(2*pi*R0);
B0 = G.bzero;

clear B G

load('wall.mat')

fig = figure;
FLAG = zeros(NR,NZ);
for rr=1:NR
    for zz=1:NZ
        a = atan2(Z(zz),R(rr)-Ro);
        if (a<0); a=a+2*pi; end
        r = sqrt( (R(rr)-Ro).^2 + Z(zz).^2 );
        
        [~,ia] = min(abs(a-angle));
        
        if ( (ia>1) && (ia~=numel(angle)) )
            d = a - angle(ia);
            dp = a - angle(ia+1);
            dm = a - angle(ia-1);
        elseif (ia == numel(angle))
            d = sign(a - angle(ia));
            dp = sign(d);
            dm = -sign(d);
        elseif (ia == 1)
            d = sign(a - angle(ia));
            if (sign(d) ~= 0)
                dp = -sign(d);
                dm = sign(d);
            else
                dp = 1;
            end
        end
        
        if ( sign(d) == sign(dp) )
            m = (radius(ia) - radius(ia-1))/(angle(ia) - angle(ia-1));
            y = m*(a - angle(ia-1)) + radius(ia-1);
        else
            m = (radius(ia+1) - radius(ia))/(angle(ia+1) - angle(ia));
            y = m*(a - angle(ia)) + radius(ia);
        end
        
        if (r > y)
            FLAG(rr,zz) = 0;
            figure(fig)
            hold on;plot(R(rr),Z(zz),'r.');hold off
        else
            FLAG(rr,zz) = 1;
            figure(fig)
            hold on;plot(R(rr),Z(zz),'k.');hold off
        end
    end
end

clear Ro a angle r radius d dm dp ia m rr zz y

load('separatrix.mat')

RSPX=flip(RSPX);
ZSPX=flip(ZSPX);

rspx = sqrt( (RSPX-R0).^2 + (ZSPX-Z0).^2 );
aspx = atan2(ZSPX-Z0,RSPX-R0);
aspx(aspx<0) = aspx(aspx<0) + 2*pi;
[~,I] = min(aspx);
aspx = [aspx(I:end)' aspx(1:I-1)'];
rspx = [rspx(I:end)' rspx(1:I-1)'];

aspx(1) = 0;
aspx = [aspx 2*pi];
rspx = [rspx rspx(1)];

[aspx,ia,ic] = unique(aspx,'stable');
rspx = rspx(ia);

% load('eq165139.mat')
load('eq165826.mat')

rho = eq.oneD.rho(1,:);
time = eq.oneD.t(:,1);

[~,it] = min(abs(timeB - time));

ne = eq.oneD.raw.ne(it,:);
Te = eq.oneD.raw.Te(it,:);
Zeff = eq.oneD.final.Zeff(it,:);

ZEFF = zeros(NR,NZ);
NE = zeros(NR,NZ);
TE = zeros(NR,NZ);
for rr=1:NR
    for zz=1:NZ
        a = atan2(Z(zz)-Z0,R(rr)-R0);
        if (a<0); a=a+2*pi; end
        r = sqrt( (R(rr)-R0).^2 + (Z(zz)-Z0).^2 );
        
        rspx_i = interp1(aspx,rspx,a,'pchip');
        
        if (r < rspx_i)
            ZEFF(rr,zz) = interp1(rho,Zeff,r/rspx_i,'pchip');
            NE(rr,zz) = interp1(rho,ne,r/rspx_i,'pchip');
            TE(rr,zz) = interp1(rho,Te,r/rspx_i,'pchip');
            figure(fig)
            hold on;plot(R(rr),Z(zz),'mo');hold off
        else
            ZEFF(rr,zz) = Zeff(end);
            NE(rr,zz) = ne(end);
            TE(rr,zz) = Te(end);
            figure(fig)
%             hold on;plot(R(rr),Z(zz),'mx');hold off
        end        
    end
end


% outputfile = 'D3D_165139_flux.h5';
outputfile = 'D3D_165826_flux.h5';
fields2hdf(R,[],Z,BR',BPHI',BZ',[],[],[],FLUX,FLAG,NE,TE,ZEFF,outputfile,B0,E0,R0,Z0)