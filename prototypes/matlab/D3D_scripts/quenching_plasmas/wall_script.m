clear all

load('tmp.mat')

da = 0.05;
Ro = 1.35;
N = 20;

n = numel(r);

va = [];
vr = [];
vR = [];
vZ = [];


for ii=1:n-1
    if( (a(ii+1)-a(ii))> da )
        DR = R(ii+1) - R(ii);
        DZ = Z(ii+1) - Z(ii);
        
        if (abs(DR)>1E-15)
            m = DZ/DR;
            yo = Z(ii); xo = R(ii);
            x = linspace(R(ii),R(ii+1),N);
            y = m*(x - xo) + yo;
        else
            y = linspace(Z(ii),Z(ii+1),N);
            x = R(ii)*ones(1,N);
        end
    else
        x = R(ii);
        y = Z(ii);
    end
    vR = [vR x];
    vZ = [vZ y];
    hold on;plot(x,y,'o');hold off
end

r = sqrt((vR-Ro).^2 + vZ.^2);
a = atan2(vZ,vR-Ro);
a(a<0) = a(a<0) + 2*pi;

[a,ia,ic] = unique(a,'stable');
r = r(ia);
vR = vR(ia);
vZ = vZ(ia);

hold on;plot(vR,vZ,'k.-');hold off
figure;plot(a,r)