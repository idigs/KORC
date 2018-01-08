function Psyn = singleParticleSpectrum(Bo,lambda,g,eta)
Psyn = zeros(size(lambda));

c = 2.9979E8; % Speed of light
q = 1.602176E-19; % Electron charge
m = 9.109382E-31; % Electron mass
ep = 8.854E-12;% Electric permitivity

v = c*sqrt(1-1/g^2);

k = q*Bo*sin(eta)/(g*m*v);

l = lambda; 
lc = 4*pi/(3*k*g^3);

z = lc./l;

BK53 = @(x) besselk(5/3,x);
IntBKv = @(nu,x) (pi/sqrt(2))*(1 - 0.25*(4*nu^2 -1))*erfc(sqrt(x)) + ...
    0.25*(4*nu^2 - 1)*sqrt(0.5*pi./x).*exp(-x);

for ii=1:numel(z)
    if (z(ii) < 0.5)
        a = (2.16/2^(2/3))*z(ii)^(1/3);
        Psyn(ii) = integral(BK53,z(ii),a) + IntBKv(5/3,a);
    elseif (z(ii) >= 0.5) && (z(ii) < 2.5)
        a = 0.72*(z(ii) + 1);
        Psyn(ii) = integral(BK53,z(ii),a) + IntBKv(5/3,a);
    else
        Psyn(ii) = IntBKv(5/3,z(ii));
    end
end

Psyn = c*q^2*Psyn./(sqrt(3)*ep*g^2*l.^3);
end