% Rejection method for a thermal distribution in three dimensions of
% velovity space.

clear all
close all
clc

c = 2.9979E8; % Speed of light
qe = 1.602176E-19; % Electron charge
me = 9.109382E-31; % Electron mass

NV = 500;

Npcls = 5E5;

Eo = 10E3; % In eV
Vth = sqrt(Eo*qe/me); % Thermal velocity

% Below this line all quantities are dimensionless
Vth = Vth;

vmax = 0.8*c;
v = linspace(0,vmax,NV);

f = @(x,y,z) (1/Vth^3)*exp(-0.5*(x.^2 + y.^2 + z.^2)/Vth^2)/(2*pi)^1.5;
%% Sampling of distribution
Nsamples = 2*Npcls;

vx = zeros(1,Nsamples);
vy = zeros(1,Nsamples);
vz = zeros(1,Nsamples);

sv = 0.05*c;

ii=2;
while (ii <= Nsamples)
    x = vx(ii-1) + random('norm',0,sv);
    y = vy(ii-1) + random('norm',0,sv);
    z = vz(ii-1) + random('norm',0,sv);
    
    ratio = f(x,y,z)/f(vx(ii-1),vy(ii-1),vz(ii-1));
    
    if ( ratio >= 1.0 )
        vx(ii) = x;
        vy(ii) = y;
        vz(ii) = z;
        ii  = ii + 1;
    elseif (ratio > random('uniform',0,1))
        vx(ii) = x;
        vy(ii) = y;
        vz(ii) = z;
        ii  = ii + 1;
    end
end

figure
hold on
histogram(vx,'Normalization','pdf','LineStyle','none')
hold off


