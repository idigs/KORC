clear all
close all
clc

nbins = 40;
s = [1, 10, 20];
m = [0, -10, 17];

h = figure;
rn = figure;

for ii=1:numel(s)
    f = random('norm',m(ii),s(ii),1,5000);
    
    x = linspace(min(f),max(f),nbins)
    fh = hist(f,x)
    
    dx = mean(diff(x));
    fh = fh/(sum(fh)*dx);
    
    fx = exp( -(x - m(ii)).^2/(2*s(ii)^2) )/(s(ii)*sqrt(2*pi));
    
    figure(h)
    hold on
    plot(x,fh,'s',x,fx,'k')
    hold off

    y = (x - m(ii))/s(ii);
    frn = s(ii)*fh;

    figure(rn)
    hold on
    plot(y,frn,'o-')
    hold off
    
end

figure(h)
box on
grid on

z = linspace(-6,6,100);
fz = exp( -0.5*z.^2 )/sqrt(2*pi);

figure(rn)
hold on
plot(z,fz)
hold off
grid on
box on