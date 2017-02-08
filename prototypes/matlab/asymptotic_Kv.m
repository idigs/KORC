% Asymptotic expression of Kv(x).
clear all
close all

x = linspace(0.1,500,1E5);

bk13 = @(x) besselk(1/3,x);
bk23 = @(x) besselk(2/3,x);
bk53 = @(x) besselk(5/3,x);

fx = @(x) sqrt(pi/2)*exp(-x)./sqrt(x);

f = @(v,x) sqrt(pi/2)*(1 + 0.125*(4*v^2-1)./x).*exp(-x)./sqrt(x);

f13 = f(1/3,x);
f23 = f(2/3,x);
f53 = f(5/3,x);

h = figure;
subplot(3,1,1)
loglog(x,abs(1-bk13(x)./fx(x)),'b',x,abs(1-bk23(x)./fx(x)),'r',x,abs(1-bk53(x)./fx(x)),'g')
axis([0 max(x) 0 1])

subplot(3,1,2)
loglog(x,abs(1-bk13(x)./f(1/3,x)),'b',x,abs(1-bk23(x)./f(2/3,x)),'r',x,abs(1-bk53(x)./f(5/3,x)),'g')
axis([0 max(x) 0 1])

Tol = 10E-2;

I13 = find(abs(1 - bk13(x)./f(1/3,x)) < Tol,1);
I23 = find(abs(1 - bk23(x)./f(2/3,x)) < Tol,1);
I53 = find(abs(1 - bk53(x)./f(5/3,x)) < Tol,1);

%% Integrals
N = 1E3;
infinity = 1000;
amax = 20;

a = linspace(x(I13),amax,N);
Intbk = zeros(1,N);
IntAsymptotic = zeros(1,N);

I = @(v,z) pi*(1 - 0.25*(4*v.^2 - 1)).*(1 - erf(sqrt(z)))/sqrt(2) + ...
    0.25*(4*v.^2 - 1)*sqrt(0.5*pi./z)*exp(-z);

for ii=1:N
    Intbk(ii) = integral(bk13,a(ii),infinity);
    IntAsymptotic(ii) = I(1/3,a(ii));
end

figure(h);
subplot(3,1,3)
semilogy(a,abs(1 - Intbk./IntAsymptotic),'b')

a = linspace(x(I23),amax,N);
Intbk = zeros(1,N);
IntAsymptotic = zeros(1,N);

I = @(v,z) pi*(1 - 0.25*(4*v.^2 - 1)).*(1 - erf(sqrt(z)))/sqrt(2) + ...
    0.25*(4*v.^2 - 1)*sqrt(0.5*pi./z)*exp(-z);

for ii=1:N
    Intbk(ii) = integral(bk23,a(ii),infinity);
    IntAsymptotic(ii) = I(2/3,a(ii));
end

figure(h);
subplot(3,1,3)
hold on
semilogy(a,abs(1 - Intbk./IntAsymptotic),'r')
hold off

a = linspace(x(I53),amax,N);
Intbk = zeros(1,N);
IntAsymptotic = zeros(1,N);

I = @(v,z) pi*(1 - 0.25*(4*v.^2 - 1)).*(1 - erf(sqrt(z)))/sqrt(2) + ...
    0.25*(4*v.^2 - 1)*sqrt(0.5*pi./z)*exp(-z);

for ii=1:N
    Intbk(ii) = integral(bk53,a(ii),infinity);
    IntAsymptotic(ii) = I(5/3,a(ii));
end

figure(h);
subplot(3,1,3)
hold on
semilogy(a,abs(1 - Intbk./IntAsymptotic),'g')
hold off

%% Spectral density of radiated power
N = 1E3;
amax = 10;
amin = 1E-3;
a = linspace(amin,amax,N);
P = zeros(N,N);
ac = zeros(1,N);
Tol = 1E-2;

I = @(v,z) pi*(1 - 0.25*(4*v.^2 - 1)).*(1 - erf(sqrt(z)))/sqrt(2) + ...
    0.25*(4*v.^2 - 1)*sqrt(0.5*pi./z)*exp(-z);

for jj=1:N-1
    for ii=jj+1:N
        P(jj,ii) = trapz(a(jj:ii),bk53(a(jj:ii))) + I(5/3,a(ii));
    end
    D = abs(1 - P(jj,jj+1:end)./P(jj,end));
    [~,index] = find(D < Tol,1);
    ac(jj) = a(index+jj);
end
