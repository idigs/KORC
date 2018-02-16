% Other code lines
clear all
% close all

BfieldsTime = 6.0120;
load('eq165139.mat')

% load('eq165826.mat')
% BfieldsTime = 5.635;

r = eq.oneD.rho(1,:);
time = eq.oneD.t(:,1);

[~,I] = min(abs(BfieldsTime - time));

ne1D = eq.oneD.raw.ne(I,:);
ne0D = eq.zeroD.raw.ne(I);
Te1D = eq.oneD.raw.Te(I,:);
Te0D = eq.zeroD.raw.Te;
Zeff1D = eq.oneD.final.Zeff(I,:);
Zeff0D = eq.zeroD.final.Zeff(I);
E_over_ECH1D = eq.oneD.final.E_over_ECH(I,:);
C = 2*E_over_ECH1D./(Zeff1D + 1);
E_CH1D = eq.oneD.calc.E_CH(I,:);
E_CH0D = eq.zeroD.calc.E_CH;
E_over_ECH0D = eq.zeroD.calc.E_phi(I)/eq.zeroD.calc.E_CH(I);
BT = abs(eq.zeroD.raw.BT0);
ET = eq.zeroD.calc.E_phi;
ET1D = E_CH1D.*E_over_ECH1D;
Tauc1D = eq.oneD.calc.tau_coll(I,:);

fig=figure;
subplot(2,3,1)
yyaxis left
plot(r,ne1D,'b',0,ne0D,'b*')
grid minor
ylabel('$n_e$ (m$^{-3}$)','Interpreter','latex');
xlabel('$r/a$','Interpreter','latex');

yyaxis right
plot(r,Te1D,'r',0,Te0D(I),'r*')
ylabel('$T_e$ (eV)','Interpreter','latex');
xlabel('$r/a$','Interpreter','latex');


subplot(2,3,2)
yyaxis left
plot(r,Zeff1D,'b',0,Zeff0D,'b*');
grid minor
ylabel('$Z_{eff}$','Interpreter','latex');
xlabel('$r/a$','Interpreter','latex');

yyaxis right
plot(r,ET1D,'r');
ylabel('$E_{\phi}$ (V/m)','Interpreter','latex');
xlabel('$r/a$','Interpreter','latex');


subplot(2,3,3)
yyaxis left
plot(time,E_CH0D,'b',time(I),E_CH0D(I),'b*');
ylabel('$E_{CH}$ (V/m)','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');

yyaxis right
plot(time,Te0D,'r',time(I),Te0D(I),'r*');
grid minor
ylabel('$T_e$ (eV)','Interpreter','latex');


subplot(2,3,4)
yyaxis left
plot(r,E_over_ECH1D,'b',0,E_over_ECH0D,'b*');
ylabel('$E/E_{CH}$','Interpreter','latex');
xlabel('$r/a$','Interpreter','latex');

yyaxis right
plot(r,E_CH1D,'r',0,E_CH0D(I),'r*');
grid minor
ylabel('$E_{CH}$ (V/m)','Interpreter','latex');
xlabel('$r/a$','Interpreter','latex');


subplot(2,3,5)
yyaxis left
plot(r,C,'b');
grid minor
ylabel('$2\bar{E}/(Z_{eff} + 1)$','Interpreter','latex');
xlabel('$r/a$','Interpreter','latex');

yyaxis right
plot(r,Tauc1D,'r');
ylabel('$\tau_{coll}$','Interpreter','latex');


subplot(2,3,6)
yyaxis left
plot(time,BT,time(I),BT(I),'r*')
ylabel('$B_{\phi}$ (T)','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');

yyaxis right
plot(time,ET,time(I),ET(I),'r*')
grid minor
ylabel('$E_{\phi}$ (V/m)','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
