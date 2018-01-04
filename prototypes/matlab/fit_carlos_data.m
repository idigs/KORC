% close all
% load('pdf_re.mat');

% xfit = xd';
% yfit = 10.^yd';

% xfit = xd(yd~=0);
% yfit = yd(yd~=0);

% w = 1./err;w = w/sum(w);
fo = fitoptions('Method','NonlinearLeastSquares');

% ft = fittype('a*x^(k-1)*exp(-x/t)/(gamma(k)*t^k)','independent','x','problem','k');
% [fitobj, gof] = fit(xfit,yfit,ft,'Weights',w,'StartPoint',[1E6,0.5],'Lower',[1 0.1],'TolX',1E-8,'problem',0.8)

% ft = fittype('a*x^(k-1)*exp(-x/t)/(gamma(k)*t^k)','independent','x');
% [fitobj, gof] = fit(xfit,yfit,ft,'Weights',w,'StartPoint',[2E5,1,0.5],'Lower',[1 0.1 0.1],'TolX',1E-8)

ft = fittype('a*x^(k-1)*exp(-x/t)/(gamma(k)*t^k)','independent','x');
[fitobj, gof] = fit(xfit,yfit,ft,'StartPoint',[7E10,8,0.5],'Lower',[1 0.01 0.01],'TolX',1E-9)
% 
% ft = fittype('a*x^(k-1)*exp(-x/t)/(gamma(k)*t^k)','independent','x');
% [fitobj, gof] = fit(xfit,yfit,ft,'StartPoint',[1E2,2,0.5],'TolX',1E-9)

% ft = fittype('a*x^(k-1)*exp(-x/t)','independent','x');
% [fitobj, gof] = fit(xfit,yfit,ft,'Weights',w,'StartPoint',[1E4,1.0,0.5])

confint = predint(fitobj,xfit,0.95,'functional','on');
figure;
hold on;
% uerr = abs(10.^(yd+err) - 10.^yd);
plot(fitobj,xfit,yfit);
% errorbar(xfit,yfit,uerr)
plot(xfit,confint,'m--')
ax = gca;
ax.YAxis.Scale = 'log';
hold off
box on
grid on
grid minor
ylabel('$\log_{10}f_e$','Interpreter','latex','FontSize',16)
xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)

% fgamma = @(x,k,t) x.^(k-1).*exp(-x/t)./(gamma(k)*t^k);