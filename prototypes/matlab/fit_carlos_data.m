load('pdf_re.mat');
fo = fitoptions('Method','NonlinearLeastSquares');
ft = fittype('x^(k-1)*exp(-x/t)','independent','x')
w = 1./err;w = w/sum(w);
[fitobj, gof] = fit(xd',10.^yd',ft,'Weights',w)
confint = predint(fitobj,xd,0.95,'functional','on');
% figure;
hold on;
plot(fitobj,xd,10.^yd);
plot(xd,confint,'m--')
hold off
box on
ylabel('$\log_{10}f_e$','Interpreter','latex','FontSize',16)
xlabel('$\mathcal{E}$ (MeV)','Interpreter','latex','FontSize',16)