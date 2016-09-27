clear all
close all

N = 1E3;

lower_integration_limit = 90;
upper_integration_limit = linspace(1E2,1E3,N);

Q = zeros(1,N);

fun = @(x) besselk(5/3,x);
for ii=1:N
    Q(ii) = integral(fun,lower_integration_limit,upper_integration_limit(ii));
end

DQ = 100*(Q - Q(1))/Q(1);

figure
plot(upper_integration_limit,DQ)
xlabel('Upper integration limit','Interpreter','latex')
ylabel('$\int K_{5/3}(X)$','Interpreter','latex')