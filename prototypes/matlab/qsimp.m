function I = qsimp(a,b)
N = 4;
Tol = 1E-5;

BK53 = @(x) besselk(5/3,x);

h = b - a;
I = 0.5*(BK53(a) + BK53(b));

ii=1;
flag = true;
while flag
    Iold = h*I;
    
    ii = ii + 1;
    npoints = 2^(ii-2);
    h = 0.5*(b-a)/npoints;
    
    for jj=1:npoints
        z = a + h + 2*(jj-1)*h;
        I = I + BK53(z);
    end
    
    Inew = h*I;
    
    err = abs((Inew - Iold)/Iold);
    flag = ~(err < Tol);
end

I = h*I;

disp(['Integral: ' num2str(I) ' Iterations: ' num2str(ii)])

IMatlab = integral(BK53,a,b);
reldiff = 100*abs((IMatlab - I)/IMatlab);

disp(['Difference: ' num2str(reldiff)])
end