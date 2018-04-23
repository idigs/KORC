function ellipseInitialCondition(sf,ro,xo,yo)
% ro = radius of circle
% sf = shear factor
% xo = xo position of radius ro
% yo = yo position of radius ro

narginchk(0,4)


if nargin == 0
    sf = 0.325;
    ro = 0.3;
    xo = 1.675;
    yo = 0.02;
end

plotCurves(sf,ro,xo,yo)

% General quadratic curve
% ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0;

a = 1;
b = -sf;
c = (1 + sf^2);
d = -xo;
f = sf*xo - yo;
g = xo^2 + yo^2 - ro^2;

D = det([a b d; b c f; d f g]);
J = det([a b; b c]);
I = a + c;

major_semi = sqrt( -2*ro^2/( sf*sqrt(sf^2+4) - (2 + sf^2) ) );
minor_semi = sqrt( 2*ro^2/( sf*sqrt(sf^2+4) + (2 + sf^2) ) );

approx_angle = atan2(1,1+sf);

exact_angle = 0.5*acot(0.5*sf);


disp(['D=' num2str(D)])
disp(['J=' num2str(J)])
disp(['I=' num2str(I)])
disp(['D/I=' num2str(D/I)])

disp(['Major semi-axis:' num2str(major_semi)])
disp(['Minor semi-axis:' num2str(minor_semi)])

disp(['Approx angle:' num2str(rad2deg(approx_angle))])
disp(['Exact angle:' num2str(rad2deg(exact_angle))])

end

function plotCurves(sf,ro,xo,yo)

t = linspace(0,2*pi,100);
x = ro*cos(t);
y = ro*sin(t);

X1 = x + sf*y;
Y1 = y;

angle = pi/2 - atan2(1,1+sf);

X2 = X1*cos(angle) - Y1*sin(angle) + xo;
Y2 = X1*sin(angle) + Y1*cos(angle) + yo;

x = x + xo;
y = y + yo;


plot(x,y,'k',X2,Y2,'r','LineWidth',2)
axis equal
end