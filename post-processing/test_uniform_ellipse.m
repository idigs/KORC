N=1E5;

angle = 30;
m = 0.5;
a = 1;
xo = 1.5;
yo = 0.0;


angle = deg2rad(angle);

r = a*sqrt(rand(1,N));
t = 2*pi*rand(1,N);

x = xo + r.*cos(t);
y = yo + r.*sin(t);

fig = figure;

subplot(2,1,1)
histogram2(x,y,'FaceColor','flat')
colormap(jet)

X = zeros(size(x));
Y = zeros(size(y));

for ii=1:N
    X(ii) = x(ii) + m*y(ii);
    Y(ii) = y(ii);
end

XX = zeros(size(x));
YY = zeros(size(y));

for ii=1:N
    XX(ii) = X(ii)*cos(angle) + Y(ii)*sin(angle);
    YY(ii) = -X(ii)*sin(angle) + Y(ii)*cos(angle);
end

figure(fig)
subplot(2,1,2)
histogram2(XX,YY,'FaceColor','flat')
colormap(jet)

%%
m = 0.325;
a = 0.3;
xo = 1.675;
yo = 0.02;

t = linspace(0,2*pi,100);
x = a*cos(t);
y = a*sin(t);

X1 = x + m*y;
Y1 = y;

angle = pi/2 - atan2(1,1+m);
% angle = pi - atan2(a,a*m) - atan2(1,1+m);

X2 = X1*cos(angle) - Y1*sin(angle) + xo;
Y2 = X1*sin(angle) + Y1*cos(angle) + yo;

x = x + xo;
y = y + yo;

% figure;
% plot(x,y,'k',X1,Y1,'r',X2,Y2,'b')
% axis equal

% hold on;
plot(x,y,'k',X2,Y2,'r','LineWidth',2)
axis equal;
