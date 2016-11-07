clear all
close all

sigma_x = 1;
mu_x = 30;
sigma_y = 30;
mu_y = 20;

BI = @(x,y) exp( -0.5*(x-mu_x).^2/sigma_x^2 ).*exp( -0.5*(y-mu_y).^2/sigma_y^2 )/(2*pi*sigma_x*sigma_y);

xmin = 0;
xmax = 50;
ymin = 0;
ymax = 50;

nx = 200;
ny = 200;

x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);


P = zeros(ny,nx);
for ii=1:nx
    for jj=1:ny
        P(jj,ii) = BI(x(ii),y(jj));
    end
end

%% Sampling

num_samples = 5E4;

x_sampled = zeros(1,num_samples);
y_sampled = zeros(1,num_samples);

x_sampled(1) = mu_x;
y_sampled(1) = mu_y;

sigma_rw_x = (xmax - xmin)/20;
sigma_rw_y = (ymax - ymin)/20;


ii=2;
while (ii <= num_samples)
    x_test = x_sampled(ii-1) + random('norm',0,sigma_rw_x);
    while (x_test < xmin) || (x_test > xmax)
        x_test = x_sampled(ii-1) + random('norm',0,sigma_rw_x);
    end
    
    y_test = y_sampled(ii-1) + random('norm',0,sigma_rw_y);
    while (y_test < ymin) || (y_test > ymax)
        y_test = y_sampled(ii-1) + random('norm',0,sigma_rw_y);
    end

%     x_test = random('unif',xmin,xmax);
%     y_test = random('unif',ymin,ymax);
    
    ratio = BI(x_test,y_test)/BI(x_sampled(ii-1),y_sampled(ii-1));
    
    if ( ratio >= 1.0 )
        x_sampled(ii) = x_test;
        y_sampled(ii) = y_test;
        
        ii  = ii + 1;
    elseif (random('uniform',0,1) < ratio)
        x_sampled(ii) = x_test;
        y_sampled(ii) = y_test;
        
        ii  = ii + 1;
    end 
end

%% Figures
h = figure;
subplot(2,1,1)
histogram2(x_sampled,y_sampled,'FaceColor','flat','Normalization','pdf','LineStyle','none')
axis([xmin xmax ymin ymax])
colormap(jet)
xlabel('$x$-axis','Interpreter','latex')
ylabel('$y$-axis','Interpreter','latex')

A = P;
levels = linspace(0,max(max(A)),10);
figure(h);
subplot(2,1,2)
contourf(x,y,A,levels,'ShowText','on')
axis([xmin xmax ymin ymax])
xlabel('$x$-axis','Interpreter','latex')
ylabel('$y$-axis','Interpreter','latex')
box on;
colormap(jet)
hc = colorbar;
ylabel(hc,'$f_{RE}(\chi,p)$','Interpreter','latex','FontSize',16)