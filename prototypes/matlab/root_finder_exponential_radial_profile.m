close all
clear all

relerr = 1E-15;

g = 28.4; % gamma factor
P = rand; %rand;
a = 0.65;
k = @(x) 1.0;%2.0768 - 16.945./x;
f = @(r,a,g,P) exp(-k(g)*r).*(1+r.*k(g)) + (1-exp(-k(g).*a).*(1+a.*k(g))).*P - 1;


%%
% 
% r = linspace(0,a,500);
% fig = figure;
% plot(r,f(r,a,g,P),'b')
% grid minor
% 
% ro = 0;
% rn = a;
% 
% for ii=1:10
%     fo = f(ro,a,g,P);
%     fn = f(rn,a,g,P);
%     
%     figure(fig)
%     hold on;
%     plot([ro rn],[fo fn],'ro-','MarkerFaceColor',[0.5,0,0]);
%     hold off
%     
%     rt = rn;
%     rn = rn - fn*(rn - ro)/(fn - fo);
%     
%     figure(fig)
%     hold on;
%     plot(rn,0,'go','MarkerFaceColor',[0,0.5,0]);
%     hold off
%     
%     ro = rt;
% end
% 
% figure(fig)
% hold on;
% plot(rn,fn,'mo','MarkerFaceColor',[1,0,1]);
% hold off

%%
rl = 0;
rm = 0;
rr = a;

fl = f(rl,a,g,P);
fr = f(rr,a,g,P);

err = abs(fr - fl)

while (err > relerr)    
    rm = 0.5*(rr-rl)+rl;
    fm = f(rm,a,g,P);
    
    if (sign(fr) == sign(fl))
        R = [rl rr];
        F = [fl fr];
        [~,I] = min(F);
        disp(['r=' num2str(R(I))])
        break
    end
    
    if (sign(fm) == sign(fr))
        rr = rm;
    else
        rl = rm;
    end
    
    fl = f(rl,a,g,P);
    fr = f(rr,a,g,P);
    err = abs(fr - fl)
end

r = linspace(0,a,500);
plot(r,f(r,a,g,P),'m',rm,fm,'ob');
grid minor;box on;