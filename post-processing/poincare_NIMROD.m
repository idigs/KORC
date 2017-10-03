close all
% clear all

ST = diagnoseKORC('../KORC-FO/outputFiles/','on',[0,7000])
X = ST.data.sp1.X;
t1 = squeeze(atan2(X(2,:,:),X(1,:,:)));
t1(t1<0) = t1(t1<0) + 2*pi;
a1 = abs(diff(t1,1,2));

t2 = squeeze(atan2(X(2,:,:),X(1,:,:)));
a2 = abs(diff(t2,1,2));

nume = size(t1,1);

% R1 = [];
% Z1 = [];
% R2 = [];
% Z2 = [];
% 
% h = figure;
% for ii=1:nume
%     I = find(a1(:,ii)<0);
%     R1 = [R1; squeeze(sqrt(X(1,ii,I).^2 + X(2,ii,I).^2))];
%     Z1 = [Z1; squeeze(X(3,ii,I))];
%     
%     I = find(a2(:,ii)<0);
%     R2 = [R2; squeeze(sqrt(X(1,ii,I).^2 + X(2,ii,I).^2))];
%     Z2 = [Z2; squeeze(X(3,ii,I))];
% end
% figure(h)
% subplot(1,2,1)
% hold on;
% plot3(R1,zeros(size(R1)),Z1,'k.','MarkerSize',4);
% hold off;axis equal
% 
% subplot(1,2,2)
% hold on;
% plot3(R2,zeros(size(R2)),Z2,'k.','MarkerSize',4);
% hold off;axis equal

h = figure;
for ii=1:nume
    I = find(a1(ii,:)>6);
    R1 = squeeze(sqrt(X(1,ii,I).^2 + X(2,ii,I).^2));
    Z1 = squeeze(X(3,ii,I));
    figure(h)
    subplot(1,2,1)
    hold on;
%     plot3(R1,zeros(size(R1)),Z1,'k.','MarkerSize',4);
    plot(R1,Z1,'k.','MarkerSize',4);
    hold off;axis equal
    
    I = find(a2(ii,:)>6);
    R2 = squeeze(sqrt(X(1,ii,I).^2 + X(2,ii,I).^2));
    Z2 = squeeze(X(3,ii,I));
    subplot(1,2,2)
    hold on;
%     plot3(R2,zeros(size(R2)),Z2,'k.','MarkerSize',4);
    plot(R2,Z2,'k.','MarkerSize',4);
    hold off;axis equal
end

% np = 125;
% p = linspace(0,5*pi/4,np);
% g = figure;
% 
% for pp=1:np
%     figure(g)
%     subplot(1,2,1)
%     hold on;
%     plot3(w(:,1)*cos(p(pp)),w(:,1)*sin(p(pp)),w(:,2),'Color',[0.75,0.75,0.75],'LineWidth',0.5)
%     hold off;
%     box on
%     axis equal
%     
%     figure(g)
%     subplot(1,2,2)
%     hold on;
%     plot3(w(:,1)*cos(p(pp)),w(:,1)*sin(p(pp)),w(:,2),'Color',[0.75,0.75,0.75],'LineWidth',0.5)
%     hold off;
%     box on
%     axis equal
% end
