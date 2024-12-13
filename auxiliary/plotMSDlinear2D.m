function [p]=plotMSDlinear(strucIn,sizeT,color)

%'.k-'
colorX = [55,126,184]/255; %yz
colorY =[77,175,74]/255;%xz
colorZ =[228,26,28]/255;%xy

%sizeT =3



kCycle = strucIn.kCycle;
rowMeanFull=strucIn.rowMeanX+strucIn.rowMeanY;
intervalsF=strucIn.intervalsF;
% set(gca, 'FontSize', 20)%,'BinWidth', 50)
% xlabel('\tau (s)')
% ylabel('MSD (\mum^2)')
% title('50 nm MSD')
% 
%  plot(log(kCycle*(1:sizeT)),log(rowMeanFull(1:sizeT)*10^-6),'ok','MarkerSize',6,'LineWidth',2,'MarkerEdgeColor','b')
%  lsline

%figure
%ylim([0,0.06]);

intervalsX=strucIn.intervalsX;
intervalsY=strucIn.intervalsY;


neg= rowMeanFull(1:sizeT)*10^-6 -intervalsX(1:sizeT,1)*10^-6 - intervalsY(1:sizeT,1)*10^-6 ;%lower
pos=-( rowMeanFull(1:sizeT)*10^-6 -intervalsX(1:sizeT,2)*10^-6 - intervalsY(1:sizeT,2)*10^-6 );%upper


errorbar((1:sizeT)*kCycle,rowMeanFull(1:sizeT)*10^-6,neg,pos,...
    color,'LineWidth',2)
%title('MSD')
set(gca, 'FontSize', 20)%,'BinWidth', 50)
xlabel('\tau (s)')
ylabel('MSD (\mum^2)')
%title('50 nm MSD')
set(gca,'xscale')
set(gca,'yscale')
%plot ((1:sizeT)*kCycle,(1:sizeT)*kCycle,'g--')


x =  (kCycle*(1:sizeT))';
y =  (rowMeanFull(1:sizeT)*10^-6);
set(gca, 'FontSize', 20)%,'BinWidth', 50)
xlabel('\tau (s)')
ylabel('MSD (\mum^2)')
p = polyfit(x,y,1)
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal
title (strcat ('MSD 40nm = Dt^{\alpha}, \alpha= ',num2str(p(1))) )
p(2)

set(gca,'FontSize',24)
set(gca,'LineWidth',2.5)