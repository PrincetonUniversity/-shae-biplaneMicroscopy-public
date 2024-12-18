function createOccupationFractionFigureErrorBarsTheoryCustom(X1, Y1, Y2, y1min, y1max,y2min, y2max,lw,labelsIn,xTit,xTicksIn)
%CREATEFIGURE(X1, Y1, Y2, YNEG1, YPOS1)
%  X1:  vector of x data
%  Y1:  vector of y data
%  Y2:  errorbar y vector data
%  YNEG1:  errorbar ynegativedelta vector data
%  YPOS1:  errorbar ypositivedelta vector data

%  Auto-generated by MATLAB on 18-Dec-2022 21:52:10. Edited by DVM for
%  function. 




% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

[d0,d1,d2]=errorBarCustom(X1,Y2,y2min,y2max,lw)
set(d0,'MarkerSize',20,'MarkerFaceColor','none',...
    'Marker','.','LineWidth',2.5,...
    'LineStyle','none',...
    'DisplayName','sim',...
    'Color',[0 0 0])
set(d1,'MarkerSize',20,'MarkerFaceColor','none',...
    'Marker','none','LineWidth',2.5,...
    'LineStyle','-',...
    'DisplayName','sim',...
    'Color',[0 0 0])

set(d2,'MarkerSize',20,'MarkerFaceColor','none',...
    'Marker','none','LineWidth',2.5,...
    'LineStyle','-',...
    'DisplayName','sim',...
    'Color',[0 0 0])



[d0,d1,d2]=errorBarCustom(X1,Y1,y1min,y1max,lw)
set(d0,'MarkerSize',20,'MarkerFaceColor','none',...
    'Marker','.','LineWidth',2.5,...
    'LineStyle','none',...
    'DisplayName','sim',...
    'Color',[1 0 0])
set(d1,'MarkerSize',20,'MarkerFaceColor','none',...
    'Marker','none','LineWidth',2.5,...
    'LineStyle','-',...
    'DisplayName','sim',...
    'Color',[1 0 0])

set(d2,'MarkerSize',20,'MarkerFaceColor','none',...
    'Marker','none','LineWidth',2.5,...
    'LineStyle','-',...
    'DisplayName','sim',...
    'Color',[1 0 0])
set(axes1,'FontName','Cambria','FontSize',25,'LineWidth',3,'XTick',xTicksIn,'YTick',[0.2:0.1:0.7])%,...


% %



% Create xlabel
xlabel(xTit);

% Create title
title({''});

% Create ylabel
ylabel({'% Nucleoid occupation time'});

xlim(axes1,[X1(1)-5, X1(3)+5])
ylim(axes1,[20 70]);

box(axes1,'on');
% Set the remaining axes properties
% set(axes1,'FontName','Cambria','FontSize',25,'FontWeight','bold','LineWidth',3,'XTick',X1,...
%    'XTickLabel',labelsIn);

set(axes1,'FontName','Cambria','FontSize',25,'LineWidth',3,'XTick',xTicksIn,'YTick',20:10:70,...
    'XTickLabel',labelsIn);

x0=10;
y0=10;
width=400;
height=400;
axis square
%set(gcf,'units','points','position',[x0,y0,width,height])
