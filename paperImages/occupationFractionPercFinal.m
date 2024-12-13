%1 is sim
%2 is exp
x1=[-18,-7,+15]
y1=[39.22,39.14,26.92]
y2=[41.19,39.98,24.57]

y1max = [39.65,39.62,27.23];
y1min = [38.82,38.68,26.64];
y2max = [42.13,40.90,25.37];
y2min = [39.91,39.06,23.97];

lweb=0.8
createOccupationFractionFigureErrorBarsTheoryCustom(x1, y1, y2, y1min, y1max,y2min, y2max,lweb,{'-3240','-1260','0','+2700'},'Particle charge (e)',[-18,-7,0,+15])

 %%

X1=[20,40,50]
Y1=[49.69,39.14,33.89]
Y2=[60.11,39.98,31.18]

y1max = [51.45,39.62,34.17];
y1min = [48.07,38.68,33.62];
y2max = [61.72,40.90,31.87];
y2min = [58.70,39.06,30.45];
createOccupationFractionFigureErrorBarsTheoryCustom(X1, Y1, Y2, y1min, y1max,y2min, y2max,lweb,{'20','30','40','50'},'Particle diameter (nm)',[20,30,40,50])
