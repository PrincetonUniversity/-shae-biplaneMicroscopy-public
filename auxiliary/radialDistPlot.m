iiFile=1
rOutTot=[]
rTot=[]
lOut=[]
for iiFile =1:length(dataVarsOut)
centRotOut=load (char(dataVarsOut{iiFile}),'centRot');
tempTrajOut=csvread(listOut{iiFile});



centZ =27-1.5;%-6% z0(iiCell)-2; 
centRotO=centRotOut.centRot;
resCentOut = [tempTrajOut(:,1)-centRotO(1),tempTrajOut(:,2)-centRotO(2),tempTrajOut(:,3)-centZ,tempTrajOut(:,4:end)];


xRes=resCentOut(:,1)*80;
yRes=-resCentOut(:,2)*80;
zRes=resCentOut(:,3)*100/0.59;

rOut=sqrt(yRes.^2+zRes.^2);
    if abs(xRes)<700
        rOutTot=[rOutTot,rOut'];
        rTot=[rTot,rOut'];
    end
    lOut(iiFile)=length(xRes);
end
figure
histogram(lOut,'BinWidth',5)

%%
rInTot=[]
lIn =[]
for iiFile =1:length(dataVarsIn)
centRotIn=load (char(dataVarsIn{iiFile}),'centRot');
tempTrajIn=csvread(listIn{iiFile})



centZ =27-1.5%-6% z0(iiCell)-2; 
centRotO=centRotIn.centRot
resCentIn = [tempTrajIn(:,1)-centRotO(1),tempTrajIn(:,2)-centRotO(2),tempTrajIn(:,3)-centZ,tempTrajIn(:,4:end)];


xRes=resCentIn(:,1)*80;
yRes=-resCentIn(:,2)*80;
zRes=resCentIn(:,3)*100/0.59;
rIn=sqrt(yRes.^2+zRes.^2)

    if abs(xRes)<700
        rInTot=[rInTot,rIn']
        rTot=[rTot,rIn'];
    end
lIn(iiFile)=length(xRes);
end


figure
histogram(lIn,'BinWidth',5)



%%
figure
H=histogram([rInTot,rOutTot],'lineWidth',2,'FaceColor','r','FaceAlpha',0.5)
%hold on
%H=histogram(rOutTot,'lineWidth',2,'FaceColor','k','FaceAlpha',0.7)
%set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);

%%
figure
hold on
subplot(3,1,1)
histogram(rTot20,'lineWidth',2,'FaceColor','g','FaceAlpha',0.5,'Normalization','probability','BinWidth',30)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
xlim([0 650])
ylim([0 0.2])
subplot(3,1,2)
histogram(rTot40,'lineWidth',2,'FaceColor','r','FaceAlpha',0.5,'Normalization','probability','BinWidth',30)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
xlim([0 650])
ylim([0 0.2])
subplot(3,1,3)
H=histogram(rTot50,'lineWidth',2,'FaceColor','b','FaceAlpha',0.5,'Normalization','probability','BinWidth',30)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
xlim([0 650])
ylim([0 0.2])

%%
figure
subplot(3,1,1)
H=histogram(rInTot18,'lineWidth',2,'FaceColor','y','FaceAlpha',0.5,'Normalization','probability','BinWidth',30)
hold on
H=histogram(rOutTot18,'lineWidth',2,'FaceColor','y','FaceAlpha',0.7,'Normalization','probability','BinWidth',30)

set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
xlim([0 650])
ylim([0 0.2])
subplot(3,1,2)
histogram(rInTot40,'lineWidth',2,'FaceColor','r','FaceAlpha',0.5,'Normalization','probability','BinWidth',30)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
hold on
histogram(rOutTot40,'lineWidth',2,'FaceColor','r','FaceAlpha',0.7,'Normalization','probability','BinWidth',30)
xlim([0 650])
ylim([0 0.2])
subplot(3,1,3)
hold on
histogram(rInTot15,'lineWidth',2,'FaceColor','m','FaceAlpha',0.5,'Normalization','probability','BinWidth',30)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
histogram(rOutTot15,'lineWidth',2,'FaceColor','m','FaceAlpha',0.7,'Normalization','probability','BinWidth',30)

xlim([0 650])
ylim([0 0.2])
