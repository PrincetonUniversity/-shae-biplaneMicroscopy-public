%This script will create simulated data for comparison with real data and
%pick finding. 

%First we are going to start by writing a few gaussians on black
%background. Overlay noise and an offset 

%Then we want to simulate the midline of the cell. 

%We

%lets generate a cell body 
%simply put a gaussian in each pixel inside the outline of the cell. 

%next step would be to run this code through the full pipeline to generate
%the localization images from the convolution 



%Let's design a cell :
xImSz =100;
yImSz =50;

xCellCent = ceil(xImSz/2)
yCellCent = ceil(yImSz/2)
rCell = 6;
lCell = 36;
cl = lCell-2*rCell;
ytop = yCellCent + rCell;
ybot =yCellCent - rCell;
xtop=xCellCent+cl/2+rCell;
xbot=xCellCent-cl/2-rCell;



maskCell = zeros(yImSz,xImSz);
maskCell (ybot:ytop,xCellCent-cl/2:xCellCent+cl/2)= 1

circle =zeros(rCell*2+1,rCell*2+1)'

[coordsCircX,coordsCircY] = meshgrid(-rCell:rCell, -rCell:rCell);
rPos = sqrt(coordsCircX.^2+coordsCircY.^2);

for iiX = 1:size(rPos,1)
    for iiY = 1:size(rPos,2)
        if rPos(iiX,iiY)<rCell
            circle(iiX,iiY)=1;
        end
    end
end

circR = zeros(size(circle))
circL = circR;
circR(:,rCell+1:end)=circle(:,rCell+1:end)
circL(:,1:rCell+1)=circle(:,1:rCell+1)

maskCell(ybot:ytop,xbot:xCellCent-cl/2+rCell)=maskCell(ybot:ytop,xbot:xCellCent-cl/2+rCell)|circL;
maskCell(ybot:ytop,xCellCent+cl/2-rCell:xtop)=maskCell(ybot:ytop,xCellCent+cl/2-rCell:xtop)|circR;


%figure; imshow(maskCell,[])


numPhots =1000;

xCellSz =xtop-xbot+1;
yCellSz =ytop-ybot+1;


numPoints =1%10000
sigma =2;
sigma = ones(numPoints,1)*sigma;

%rng(3.5)
rng(409457)
xCent = xbot+(xtop -xbot).*rand(numPoints,1)
yCent = ybot+(ytop -ybot).*rand(numPoints,1)


iMatrix = zeros(yImSz,xImSz);
xtraj = xCellCent+Xprime;
ytraj = yCellCent+Yprime;

for iiFrame = 1:1000%:20%tsteps
    iiFrame
    iMatrix(:,:,iiFrame) = 0.2*numPhots*DVMgetMatGauss2D(xImSz,yImSz,sigma,xtraj(iiFrame),ytraj(iiFrame));%xCent(1)+Xprime,yCent(1)+Yprime);
    %figure;
    %imshow(iMatrix(:,:,iiFrame),[])
end

hold on
plot(xCent(1),yCent(1),'*r')
plot(xtraj(1000),ytraj(1000),'g*')
plot(xCellCent+Xprime,yCellCent+Yprime)
figure
plot(xCellCent+Xprime,yCellCent+Yprime)


%% Up until just designing the cell and the random walk simulation. 
for iiFrame = 1:1000

    maskCellG = imgaussfilt(maskCell,2);
    nMatrix(:,:,iiFrame) = poissrnd(numPhots,yImSz,xImSz);
    nMatrix(:,:,iiFrame) =2*maskCellG.* nMatrix(:,:,iiFrame) ;%+nMatrix(:,:,iiFrame);
%   nMatrix(:,:,iiFrame) =0.8*max(iMatrix(:))*maskCellG;%.* nMatrix(:,:,iiFrame) ;%+nMatrix(:,:,iiFrame);
    
    imMatrix(:,:,iiFrame) = iMatrix(:,:,iiFrame)+ nMatrix(:,:,iiFrame);

end

figure;imshow(imMatrix(:,:,1),[])


%%
stF =1
endF=1000

tFiltTop = temporalFiltering(imMatrix,stF,endF);
figure; imshow(tFiltTop(:,:,1),[0,max(imMatrix(:))])
hold on
plot(xCent,yCent,'g*')
figure; imshow(imMatrix(:,:,1),[0,max(imMatrix(:))])
hold on
plot(xCent,yCent,'g*')

%%    
for iiPlane = 1:size(psf,3)
    psfG(:,:,iiPlane) = imgaussfilt(psf(:,:,iiPlane),1);
end

%%
for iiFrame =1:1000
    iiFrame
    cMatrixFull(:,:,:,iiFrame) = convn(imMatrix(:,:,iiFrame),psfG);
end
%%

clear peaksOut6
clear peaksOut3
clear roughPeaks2
clear error
for iiFrame = 1:size(imMatrix,3)
   try
iiFrame
    if iiFrame ==1
        cMatrix = imMatrix(:,:,iiFrame);
        imageIn = cMatrix(:,:);
        a0(1) = 25; %3um
        a0(2) = 6; %1um
        a0(3) = size(imageIn,2)/2;
        a0(4) = size(imageIn,1)/2;
        a0(5) = 0.7;
        a0(6) = 0;
 %   fittedCell = cellFit2Stage(cMatrix(:,:,26), a0);


[fittedCell,th,x0,y0,cellMask1] = cellFit2StageNoConvRot(imageIn,a0);


       % fittedCell = cellFit2Stage(cMatrix(:,:,26), roughPeaks2);
    end

    %figure; imshow(imMatrix(:,:,iiFrame),[])
    cMatrix =  imMatrix(:,:,iiFrame);
    [roughPeaks2,intP2] = pkfinderDVM(cMatrix(:,:),max(cMatrix(:))*0.85,15) ;% be generous in locating peaks, the ones that are just noise will be removed later
   
    resMat = cMatrix -max(imageIn(:))*fittedCell; 
     [roughPeaks3,intP3] = pkfinderDVM(resMat(:,:),max(resMat(:))*0.85,15) ;% be generous in locating peaks, the ones that are just noise will be removed later

    %figure;imshow(cellMat,[])

     figure;imshow(resMat,[])
    peaksOut6(:,:,iiFrame) = centfinderDVM (resMat, roughPeaks3(1,:),9);
    peaksOut7(:,:,iiFrame) = centfinderDVM (resMat, [xtraj(iiFrame),ytraj(iiFrame)],9);
    peaksOut8(:,:,iiFrame) = centfinderDVM (cMatrix, roughPeaks2(1,:),9);
    peaksOut9(:,:,iiFrame) = centfinderDVM (cMatrix, [xtraj(iiFrame),ytraj(iiFrame)],9);
    peaksOut3 = centfinderDVM (cMatrix(:,:), roughPeaks2,7);
   peaksTemp =peaksOut6(:,:,iiFrame);
   peaksTemp2 =peaksOut8(:,:,iiFrame);
      figure ; imshow(cMatrix(:,:),[])
      hold on
      coralColor =1/255*[255,138,138];
      grayColor =1/255*[115,115,115];
%      plot(peaksOut3 (:,1), peaksOut3(:,2),'r*')
      plot(peaksTemp2 (:,1), peaksTemp2(:,2),'*','Color', grayColor)
      plot(peaksTemp (:,1), peaksTemp(:,2),'*','Color', coralColor)
      plot(xtraj(iiFrame),ytraj(iiFrame),'*g')
     hold on
     plot(roughPeaks2 (:,1), roughPeaks2(:,2),'m*')

   catch
         error(iiFrame)=1
         warning('error')
     end
    
end

%%
%peaksTraj = squeeze(peaksOut6)';

for iiFrame =1:size(peaksOut6,1)
    temp = peaksTraj(iiFrame,:);
    if (sum(temp ==[0 , 0])>1)
        5
        peaksTraj(iiFrame,:)=[NaN,NaN]
    end 
end 
%%
getRidOfZerosSimData(peaksTraj(:,1))
%%
peaksTraj = squeeze(peaksOut6)';
errorX =xtraj(1:1000) - peaksTraj((1:1000),1);
errorY =ytraj(1:1000) - peaksTraj((1:1000),2);

errorXc=errorX;
errorYc=errorY;


        errorXc(error==1)=[];
        %errorXc(abs(errorXc)>5)=[]; 
        errorYc(error==1)=[];
        %errorYc(abs(errorYc)>5)=[];


peaksTraj2 = squeeze(peaksOut8)';
errorX2 =xtraj(1:1000) - peaksTraj2((1:1000),1);
errorY2 =ytraj(1:1000) - peaksTraj2((1:1000),2);

errorXc2=errorX2;
errorYc2=errorY2;

        errorXc2(error==1)=[];
      %  errorXc(abs(errorXc)>5)=[]; 
        errorYc2(error==1)=[];
       % errorYc(abs(errorYc)>5)=[];

figure


subplot(1,2,1)
histogram(errorXc*80,[-1000:50:1000],'lineWidth',2,'FaceColor','r','FaceAlpha',0.5)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
axis square
%xlim([0 8])
ylim([0 400])
title('x')
hold on
histogram(errorXc2*80,[-1000:50:1000],'lineWidth',2,'FaceColor','k','FaceAlpha',0.5)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
subplot(1,2,2)
histogram(errorYc*80,[-1000:50:1000],'lineWidth',2,'FaceColor','r','FaceAlpha',0.5)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
hold on
histogram(errorYc2*80,[-1000:50:1000],'lineWidth',2,'FaceColor','k','FaceAlpha',0.5)
set(findobj(gcf,'type','axes'),'FontSize',25,'FontWeight','Bold', 'LineWidth', 3);
title('y')
axis square
%xlim([0 8])
ylim([0 400])


mean(errorXc*80)
mean(errorYc*80)
std(errorXc*80)
std(errorYc*80)

std(errorXc*80)/sqrt(size(errorXc,1))
std(errorYc*80)/sqrt(size(errorXc,1))


mean(errorXc2*80)
mean(errorYc2*80)
std(errorXc2*80)
std(errorYc2*80)

std(errorXc2*80)/sqrt(size(errorXc2,1))
std(errorYc2*80)/sqrt(size(errorXc2,1))



%%
figure
plot(peaksTraj(:,1),peaksTraj(:,2),'*')
hold on
plot(xtraj,ytraj,'g*')

       