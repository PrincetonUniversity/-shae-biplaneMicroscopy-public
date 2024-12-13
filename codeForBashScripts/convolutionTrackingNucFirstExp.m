function [] = convolutionTrackingNucFirstExp(imageStructPath,fileStructPath,jFile,iiCell) 

%%
%First lets start by defining the basic PSF: 
%PSF = zeros() 


%% constants

BOX_SIZE_PEAK = 9;
CAMERA_DIVIDER_PTS = [[1,1,512,512,1]',[1,270,270,1,1]'];
IMAGE_PIXEL_SIZE = 512;
PEDESTAL = 99;


%% Here we are getting the PSFs 
%parentDir = uigetdir('','Please choose the parent directory');

%[psfTop,psfBot,psfGTop,psfGBot]= getAggregatedPSF (parentDir);
%save('calibrations20210622.mat','-append','psfTop','psfBot','psfGTop','psfGBot')
% psfFull = psfTop+psfBot;
% 
% 
% 
% sizeA = size(psfFull);
% 
% Afull = reshape(psfFull(floor(sizeA(1)/2),:,:),[sizeA(1),sizeA(3)]);
% figure
% 
% imshow(Afull,[])
%%
%We are going to load everything we need, in the real code this paths will
%come from the cluster but it's good to think about the loading since now. 
% %%
%[imageStructFilename,imageStructPath] = uigetfile();
%load(fullfile(imageStructPath,imageStructFilename));

%[fileStructFilename,fileStructPath] = uigetfile();
%load(fullfile(fileStructPath,fileStructFilename));

%[calFilename,calPath] = uigetfile();
% load(fullfile(calPath,calFilename));
load('calibrations20210622.mat');
load('PSFcalibrations20210622.mat');
load(fileStructPath)
load(imageStructPath)
%%

%%
%jFile =1;% will come from the cluster . 
tiffName = char(fileStruct.imagePath(jFile));
load(char(fileStruct.variableList(jFile)),'savedTop');
load(char(fileStruct.variableList(jFile)),'savedBot');
%%
limits = [[1,1,512,512,1]',[1,270,270,1,1]']; %This limits should be adjusted for the new images. 
imageSize =IMAGE_PIXEL_SIZE;
[pathMD,nameMD,ext] = fileparts(tiffName);
%%
clear nucleoidList24
clear nucleoidList25
for iiNucImages = 1:size(fileStruct.nucleoidParent') 
    imageDirPath =char(fileStruct.nucleoidParent(iiNucImages));
    matFiles = dir(fullfile(imageDirPath,'*.tif'));
    nucleoidList24{iiNucImages}=fullfile(matFiles(24).folder,matFiles(24).name);
    nucleoidList25{iiNucImages}=fullfile(matFiles(25).folder,matFiles(25).name);
    nucleoidList21{iiNucImages}=fullfile(matFiles(21).folder,matFiles(21).name);
    
end 


% clear membraneList18
% clear membraneList19
% for iiMemImages = 1:size(fileStruct.membraneParent') 
%     imageDirPath =char(fileStruct.membraneParent(iiMemImages));
%     matFiles = dir(fullfile(imageDirPath,'*.tif'));
%     membraneList18{iiMemImages}=fullfile(matFiles(24).folder,matFiles(24).name);
%     membraneList19{iiMemImages}=fullfile(matFiles(25).folder,matFiles(25).name);
% %     membraneList18{iiMemImages}=fullfile(matFiles(18).folder,matFiles(18).name);
% %     membraneList19{iiMemImages}=fullfile(matFiles(19).folder,matFiles(19).name);
% %     for i=1:size(matFiles)
% %         Imem = imread(fullfile(matFiles(i).folder,matFiles(i).name));
% %         figure; imshow(Imem,[])
% %     end
% %     
% end

close all
clear membraneList24
clear membraneList25
clear membraneList26
for iiMemImages = 1:size(fileStruct.membraneParent') 
    imageDirPath =char(fileStruct.membraneParent(iiMemImages));
    matFiles = dir(fullfile(imageDirPath,'*.tif'));
 
    membraneList23{iiMemImages}=fullfile(matFiles(23).folder,matFiles(24).name);
    membraneList24{iiMemImages}=fullfile(matFiles(24).folder,matFiles(25).name);
    membraneList25{iiMemImages}=fullfile(matFiles(25).folder,matFiles(26).name);
 
    
end 




%%

Ibright=imread(fileStruct.transList{jFile});
figure;imshow(Ibright,[])


%% Preparing some file names, translating membrane image, etc. 
%I=savedComp(limits(2,2)+1:imageSize ,limits(1,1):limits(3,1),1);
%Ifull=savedComp(limits(2,2)+1:imageSize ,limits(1,1):limits(3,1),:);
Inuc = imread(char(nucleoidList24(jFile)))+imread(char(nucleoidList25(jFile)));
Inuc = Inuc -2*PEDESTAL;
Imem = imread(char(membraneList23(jFile)))+imread(char(membraneList24(jFile)));
Imem = Imem-2*PEDESTAL;
Imem = imtranslate(Imem,[4, 7]); %empirically seen translation to match displacement of blue channel. 

Imem23 = imread(char(membraneList23(jFile)))-PEDESTAL;
Imem24 = imread(char(membraneList24(jFile)))-PEDESTAL;
Imem25 = imread(char(membraneList25(jFile)))-PEDESTAL;


Imem23 = imtranslate(Imem23,[4, 7]); %empirically seen translation to match displacement of blue channel. 
Imem24 = imtranslate(Imem24,[4, 7]); %empirically seen translation to match displacement of blue channel. 
Imem25 = imtranslate(Imem25,[4, 7]); %empirically seen translation to match displacement of blue channel. 



%imshow(Imem,[])
for iiRow=1:7
    Imem(iiRow,:) = Imem(8,:);
    Imem23(iiRow,:) = Imem23(8,:);
    Imem24(iiRow,:) = Imem24(8,:);
    Imem25(iiRow,:) = Imem25(8,:);
end
for iiCol =1:4
    Imem(:,iiCol) = Imem(:,5);
    Imem23(:,iiCol) = Imem23(:,5);
    Imem24(:,iiCol) = Imem24(:,5);
    Imem25(:,iiCol) = Imem25(:,5);
end
%Inuc = imageStruct(jFile).imNuc; % Use when saved image instead of path

fileNameTraj = eraseBetween(char(tiffName),1,'2');
[filepath,name,ext]=fileparts(fileNameTraj);
mkdir(filepath);
rotTrajFolder = strcat (filepath);
mkdir(rotTrajFolder);
saveFileNameTraj = fullfile(rotTrajFolder,'traj_cell_');
saveFileNameRot = fullfile(rotTrajFolder,'rotated_cell_');
saveFileNameCropped = fullfile(rotTrajFolder,'cropped_cell_');
saveFileNameMeasurementsP = fullfile(rotTrajFolder,'measurementsP_cell_');
saveFileGrainRot = strcat(saveFileNameRot,num2str(iiCell));
saveFileNameGrainCropped = strcat(saveFileNameCropped,num2str(iiCell));

C = imfuse(Inuc,Imem,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
figure;imshow(C,[])

%% Load the cells drawn by user. 

clear a;

a=imageStruct(jFile).cells;

savedTop= squeeze(savedTop);
savedBot= squeeze(savedBot);

[maskCell,xBot,yBot] = makeMaskCell (a,iiCell,2,2);

figure;
imshow(double(maskCell ).* double(Inuc), [])

figure;
imshow(double(maskCell ).* double(Imem), [])
%%
boxBot = [xBot,yBot];
boxTop = fnval(stprime,boxBot')';

rndTop = round(boxTop);
rndBot = round (boxBot);

maxBot=max(rndBot)+10;
minBot=min(rndBot)-10;

fileNameTrajMask = char(tiffName);
[filepathMask,~,~]=fileparts(fileNameTrajMask);
mkdir(filepathMask);
rotTrajFolderMask = strcat (filepathMask);
mkdir(rotTrajFolderMask); 
pathCells =fullfile(rotTrajFolderMask,num2str(iiCell))
pathNetCells = fullfile(rotTrajFolderMask,num2str(iiCell),strcat('mask','.mat'))
netMask = load(pathNetCells);

figure
subplot(1,3,1)
ImNetMaskIn= squeeze(netMask.pred(1,1,:,:));
imshow(ImNetMaskIn,[])
subplot(1,3,2)
ImNetMaskB= squeeze(netMask.pred(1,2,:,:));
imshow(1-ImNetMaskB,[])
subplot(1,3,3)
ImNetMask = ImNetMaskIn.*(1-ImNetMaskB);
imshow(ImNetMask)
fullNetMask = zeros(512,512)
fullNetMask(minBot(2):maxBot(2),minBot(1):maxBot(1))=ImNetMask;
figure; imshow(fullNetMask,[])
%%

%Finding original box orientation from user selected cells
a=imageStruct(jFile).cells;
x= cat(1,a(1:2,1,iiCell),a(3:4,1,iiCell));
y= cat(1,a(1:2,2,iiCell),a(3:4,2,iiCell));
bw=poly2mask(x,y,512,512);
%figure; imshow(bw,[])
se90 = strel('line',2,90);
se0 = strel('line',2,0);
bw = imdilate(bw,[se90 se0]); %close edges
figure; imshow(bw,[])


I2 = uint8(255 * mat2gray(fullNetMask));
figure; imshow(double(I2).*double(bw),[])

[~,thresh] =edge(I2, 'log');
edges =edge(I2, 'log',0.6*thresh); %seems to work best for nucleoids
edges = bw.*edges;
figure; imshow(edges,[]) 

se90 = strel('line',1,90);
se0 = strel('line',1,0);
BWsdil4 = imdilate(edges,[se90 se0]); %close edges
figure; imshow(BWsdil4,[])

BWsdil2 =filledgegaps(BWsdil4,3);
figure; imshow(BWsdil2,[])
BWdfill = imfill(BWsdil2,'holes'); %fill contours
figure; imshow(BWdfill,[])

masked =BWdfill;%bw.*BWdfill;
   
bw = bwareaopen(masked, 15);     %Removes noise 
figure;imshow(bw,[])


   
measurements =regionprops(bw,'Image','Perimeter','Area','centroid','minorAxisLength','majorAxisLength','Extrema','Orientation','BoundingBox');
area =[ measurements(:).Area]
[~,idx]=max(area)
measurements = measurements(idx)

figure; imshow(bw,[])
netCell=zeros(512,512);
netCell(round(measurements.BoundingBox(2)):round(measurements.BoundingBox(2))+round(measurements.BoundingBox(4))-1,round(measurements.BoundingBox(1)):round(measurements.BoundingBox(1))+round(measurements.BoundingBox(3)-1))=measurements.Image;


figure; imshow(netCell,[])
hold on
plot(measurements.Centroid(1),measurements.Centroid(2),'b*')

se90 = strel('line',4,90);
se0 = strel('line',4,0);
bw = imdilate(netCell,[se90 se0]); %close edges
%figure; imshow(bw,[])


hold on
plot(measurements.Centroid(1),measurements.Centroid(2),'b*')


 BWoutline = bwperim(bw);
    Segout = Imem; 
    Segout(BWoutline) = 255; 
   figure, imshow(Segout,[],'initialMagnification','fit'), title('outlined original image');
    hold on
    
measurements =regionprops(bw,'Image','Perimeter','Area','centroid','minorAxisLength','majorAxisLength','Extrema','Orientation','BoundingBox');


plot(measurements.Centroid(1),measurements.Centroid(2),'b*')




%%
figure
imshow(Imem,[])
%%
boxBot = [xBot,yBot];
boxTop = fnval(stprime,boxBot')';

rndTop = round(boxTop);
rndBot = round (boxBot);



% maxBot=max(rndBot);
% minBot=min(rndBot);
% 
% figure;imshow(double(maskCell).*double(Inuc),[])
% maskCrop= maskCell(minBot(2):maxBot(2),minBot(1):maxBot(1));
% imageIn = double(bw(minBot(2):maxBot(2),minBot(1):maxBot(1)));
% 
% a0(1) = round(measurements.MajorAxisLength)+2; %3um
% a0(2) = round(measurements.MinorAxisLength/2)+1; %1um
% a0(3) = measurements.Centroid(1)-minBot(1)+1;%remember matlab counts from 1 :'(
% a0(4) = measurements.Centroid(2)-minBot(2)+1;
% a0(5) = 1;
% a0(6) = measurements.Orientation;
% 
% 
%    figure; imshow(imageIn,[0,max(imageIn(:))])
%     title('Imagein')
%     hold on
%     plot(a0(3),a0(4),'r*')
% 
%     figure; imshow(fitCellonlyNoConvRot(a0,imageIn),[])
%     figure; imshow(fitCellonlyNoConvRot(a0,imageIn)-imageIn,[0,1])
%     title('Initial guess')
% cellParams = fitCellScriptNoConvRot (a0,imageIn) ;
% fitCell =fitCellonlyNoConvRot(cellParams,imageIn);
% figure;imshow(fitCell-imageIn,[0,1])
% 
% [fittedCellNoConv,th,x0,y0,cellMask1] = cellFit2StageNoConvRot(imageIn,a0);
% 
% cellMask=imbinarize(cellMask1);
% figure; imshow(cellMask)



maxBot=max(rndBot)+10;
minBot=min(rndBot)-10;







%%
clear cropCtopImFull
clear cropCbotImFull


for iiFrame =1:size(savedTop,3)
    cropCtopIm = convn(savedTop(minBot(2):maxBot(2),minBot(1):maxBot(1),iiFrame),psfGTop);
    cropCbotIm = convn(savedBot(minBot(2):maxBot(2),minBot(1):maxBot(1),iiFrame),psfGBot);
  
    cropCtopImFull(:,:,:,iiFrame)= cropCtopIm(11:end-10,11:end-10,:);
    cropCbotImFull(:,:,:,iiFrame) = cropCbotIm(11:end-10,11:end-10,:);
end


%%

cropImTop= savedTop(minBot(2):maxBot(2),minBot(1):maxBot(1),:);
cropImBot= savedBot(minBot(2):maxBot(2),minBot(1):maxBot(1),:);
 

% %Let's try normalizing the images so that both are in a 0 to 1 scale. 
% normImTop = croppedImTop/(max(croppedImTop(:)));
% normImBot = croppedImBot/(max(croppedImBot(:)));
% figure
% subplot(1,2,1) ; imshow(normImTop(:,:,1)+ normImBot(:,:,1),[])
% subplot(1,2,2) ;imshow((croppedImTop(:,:,1)+croppedImBot(:,:,1)),[])

imageIn =bw(minBot(2):maxBot(2),minBot(1):maxBot(1));%cropImTop(:,:,1)+ cropImBot(:,:,1);
a0(2) = round(measurements.MinorAxisLength/2); %1um
a0(1) = round(measurements.MajorAxisLength)-18; %3um
a0(2) = round(measurements.MinorAxisLength/2); %1um
a0(3) = measurements.Centroid(1)-minBot(1);
a0(4) = measurements.Centroid(2)-minBot(2);
a0(5) = 0.7;
a0(6) = measurements.Orientation;

[fittedCellNoConv,th,x0,y0,cellMask1] = cellFit2StageNoConvRot(imageIn,a0);




% %imageIn =cropImTop(10:end-10,10:end-10,1)+ cropImBot(10:end-10,10:end-10,1);
% 
% a0(1) = round(measurements.MajorAxisLength)-round(measurements.MinorAxisLength); %3um
% a0(2) = round(measurements.MinorAxisLength/2); %1um
% a0(3) = x0;
% a0(4) = y0;
% a0(5) = 0.7;
% a0(6) = 90; %%added this from when cellFit2StageNoConvRot commented. 
% a0(7) = th;%before a0(6) = measurements.Orientation;
% 
% 
% %[fittedCellNoConv,th,x0,y0,cellMask1] = cellFit2StageNoConvRot(imageIn,a0);
% 
% [fittedCellNoConv,cellParams,cellMask1] = copycellFit2StageRotBen(imageIn,a0);
% 
% th = cellParams(7);
% y0 =cellParams(4);
% x0 =cellParams(3);



imageIn =cropImTop(:,:,1)+ cropImBot(:,:,1);
cellMask=imbinarize(cellMask1);
resMatNoConv = imageIn -0.25*max(imageIn(:))*fittedCellNoConv; 
figure;imshow(resMatNoConv,[])


figure; imshow(cellMask,[])
figure; imshow(imageIn,[])
figure; imshow(double(Imem(minBot(2):maxBot(2),minBot(1):maxBot(1),:)),[])%.*double(cellMask),[])
 BWoutline = bwperim(cellMask);
    Segout = Inuc(minBot(2):maxBot(2),minBot(1):maxBot(1),:); 
    Segout(BWoutline) = 512; 
   figure, imshow(Segout,[],'initialMagnification','fit'), title('outlined original image');
 
C = imfuse(Inuc(minBot(2):maxBot(2),minBot(1):maxBot(1),:),Imem(minBot(2):maxBot(2),minBot(1):maxBot(1),:),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
figure; imshow(C)
    Segout = C; 
    Segout(BWoutline) = 255; 
   figure, imshow(Segout,[],'initialMagnification','fit'), title('outlined original image');
 

figure
subplot(1,2,1)
imshow(imageIn.*cellMask,[0,max(max(imageIn.*cellMask))])
title('before fitting')
subplot(1,2,2)
imshow(resMatNoConv.*cellMask,[0,max(max(imageIn.*cellMask))])
title('after subtracting fitting')


scale = max(max(max((cropImTop+cropImBot))));
figure
subplot(1,2,1)
imshow(cropImTop(:,:,1),[0,scale])
title('top plane')
subplot(1,2,2)
imshow(cropImBot(:,:,1),[0,scale])
title('bottom plane')

figure
subplot(1,2,1)
imshow(cropImTop(:,:,1)+cropImBot(:,:,1),[0,scale])
title('raw added 2 channels')
subplot(1,2,2)
imshow(resMatNoConv.*cellMask,[0,scale])
title('post processing')


%%
clear peaksOutNoConv
clear maskedCellFull
clear maskedCellTop
clear maskedCellBot 
clear resPeaksFull

indNC =1;
5
figure
for iiFrame = 1:1000
    try
    iiFrame
    
    tempTop = cropImTop(:,:,iiFrame);
    tempBot = cropImBot(:,:,iiFrame);
    imageInNoConv= tempTop+tempBot;

    resMatNoConv = imageInNoConv -0.6*max(imageInNoConv(:))*fittedCellNoConv; 

    maskedCellNoConv =resMatNoConv.*cellMask;
    maskedCellFull(:,:,iiFrame)=maskedCellNoConv;
    %figure;imshow(maskedCellFull(:,:,iiFrame),[])
    maskedCellTop(:,:,iiFrame)=tempTop.*cellMask;
    maskedCellBot(:,:,iiFrame)=tempBot.*cellMask;
    scale= max(maskedCellNoConv(:));
    
  %subplot(1,3,1);imshow(tempTop.*cellMask,[0,scale])
  %subplot(1,3,2);imshow(tempBot.*cellMask,[0,scale])
  %subplot(1,3,3);imshow(resMatNoConv.*cellMask,[0,scale])
  
    imshow(resMatNoConv.*cellMask,[0,scale])
    [roughPeaks2,~] = pkfinderDVM(maskedCellNoConv,max(maskedCellNoConv(:))*0.85,15) ;% be generous in locating peaks, the ones that are just noise will be removed later
    roughPeaks2
    hold on
    plot(roughPeaks2(:,1),roughPeaks2(:,2),'.r')
    
    [tempPeaksNC,resPeaks]= centfinderDVM (maskedCellNoConv, roughPeaks2,9);
    
    
    idxPeaks =(tempPeaksNC(:,1)==0&tempPeaksNC(:,2)==0)
    

    tempPeaksNC(idxPeaks,:)=[]

    resPeaks(idxPeaks)=[]
   
    timeVecNC =zeros(size(tempPeaksNC,1),1)+iiFrame;
    lengthTNC =size(tempPeaksNC,1)
    peaksOutNoConv(indNC:indNC+lengthTNC-1,:)=[tempPeaksNC,timeVecNC];
    peaksTemp = peaksOutNoConv(peaksOutNoConv(:,3)==iiFrame,1:2);

    plot(tempPeaksNC(:,1),tempPeaksNC(:,2),'pb')
    drawnow
    resPeaksFull(indNC:indNC+lengthTNC-1,:)=resPeaks
    %  hold on
     % plot(peaksTemp(:,1),peaksTemp(:,2),'r*')
%      plot(roughPeaks2(:,1),roughPeaks2(:,2),'b*')
    


indNC = indNC + lengthTNC;
%indC = indC + lengthTC;
    catch

    end
end
6
%%
% x0
% y0
% 
% sizex = size(imageInNoConv, 2);
% sizey = size(imageInNoConv, 1);
% 
% isOddx = mod(sizex,2);
% isOddy = mod(sizey,2);
% 
% if isOddx
%     newSzX= sizex+1;
% else
%     newSzX =sizex;
% end
% if isOddy 
%     newSzY = sizey+1;
% else
%     newSzY = sizey;
% end
% 
% oddSzdImage = zeros(SzY,SzX)
% figure = 

%%
Irot = imrotate(maskedCellFull,-th);
I = imageInNoConv;
peaksOutZero= [peaksOutNoConv(:,1)-x0,peaksOutNoConv(:,2)-y0,peaksOutNoConv(:,3)];
clear posXYrot
    resXYtrans =[peaksOutNoConv(:,1),peaksOutNoConv(:,2)];
    % plot(resXYtrans(:,1),resXYtrans(:,2), 'w.')

    %Put the origin of coordinates in center of image
    resXYtrans2 = [resXYtrans(:,1)-size(I,2)/2 , resXYtrans(:,2)-size(I,1)/2];
    %plot(resXYtrans2(:,1),resXYtrans2(:,2), 'b.')
  
    rotMatrix = [cosd(-th), -sind(-th); sind(-th), cosd(-th)];
    for ii=1:size(peaksOutZero,1)
        posXYrot(ii,:) =resXYtrans2(ii,:) * rotMatrix;
    end
%figure;plot(posXYrot(:,1),posXYrot(:,2),'r*')


    posXYrot = [posXYrot(:,1)+size(I,2)/2 , posXYrot(:,2)+size(I,1)/2];
    posXYrot = [posXYrot(:,1)+ (size(Irot,2)-size(I,2))/2 , posXYrot(:,2)+ (size(Irot,1)-size(I,1))/2];
   posXYrot = [posXYrot(:,1),posXYrot(:,2),peaksOutNoConv(:,3)];

   
    for iiFrame =1:3%:size(maskedCellFull,3)
        
        figure; imshow(Irot(:,:,iiFrame),[])
        hold on
        posTemp = posXYrot(posXYrot(:,3)==iiFrame,:);
         plot(posTemp(:,1),posTemp(:,2),'b.')
         figure ; imshow(maskedCellFull(:,:,iiFrame))
         hold on 
         posTemp = peaksOutNoConv(peaksOutNoConv(:,3)==iiFrame,:);
         plot(posTemp(:,1),posTemp(:,2),'b.')

    end 
    %plot(posXYrot(1,1),posXYrot(1,2),'rp')
    %posRotatedInNM=[posXYrot(:,1)*80, posXYrot(:,2)*80, posGrain(:,3:8)];
    %posRotated=[posXYrot(:,1), posXYrot(:,2), posGrain(:,3:8)];

centZero = [x0-size(I,2)/2,y0-size(I,1)/2];
centRot =centZero * rotMatrix;
centRot = [centRot(:,1)+size(I,2)/2 , centRot(:,2)+size(I,1)/2];
centRot = [centRot(:,1)+ (size(Irot,2)-size(I,2))/2 , centRot(:,2)+ (size(Irot,1)-size(I,1))/2];
   
plot(centRot(1),centRot(2),'r*')

%%
%%
%First let's reorient the cells 
a0(1) = 30; %3um
a0(2) = 6; %1um
a0(3) = size(cellMask,2)/2;
a0(5) = 0.7;
a0(6) = th;
[fittedCellRot,params] = cellFit2StageRot(cropCbotImFull(:,:,23,1).* cellMask,a0);

rotCcell =imrotate(fittedCellRot,-th,'bicubic');
figure;imshow(rotCcell,[])

rotCbotPreMasked =imrotate(cropCbotImFull(:,:,:,:),-th,'bicubic');
rotCtopPreMasked =imrotate(cropCtopImFull(:,:,:,:),-th,'bicubic');
rotCMask =imrotate(cellMask,-th,'bicubic');
figure;imshow(rotCMask,[])
   
se90 = strel('line',1,90);
se0 = strel('line',2,0);
rotCMask = imdilate(rotCMask,[se90 se0]); %close edges

rotCtop= rotCtopPreMasked .* rotCMask;
rotCbot= rotCbotPreMasked .* rotCMask;
figure;imshow(rotCbot(:,:,23,1),[]);


figure; imshow(rotCbot(:,:,23,1),[]);
figure; imshow(rotCtop(:,:,29,1),[]);

maxCbot = max(rotCbot(:));
maxCtop = max(rotCtop(:));


% %First let's reorient the cells 
% cropCtopImFull;
% a0(1) = 30; %3um
% a0(2) = 6; %1um
% a0(3) = size(cellMask,2)/2;
% a0(4) = size(cellMask,1)/2;
% a0(5) = 0.7;
% a0(6) = th;
% [fittedCellRot,params] = cellFit2StageRot(cropCbotImFull(:,:,23,1).* cellMask,a0);
% 
% rotCcell =imrotate(fittedCellRot,-th,'bicubic');
% figure;imshow(rotCcell,[])
% 
% rotCbotPreMasked =imrotate(cropCbotImFull(:,:,:,:),-th,'bicubic');
% rotCMask =imrotate(cellMask,-th,'bicubic');
% rotCbot= rotCbotPreMasked .* rotCMask;
% figure;imshow(rotCbot(:,:,23,1),[]);
% rotCtop =imrotate(cropCtopImFull(:,:,:,:).* cellMask,-th,'bicubic');
% 
% figure; imshow(rotCbot(:,:,23,1),[]);
% figure; imshow(rotCtop(:,:,29,1),[]);
% 
% maxCbot = max(rotCbot(:));
% maxCtop = max(rotCtop(:));
%%
% a0(1) = 30; %3um
% a0(2) = 6; %1um
% a0(3) = size(rotCbot,2)/2;
% a0(4) = size(rotCbot,1)/2;
% a0(5) = 0.7;
% a0(6) = 0;

%[fittedCellRot,params] = cellFit2StageRot(rotCbot(:,:,23,1),a0);
 %%
 intFTopB=[];
 intFBotB=[];
 intFTopAll=[];
 intFBotAll=[];
 tIndAll =[];
 7
     for iiFrame =1:size(maskedCellFull,3)
         iiFrame
         imgTemp =maskedCellFull(:,:,iiFrame);
       
         posTemp = peaksOutNoConv(peaksOutNoConv(:,3)==iiFrame,:);
%         figure
%         subplot(1,3,1)
%         imshow(maskedCellFull(:,:,iiFrame),[])
%         hold on 
%         plot(posTemp(:,1),posTemp(:,2),'*r')
%         subplot(1,3,2)
%         imshow(maskedCellTop(:,:,iiFrame),[])
%         hold on 
%         plot(posTemp(:,1),posTemp(:,2),'*r')
%         subplot(1,3,3)
%         imshow(maskedCellBot(:,:,iiFrame),[])
%         hold on 
%         plot(posTemp(:,1),posTemp(:,2),'*r')
        try 
      [intTopBright,intBotBright,~,~] =  intensityBead(maskedCellTop(:,:,iiFrame),maskedCellBot(:,:,iiFrame),posTemp(1,:),9);
      [intTop,intBot,~,~] =  intensityBead(maskedCellTop(:,:,iiFrame),maskedCellBot(:,:,iiFrame),posTemp(:,:),9);
        
        intFTopB(iiFrame)=intTopBright;
        intFBotB(iiFrame)=intBotBright;
        intFTopAll=cat(2,intFTopAll,intTop);
        intFBotAll=cat(2,intFBotAll,intBot);
        ttemp = ones(length(intTop),1)*iiFrame;
        tIndAll =cat(1, tIndAll,ttemp);
        tIn(iiFrame) = iiFrame;
        catch
        intFTopB(iiFrame)=NaN;
        intFBotB(iiFrame)=NaN;  
        tIn(iiFrame) = NaN
        end
      
     
     end 

 %%

xInt = intFTopB./intFBotB
minIratio = min(xInt);
maxIratio = max(xInt);
xWeight =(xInt -min(xInt))./(max(xInt)-min(xInt))
 figure
 plot(tIn,xWeight,'*') 
sum(xWeight<0.5)
%Now try it out on everyone
xInt = intFTopAll./intFBotAll;
minIratio = min(xInt);
maxIratio = max(xInt);
xAll =(xInt -minIratio)./(maxIratio-minIratio);
 figure
 plot(tIndAll,xAll,'*') 
sum(xAll<0.5)

intVec=[intFTopAll',intFBotAll', tIndAll];
 
%%
clear fit1f
clear fit2f
clear fitCf

for iiFrame =1:5
    try
    8
    iiFrame
roughXY =posXYrot;
rndposX = round(roughXY(roughXY(:,3)==iiFrame,1));
rndposY = round(roughXY(roughXY(:,3)==iiFrame,2));

    
tempRotCbot = rotCbot(:,:,:,iiFrame);%-0.45*max(rotCbot(:))*(fittedCellRot);
tempRotCtop = rotCtop(:,:,:,iiFrame);%-0.45*max(rotCtop(:))*(fittedCellRot);
tempRotCbot = tempRotCbot/(max(rotCbot(:)));
tempRotCtop = tempRotCtop/(max(rotCtop(:)));


tempRotCbot =tempRotCbot.* double(tempRotCbot>0);
tempRotCtop =tempRotCtop.*double(tempRotCtop>0);
%tempRotCtop= tempRotCtop/max(tempRotCtop(:));
%tempRotCbot= tempRotCbot/max(tempRotCbot(:));
intToptemp =intVec(intVec(:,3) ==iiFrame,1);
intBottemp=intVec(intVec(:,3) == iiFrame,2);
xInt = intToptemp./intBottemp;
xtemp =(xInt -minIratio)./(maxIratio-minIratio);

iiElement=1;
tempRotCfull= (1-xtemp(iiElement))*tempRotCbot + xtemp(iiElement)*tempRotCtop;
  

sizeC = size(tempRotCtop);
newfull= reshape(tempRotCfull(rndposY(iiElement),:,:),[sizeC(2),sizeC(3)])';
   
fitP(1)= size(newfull,2)-12;
fitP(2)= 1;
fitP(3) = size(newfull,2)/2;
fitP(4) =23;
fitP(5) =  size(newfull,2)/2;
fitP(6) = 29;
fitP(7) =10;
fitP(8) = 0.7;
%fitP(9) = 0.7;
 fitZConvData(fitP,newfull);
% 
    try
         [fittedConv, params,fit1,fit2]=cellFitZConv (newfull,fitP);

         fit1f(:,:,iiFrame)=fit1;
         fit2f(:,:,iiFrame)=fit2;
         fitCf(:,:,iiFrame)=fittedConv;
    catch
        if iiFrame ==1
            
            
            
        else
        fit1f(:,:,iiFrame)=nan(size(fit1f(:,:,iiFrame)));
        fit2f(:,:,iiFrame)=nan(size(fit2f(:,:,iiFrame)));
        fitCf(:,:,iiFrame)=nan(size(fitCf(:,:,iiFrame)));
        end
    end 
    catch
        warning('missed Frame')
    end
end
 
fit1=mean(fit1f,3);
figure;imshow(fit1,[])
fit2=mean(fit2f,3);
figure; imshow(fit2,[])
fittedConv=mean(fitCf,3);
figure;imshow(fittedConv,[])

%%
clear fit1save
clear fit2save
clear fitConv
clear newTopFull
clear newBotFull 
clear newTopFullSub
clear newBotFullSub 
clear newAddFull 
clear newPreAddFull


fittedCellRot=rotCcell;
% clear newTopFull
% clear newBotFull
% clear newAddFull 
% clear newPreAddFull
9
for iiFrame =1:1000
    try
    
   iiFrame
tempRotCbot = rotCbot(:,:,:,iiFrame);%-0.45*max(rotCbot(:))*(fittedCellRot);
tempRotCtop = rotCtop(:,:,:,iiFrame);%-0.45*max(rotCtop(:))*(fittedCellRot);
tempRotCbot = tempRotCbot/(max(rotCbot(:)));
tempRotCtop = tempRotCtop/(max(rotCtop(:)));

tempRotCbot =tempRotCbot.* double(tempRotCbot>0);
tempRotCtop =tempRotCtop.*double(tempRotCtop>0);
%tempRotCtop= tempRotCtop/max(tempRotCtop(:));
%tempRotCbot= tempRotCbot/max(tempRotCbot(:));
intToptemp =intVec(intVec(:,3) ==iiFrame,1);
intBottemp=intVec(intVec(:,3) == iiFrame,2);
xInt = intToptemp./intBottemp;
xtemp =(xInt -minIratio)./(maxIratio-minIratio);


%tempRotCfull =tempRotCfull.*double(tempRotCfull>0);

roughXY =posXYrot;
rndposX = round(roughXY(roughXY(:,3)==iiFrame,1));
rndposY = round(roughXY(roughXY(:,3)==iiFrame,2));

for iiElement =1:size(rndposX)
    tempRotCfull= (1-xtemp(iiElement))*tempRotCbot + xtemp(iiElement)*tempRotCtop;
  
    sizeC = size(tempRotCtop);
     newtop=xtemp(iiElement)* reshape(tempRotCtop(rndposY(iiElement),:,:),[sizeC(2),sizeC(3)])';
    % figure; imshow(newtop,[])
     newbot= (1-xtemp(iiElement))*reshape(tempRotCbot(rndposY(iiElement),:,:),[sizeC(2),sizeC(3)])';
   % figure; imshow(newbot,[])
     % newfull = newtop+newbot;
     newfull= reshape(tempRotCfull(rndposY(iiElement),:,:),[sizeC(2),sizeC(3)])';
    % figure;imshow(newfull,[])
    %max(newtop(:))
     %     max(newbot(:))
     %max(newfull(:))


%figure; imshow(newtop,[])


 
%  figure; imshow(fit1,[])
 newfull2=newtop-max(newtop(:))*fit1 +newbot -max(newbot(:))*fit2;

 %figure ; imshow(newfull2<0,[])
 
fit1save(:,:,iiElement,iiFrame)=fit1;
fit2save(:,:,iiElement,iiFrame)=fit2;
fitConv(:,:,iiElement,iiFrame)=fittedConv;
newTopFull (:,:,iiElement,iiFrame)= newtop;
newBotFull (:,:,iiElement,iiFrame)= newbot;
newTopFullSub (:,:,iiElement,iiFrame)= newtop-max(newtop(:))*fit1;
newBotFullSub (:,:,iiElement,iiFrame)= newbot-max(newbot(:))*fit2;
newAddFull (:,:,iiElement,iiFrame)=newfull2;
newPreAddFull (:,:,iiElement,iiFrame)= newfull-max(newfull(:))*fittedConv;

% 
% 
%      figure
%     subplot (1,4,1)
%      imshow(newtop-fit1*max(newtop(:)),[0,max(newfull2(:))])
%      title('Top')
%       
%  subplot (1,4,2)
%      imshow(newbot-fit2*max(newbot(:)),[0,max(newfull2(:))])
%      title('Bot')
%      
%  subplot (1,4,3)
%      imshow(newfull2,[0,max(newfull2(:))])
%      title('Sum')
%      hold on
% %      plot(tempPeaksFull(:,1),tempPeaksFull(:,2),'r*')
%  subplot (1,4,4)
%      imshow(newfull-max(newfull(:))*fittedConv,[0,max(newfull2(:))])
%      title('Sum')
%      hold on
%     plot(tempPeaksFull(:,1),tempPeaksFull(:,2),'r*')
%     plot(tempPeaksFull2(:,1),tempPeaksFull2(:,2),'b*')
%figure; imshow(rotCtopFilt(:,:,23,1),[])


% a =newtop%-0.7*fit1
% b =newbot%-0.7*fit1
% c= a+b;
% d=newtop+newbot-fittedConv
% 
% scale = max(max(a(:)),max(b(:)))
%      figure
%     subplot (1,4,1)
%      imshow(a,[0, scale])
%      title('Top')
%       
%  subplot (1,4,2)
%      imshow(b,[0,scale])
%      title('Bot')
%      
%  subplot (1,4,3)
%      imshow(c/max(c(:)),[0,1])
%      title('Sum')
%      hold on
%     
%  subplot (1,4,4)
%      imshow(d/max(d(:)),[0,1])
%      title('Sum')
%      hold on
%   







end

    catch
    end
end

%%


index=1
clear posXYZ
10
for iiFrame =1:1000
    iiFrame
try
    tempTop =newTopFullSub(:,:,:,iiFrame);%/max(newTopFullSub(:));
    tempBot =newBotFullSub(:,:,:,iiFrame);%/max(newBotFullSub(:));
    tempAdd =tempTop+tempBot;
    tempPreAdd =newPreAddFull(:,:,:,iiFrame);
    roughXY =posXYrot;
    posX = (roughXY(roughXY(:,3)==iiFrame,1));
    posY = (roughXY(roughXY(:,3)==iiFrame,2));
    rndposX = round(roughXY(roughXY(:,3)==iiFrame,1));
    rndposY = round(roughXY(roughXY(:,3)==iiFrame,2));
    time= round(roughXY(roughXY(:,3)==iiFrame,3));

    %figure; imshow (Irot (:,:,iiFrame),[])
    %hold on 
    %plot(rndposX,rndposY, 'b.')
    clear zCom
    clear varZcom
    clear posZ
    
    
    for iiElement=1:size(rndposX)
        try
        iiElement
        %figure; imshow(tempTop(:,:,iiElement),[])
        tempC = tempPreAdd(:,:,iiElement);    
        tempC= double(tempC>0).*tempC;
        
%tempC = newPreAddFull(:,:,:,iiFrame);

        sliverX = tempC(:,rndposX(iiElement));
        
        
        %figure; imshow(sliverX,[])
        %roughSliverY = reshape(tempC(ind2(iiElement),ind1(iiElement)-3:ind1(iiElement)+3,:),[7,size(tempC,3)])';

        [~,midx]=max(sliverX)
        cropSliver = sliverX(midx-5:midx+5);
        zIdxPos = 1:size(tempC,1);
        cropZidx = zIdxPos(midx-5:midx+5);
        
        [fitG,gofG]=createGaussianFit1D(cropZidx, cropSliver);
        y = normpdf(cropZidx,fitG.b1,fitG.c1);
        
         if fitG.b1 <20 | fitG.b1>30
            posZ(iiElement) =NaN;
         else
             posZ(iiElement)= fitG.b1;
         end
        
        
%         figure; plot(cropZidx,cropSliver)
% 
%         figure; 
%         subplot (1,5,1); imshow(tempTop(:,:,iiElement),[])
%         hold on
%         plot(posX(iiElement),fitG.b1,'r*')
%         plot(posX(iiElement),26,'g*')
%         subplot (1,5,2); imshow(tempBot(:,:,iiElement),[])
%         hold on
%         plot(posX(iiElement),fitG.b1,'r*')
%         plot(posX(iiElement),26,'g*')
%         subplot (1,5,3); imshow(tempPreAdd(:,:,iiElement),[])
%         hold on
%         plot(posX(iiElement),fitG.b1,'r*')
%         plot(posX(iiElement),26,'g*')
%         subplot (1,5,4); imshow(tempC,[])
%         hold on
%         plot(posX(iiElement),fitG.b1,'r*')
%         plot(posX(iiElement),26,'g*')
%         subplot (1,5,5); imshow(Irot(:,:,iiFrame),[])
%         hold on
%         plot(posX(1),posY(1),'b*')
        
        
        indexes= flip(cropZidx);
        %zCom(iiElement) = sum(indexes'.*cropSliver)/sum(cropSliver)
        %std(indexes, cropSliver)
        %varZcom (iiElement) =sqrt(sum( (indexes.^2)'.*cropSliver )/sum(cropSliver)- zCom(iiElement).^2);

       % hold on
        %plot (posX(1), 26+zeros(size(posX(1))),'b.')
        
       
            
        
        catch
            warning(strcat('iiElement',num2str(iiElement),'in frame',num2str(iiFrame) ,' did not work'))
            posZ(iiElement)=NaN;
        end
    end
    try
        intToptemp =intVec(intVec(:,3) ==iiFrame,1);
        intBottemp=intVec(intVec(:,3) == iiFrame,2);
        xOg=resXYtrans(intVec(:,3) ==iiFrame,1);
        yOg=resXYtrans(intVec(:,3) ==iiFrame,2);
        posXYZrot(index:index+size(rndposX,1)-1,:)=[posX,posY,posZ',intToptemp,intBottemp,time];
        posXYZ(index:index+size(rndposX,1)-1,:)=[xOg,yOg,posZ',intToptemp,intBottemp,time];
        
        
        index = index + size(rndposX,1);%
        %posXYZ(index:index+size(rndposX(1),1)-1,:)=[posX(1),posY(1),posZ',iiFrame];
        %index = index + 1;%
    catch
    posXYZrot(index,:)=[NaN,NaN,NaN,NaN,NaN,iiFrame] ;   
    posXYZ(index,:)=[NaN,NaN,NaN,NaN,NaN,iiFrame] ;   
    index=index+1;
    end
    
    
catch
    warning('missing frame')
end 
end
%%
posXYZrot(posXYZrot(:,4)==0,:)=[];
posXYZrot(isnan(posXYZrot(:,3)),:)=[];
posXYZ(posXYZ(:,4)==0,:)=[];
posXYZ(isnan(posXYZ(:,3)),:)=[];

% figure;plot((posXYZ(:,2)-min(posXYZ(:,2)))*80,((posXYZ(:,3)-min(posXYZ(:,3)))*100),'.r')
% figure;plot((posXYZ(:,3)*100),(posXYZ(:,2))*80,'.')

dlmwrite(strcat(saveFileGrainRot,'.csv'),posXYZrot,'delimiter',',','precision',32);
dlmwrite(strcat(saveFileNameGrainCropped,'.csv'),posXYZ,'delimiter',',','precision',32);
  
   
%%
save(strcat(saveFileNameGrainCropped,'_vars.mat'))



%%

% xPos=(posXYZ(:,1));%-min(posXYZ(:,1)));
% yPos=(posXYZ(:,2));%-min(posXYZ(:,2)));
% zPos=(posXYZ(:,3));%-min(posXYZ(:,3)));
% 
% 
% xyzPos=[xPos,yPos,zPos,posXYZ(:,4)];
%  stitched=stitchingTrajConv(xyzPos);
% figure;plot(stitched(1:100,1)*80,stitched(1:100,2)*80,'.')
% figure;plot(stitched(1:100,2)*80,stitched(1:100,3)*100,'.')

% for iiFrame=1:30%length(posXYZ)
% figure; imshow(Irot(:,:,iiFrame),[])
%         hold on
%         
%         posTemp = stitched(stitched(:,4)==iiFrame,:);
%          plot(posTemp(:,1),posTemp(:,2),'b.')
% end







%%

% 
%   
% for iiElement = 1:length(rndposX)
% 
% %might be worth to do this for the intensity calcs.   
% %tempBead = temp (roughPeaks1(iiElement,2)-BOX_SIZE:roughPeaks1(iiElement,2) +BOX_SIZE,roughPeaks1(iiElement,1)-BOX_SIZE:roughPeaks1(iiElement,1) +BOX_SIZE);
% 
% 
% sliverX = tempC(:,rndposX(iiElement));
% figure; imshow(sliverX,[])
% %roughSliverY = reshape(tempC(ind2(iiElement),ind1(iiElement)-3:ind1(iiElement)+3,:),[7,size(tempC,3)])';
% 
% indexes= flip(1:size(tempC,1));
% zCom(iiElement) = sum(indexes'.*sliverX)/sum(sliverX)
% std(indexes, sliverX)
% varZcom (iiElement) =sqrt(sum( (indexes.^2)'.*sliverX )/sum(sliverX)- zCom(iiElement).^2);
% end
% 
% % if (length(com))>1
% %     plot(((com)*100/0.59),(subPeaks1(:,2))*80,'.r')
% % else
% %     plot(((com-min(com))*100/0.59),(subPeaks1(:,2)-min(subPeaks1(:,2)))*80,'.b')
% % end
% %figure;imshow(roughSliverY,[]);hold on;plot(4,com,'r*');
% posXYZ(index:index+size(rndposX,1)-1,:)=[rndposX,rndposY,zCom',varZcom',time];
% index = index + size(rndposX,1);%
% end

%%

