function [psfTop,psfBot,psfGTop,psfGBot] = getAggregatedPSF (parentDir)

BOX_SIZE_PEAK = 9;
CAMERA_DIVIDER_PTS = [[1,1,512,512,1]',[1,270,270,1,1]'];
IMAGE_PIXEL_SIZE = 512;
PEDESTAL = 99;

%%
%First get file names
subDirs = getDirectoryNames(parentDir);
exts = '.tif';
getNewTransform = false;
displayImage =false;

filepath = eraseBetween(char(parentDir),1,'2');
mkdir(filepath);

for iiDirs= 1:numel(subDirs)
    fullPathTif= cell2mat(fullfile(parentDir, subDirs(iiDirs),'*.tif'));
    fullPath= cell2mat(fullfile(parentDir, subDirs(iiDirs)));
    dd = dir(fullPathTif);
    tifListTemp = {dd.name};
    tifList(:,iiDirs) = fullfile(fullPath,tifListTemp);
end

%%

for jjDirs = 1:numel(subDirs)
    jjDirs

    clear aggregatedImage
    clear aggregatedFiltImage
    for kkFile = 1:numel(tifList(:,jjDirs))
       tempImage = imread(tifList{kkFile,jjDirs},1); 
       tempImage= tempImage-PEDESTAL;
       %figure; imshow(tempImage,[])
       aggregatedFiltImage(:,:,kkFile) = bpass(double(tempImage),0.5,15); % spatial bandpass filter 
       
    end
    maxProj = max(aggregatedFiltImage,[],3); % z project
     if(jjDirs==1)
          figure
             imshow(maxProj,[],'initialMagnification','fit');
             hold on
    end

%% find peaks


    if jjDirs ==1
        if exist('tformOut2','var')
            [cpPairs,tformOut2]= extractPairs(tifList(:,jjDirs),tformOut2, getNewTransform ,displayImage);
        else 
            [cpPairs,tformOut2]= extractPairs(tifList(:,jjDirs),'tformOut2', true ,displayImage);
        end
    else
        [cpPairs,tformOut2]= extractPairs(tifList(:,jjDirs),tformOut2, false ,displayImage);
    end

    peaksBot = cpPairs(:,3:4);
    clear peaksTop;
    [peaksTop(:,1),peaksTop(:,2)]= transformPointsInverse(tformOut2,...
        peaksBot(:,1),peaksBot(:,2));

    if (jjDirs==1)
     plot(peaksTop(:,1),peaksTop(:,2),'.r')
     hold on
     plot(peaksBot(:,1),peaksBot(:,2),'b.')
    end

%%
    clear tempBeadRed
    clear tempBeadBottom
    for iiBead = 1:size(peaksBot,1)
         
        tempBeadTop(:,:,:,iiBead) =  aggregatedFiltImage(...
            round(peaksTop(iiBead,2))-BOX_SIZE_PEAK:round(peaksTop(iiBead,2))+BOX_SIZE_PEAK,...
            round(peaksTop(iiBead,1))-BOX_SIZE_PEAK:round(peaksTop(iiBead,1))+BOX_SIZE_PEAK,:);
        tempBeadBot(:,:,:,iiBead) =  aggregatedFiltImage(...
            round(peaksBot(iiBead,2))-BOX_SIZE_PEAK:round(peaksBot(iiBead,2))+BOX_SIZE_PEAK,...
            round(peaksBot(iiBead,1))-BOX_SIZE_PEAK:round(peaksBot(iiBead,1))+BOX_SIZE_PEAK,:);
    end
%%
%      clear avgBeadTop
%      clear avgBeadBottom
    avgBeadTop(:,:,:,jjDirs) = sum(tempBeadTop,4);
    avgBeadBottom(:,:,:,jjDirs) = sum(tempBeadBot,4);
    THRESHOLD = 0;

    sizeA = size(avgBeadTop)
    A= reshape(avgBeadTop(floor(sizeA(1)/2),:,:,jjDirs),[sizeA(1),sizeA(3)]);
%     figure
%     subplot (1,2,1)
%      imshow(A,[])
%      title('Top')
%       A= reshape(avgBeadBottom(floor(sizeA(1)/2),:,:,jjDirs),[sizeA(1),sizeA(3)])
%  subplot (1,2,2)
%      imshow(A,[])
%      title('Bot')
     
     
 end
 
 
 
 %%
avgBeadTopComp=sum(avgBeadTop,4);
avgBeadBotComp=sum(avgBeadBottom,4);
sizeA = size(avgBeadTopComp);
Atop= reshape(avgBeadTopComp(floor(sizeA(1)/2),:,:),[sizeA(1),sizeA(3)]);
figure
subplot (1,2,1)
imshow(Atop,[])
title('Top')
Abot= reshape(avgBeadBotComp(floor(sizeA(1)/2),:,:),[sizeA(1),sizeA(3)]);
subplot (1,2,2)
imshow(Abot,[])
title('Bot')
     

%%
figure
C = imfuse(Atop,Abot,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
imshow(C)


%%

%We want our kernel for convolution to be 21 by 21 by 51 so that it's more
%or less padded by 1 psf. Based on looking at the images found that the
%center should be at pixel 26. See page 3 of lab book. 
psfTop = zeros(21,21,51);
psfBot = zeros(21,21,51);


psfTop (2:20,2:20,6:45)=avgBeadTopComp;
psfBot (2:20,2:20,6:45)=avgBeadBotComp;

sigma = 1;
for iiPlane = 1:size(psfTop,3)
psfGxyTop(:,:,iiPlane) = imgaussfilt(psfTop(:,:,iiPlane),1);
psfGxyBot(:,:,iiPlane) = imgaussfilt(psfBot(:,:,iiPlane),1);
end

for iiRow = 1:size(psfTop,1)
    for iiCol =1:size (psfTop,2)
        sliceTop =psfGxyTop(iiRow,iiCol,:);
        sliceTop = imgaussfilt (sliceTop(:),2);
        sliceBot =psfGxyBot(iiRow,iiCol,:);
        sliceBot = imgaussfilt (sliceBot(:),2);
        psfGTop(iiRow,iiCol,:) = reshape(sliceTop, [1,1,length(sliceTop)]);
        psfGBot(iiRow,iiCol,:) = reshape(sliceBot, [1,1,length(sliceBot)]);
    end
end
% figure
% subplot(1,2,1)
% imshow(psfGTop(:,:,23),[]);
% subplot(1,2,2)
% imshow(psfGBot(:,:,29),[]);


sizeA = size(psfGTop);
Atop= reshape(psfGTop(floor(sizeA(1)/2),:,:),[sizeA(1),sizeA(3)])
figure
subplot (1,2,1)
imshow(Atop,[])
title('Top')
Abot= reshape(psfGBot(floor(sizeA(1)/2),:,:),[sizeA(1),sizeA(3)])
subplot (1,2,2)
imshow(Abot,[])
title('Bot')