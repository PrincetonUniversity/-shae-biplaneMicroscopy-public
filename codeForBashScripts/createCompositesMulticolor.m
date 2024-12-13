function [kk,positionList] = createCompositesMulticolor(f, st, stprime, tiffList, jjImage)
%%
clc
%clear all
% close all

%% constants
BOX_SIZE_CUTOUT = 3;
BOX_SIZE_PEAK = 7;
CAMERA_DIVIDER_PTS = [[1,1,512,512,1]',[1,270,270,1,1]'];
%CAMERA_DIVIDER_PTS = [[1,1,414,414,1]',[1,220,220,1,1]'];
%MAPPING_ERROR_MAX = 2;
IMAGE_PIXEL_SIZE = 512;

jj=1;

%%
    fileNameDisp = eraseBetween(char(tiffList{jjImage}),1,'2');
    [filepath,name,~]=fileparts(fileNameDisp);
    mkdir(filepath);
    filename=strcat(fullfile(filepath,name),'_Positions_test.csv');
    warning('No fileName input, pairs will be saved to %s', filename);
    tiffName = tiffList{jjImage};
    info = imfinfo(tiffName);
    numImages =  numel(info);

    %%
    savedComp=zeros(IMAGE_PIXEL_SIZE,IMAGE_PIXEL_SIZE,numImages);

    %%
    limits = CAMERA_DIVIDER_PTS;
    clear peaksRedNew2
    PEDESTAL =99;

    index =1;

    
  
    for kk = 1:numImages
          kk
     %   try
        tempImage = imread(tiffName, kk, 'Info', info);
        tempImage= tempImage-PEDESTAL; %this subtracts the pedestal value
       % imshow(tempImage,[])
       % hold on
       [pathMD,~,~] = fileparts(tiffName);
       pathMD;
        %tempImage = resizeImageMulticolor(tiffName,pathMD, tempImage,IMAGE_PIXEL_SIZE);
        filteredImage = bpass(double(tempImage),0.1,15);
%         if (kk==1)
%             figure
%             imshow(filteredImage(:,:),[],'initialMagnification','fit');
%             hold on
%         end
        if sum(savedComp(:))~=0
            compositeImage =savedComp(:,:,index); 
        end
        
       
        [compositeImage,newImage] = biplaneCompImage(tempImage(:,:), st, limits, IMAGE_PIXEL_SIZE);
        filteredComp = bpass(double(compositeImage),0.30,10);
        savedComp(:,:,index)=compositeImage;
        savedNew(:,:,:,index)=newImage;
        index = index +1;
      

    end
  
varFileName=strcat(filepath,'/',name,'_variables_Voll5_size7.mat');
save(varFileName);


end