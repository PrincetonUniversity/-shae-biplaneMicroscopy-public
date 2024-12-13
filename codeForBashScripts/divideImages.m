%expDir = uigetdir();

function divideImages(expDir, iiFile)
%%
subExp = getDirectoryNames(expDir);

for iiExp= 1:numel(subExp)
    iiExp
    zStackPath = fullfile (expDir,char(subExp(iiExp)))

    imagePathTxt =fullfile(zStackPath,'*.tif');
    imageName = dir(imagePathTxt);

    for iiStep =1:length(imageName)
        if strfind(imageName(iiStep).name,'GEM')
             tifList (iiExp) = {fullfile(imageName(iiStep).folder,imageName(iiStep).name)};
        end
    
        if strfind(imageName(iiStep).name,'trans')
             transList(iiExp) = {fullfile(imageName(iiStep).folder,imageName(iiStep).name)};
        end

    end
    
    imagePathMD =fullfile(zStackPath,'*.txt');
    imageMD = dir(imagePathMD);
    
    for iiStep = 1:length(imageMD)
        if strfind(imageMD(iiStep).name,'GEM')
            MDList (iiExp) = {fullfile(imageMD(iiStep).folder,imageMD(iiStep).name)};
        end
    end
        
end    

%%
%for iiFile =1: size(tifList,2) %1: size(tifList,2)
    iiFile
    fullTiffName = char(tifList(iiFile));
    mdName = char(MDList(iiFile));
    filepath = eraseBetween(fullTiffName,1,'2');
    
   
    info = imfinfo(fullTiffName);
    numImages =  numel(info);
    [pathI,nameI,extI]=fileparts(fullTiffName);
    imagePath = fullfile(pathI,nameI);
    mkdir(imagePath)
  
    for  iisubFiles= 1:ceil(numImages/100)
        iisubFiles
       
        
        for iiImages =(iisubFiles-1)*100+1:iisubFiles*100
            tempI = imread(fullTiffName, iiImages);
            
            if (iiImages == (iisubFiles-1)*100+1) 
                copyfile (mdName ,fullfile(imagePath, strcat( num2str(iisubFiles),'.txt')))
                imwrite(tempI, fullfile(imagePath, strcat(num2str(iisubFiles),'.tif')))
            else
                copyfile (mdName,fullfile( imagePath,strcat( num2str(iisubFiles),'.txt')))
                imwrite(tempI,fullfile(imagePath, strcat(num2str(iisubFiles),'.tif')),'WriteMode','append')
            end
        end 
      
        
    end  
%end 
end
%%