function mergeBackImages(parentDir,imgSubDirs,indexDir)
    %parentDir = uigetdir();
    %%
    %imgSubDirs=getDirectoryNames(parentDir)

    imgSubDirs2= imgSubDirs
    
    %%
    %%for indexDir = 1:length(imgSubDirs)
        imageDirPath = fullfile(parentDir, char(imgSubDirs(indexDir)));
        imageDirPath= fullfile(imageDirPath ,char(getDirectoryNames(imageDirPath)));
        matFiles2 =dir (fullfile(imageDirPath,'*.mat'));
        savedCompFull= [];
        savedNewFull = [];
        for iiMatFiles=1:size(matFiles2)
            matNames(iiMatFiles) ={matFiles2(iiMatFiles).name}
        end
        sortIdx = numericalIndexSorting(matNames,'\d*_variables')
        matFiles2=matFiles2(sortIdx)

        for iiMatFiles=1:size(matFiles2)

            if iiMatFiles ==1
                load(fullfile(matFiles2(iiMatFiles).folder,matFiles2(iiMatFiles).name))
                2
            else 
                load (fullfile(matFiles2(iiMatFiles).folder,matFiles2(iiMatFiles).name),'savedComp')
                load (fullfile(matFiles2(iiMatFiles).folder,matFiles2(iiMatFiles).name),'savedNew')
                56
            end
            savedCompFull = cat(3,savedCompFull,savedComp);
            savedNewFull = cat(4,savedNewFull,savedNew);
        end
        savedComp = savedCompFull;
        savedTop = savedNewFull(:,:,2,:);
        savedBot = savedNewFull(:,:,1,:);
        %clear savedCompFull

        imgSubDirs2
        saveFolder = parentDir;
        cropInd = strfind(saveFolder,fullfile('code','Diana'));
        saveFolder(cropInd-1:cropInd+9)=[];
        imgSubDirs2
        fullfile(saveFolder,char(imgSubDirs2(indexDir)))
        save(fullfile(saveFolder,char(imgSubDirs2(indexDir)),'variables.mat'))

   
end
%%