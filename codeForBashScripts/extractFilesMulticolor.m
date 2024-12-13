% This function finds all the files required for tracking.
% in: expDir 
% out : fileStruct
%       fileStruct.imagePath
%       fileStruct.nucleoidParent
%       fileStruct.membraneParent
%       fileStruct.variableList
%       fileStruct.MDList 
%       fileStruct.transList 

function fileStruct = extractFilesMulticolor (expDir)
expDir
    subExp = getDirectoryNames(expDir)

    for iiExp= 1:numel(subExp)
        iiExp
        zStackPath = fullfile (expDir, char(subExp(iiExp)))

        zStacks = getDirectoryNames(zStackPath)
        try
        a=strfind(zStacks,'membrane');
        for iiStep =1:length(a)
             if ~isempty(a{iiStep})
                  membraneParent(iiExp) = fullfile(zStackPath, zStacks(iiStep)) ;
             end
        end
        catch
            warning("No membrane image")
        end

        a=strfind(zStacks,'nucleoid');
        for iiStep =1:length(a)
             if ~isempty(a{iiStep})
                  nucleoidParent(iiExp) = fullfile(zStackPath, zStacks(iiStep)) ;
             end
        end

        imagePathTxt =fullfile(zStackPath,'*.tif');
        imageName = dir(imagePathTxt);

        for iiStep =1:length(imageName)
            if strfind(imageName(iiStep).name,'GEM')
                 tifList (iiExp) = {fullfile(imageName(iiStep).folder, imageName(iiStep).name)};
            end

            if strfind(imageName(iiStep).name,'trans')
                 transList(iiExp) = {fullfile(imageName(iiStep).folder, imageName(iiStep).name)};
            end

        end

        imagePathMD =fullfile(zStackPath,'*.txt');
        imageMD = dir(imagePathMD);

        for iiStep = 1:length(imageMD)
            if strfind(imageMD(iiStep).name,'GEM')
                MDList (iiExp) = {fullfile(imageMD(iiStep).folder, imageMD(iiStep).name)};
            end
        end


        imagePathVars =fullfile(zStackPath,'*.mat');
        imageVar = dir(imagePathVars);

        for iiStep = 1:length(imageVar)
            if strfind(imageVar(iiStep).name,'variables')
                variableList (iiExp) = {fullfile(imageVar(iiStep).folder, imageVar(iiStep).name)};
            end
        end

    end    

    %%
    fileStruct.imagePath = tifList;
    fileStruct.nucleoidParent = nucleoidParent;
    try
    fileStruct.membraneParent = membraneParent;
    catch
        warning('No membrane files available')
    end
    
    fileStruct.MDList = MDList;
    try
        fileStruct.variableList = variableList;
    catch
        warning('No variable files available')
    end
    
    try
        fileStruct.transList = transList;
    catch
        warning('No transmitted files available')
    end


    %%
    fullfile(expDir,'filePaths.mat')
    save(fullfile(expDir,'filePaths.mat'),'fileStruct')

end

%expDir = uigetdir();