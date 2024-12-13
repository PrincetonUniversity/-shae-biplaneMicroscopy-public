function [varsList,trajFilesListProb,trajFilesList] = trackingParticlesprobConvFunction(saveDirParent)
    %saveDirParent = uigetdir(); %code directory
    %saveDirParent= 'Z:\dsmendez\diana\code\Diana\20211013_newFilterOvnights\AqLS\100msCopy'
    %%
    saveDirParent
    saveDirs = getDirectoryNames(saveDirParent)
    counter =1;
    clear varsList
    clear saveDirList

    for iiExp= 1:numel(saveDirs) 
        saveDirTemp = fullfile(saveDirParent,char(saveDirs(iiExp)))
        fileNameStr = dir(saveDirTemp)

        for iiFile =1:length(fileNameStr)
            if strfind(fileNameStr(iiFile).name,'vars.mat')
                iiFile
                saveDirList (counter) = {fileNameStr(iiFile).folder};
                varsList (counter) = {fullfile(fileNameStr(iiFile).folder,fileNameStr(iiFile).name)};
                cropInd = strfind(saveDirTemp,'code');
                dataDir = saveDirTemp;
                dataDir(cropInd:cropInd+10)=[];
                counter = counter +1;
            end

        end

    end    
    %%
    clear order1
    %pattern='cell_int_\d*'
    pattern='cell_\d*';

    order=regexp(varsList, pattern, 'match');
    order=string(order);
    order=regexp(order, '\d*', 'match');
    for i=1:length(varsList);
        order1(i)=str2num(char(order{i}));
    end
    [~,indexx]=sort(order1);
    varsList=varsList(indexx);
    saveDirList = saveDirList(indexx);
    %the code above teaches matlab how to count like this 1,2,3...
    %:)
    %clear trajFilesListProb


    %%
    counterBad =zeros(length(varsList),1);
    length(varsList)
    for iiCell=  1:length(varsList)

        %try
        iiCell;
        [~,tempStr,~]=fileparts(char (varsList(iiCell)))
        indx = regexp (tempStr,'\d');
        cellNum= tempStr(indx);

        filepath = char(saveDirList(iiCell));
        mkdir(filepath);
        rotTrajFolder = strcat (filepath);

        saveFileNameTrajProb = fullfile(rotTrajFolder,'traj_cell_prob');
        saveFileNameTraj = fullfile(rotTrajFolder,'traj_cell_');
        saveFileNameRot = fullfile(rotTrajFolder,'rotated_cell_');
        saveFileGrainRot = strcat(saveFileNameRot,num2str(cellNum));
        posRotated = csvread(strcat(saveFileGrainRot,'.csv'));
        saveFileGrainTrajProb = strcat(saveFileNameTrajProb,num2str(cellNum));
        saveFileGrainTraj = strcat(saveFileNameTraj,num2str(cellNum));
        
        trajFilesListProb (iiCell) = {saveFileGrainTrajProb};
        trajFilesList (iiCell) = {saveFileGrainTraj};
     
        stitchedProb=stitchingTrajProbConvFit(posRotated);
        resProb= stitchedProb (stitchedProb(:,4)~=0,:);
        dlmwrite(strcat(saveFileGrainTrajProb,'.csv'),resProb,'delimiter',',','precision',32);

        stitched=stitchingTrajConv(posRotated);
        res= stitched (stitched(:,4)~=0,:);
        dlmwrite(strcat(saveFileGrainTraj,'.csv'),res,'delimiter',',','precision',32);
    end
    
    
fullfile(saveDirParent,'files2Fix.mat')
save(fullfile(saveDirParent,'files2Fix.mat'),'varsList','trajFilesListProb','trajFilesList')
 
    
end


