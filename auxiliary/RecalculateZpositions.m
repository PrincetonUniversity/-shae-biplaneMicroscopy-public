%   positionFiles = uipickfiles();
%   varList =uipickfiles();

%positionFiles = trajFilesListProb;
iiFile=3
%traj = csvread(strcat(char(positionFiles(iiFile)),'.csv')); 
traj = csvread(strcat(char(positionFiles(iiFile))));%,'.csv')); 


%for iiFile =1:length(positionFiles)
%traj = csvread(char(positionFiles(iiFile))); 
    


Ivars =load(char(varList(iiFile)),'Irot','cropImBot','cropImTop','th',...
    'saveFileNameTraj','iiCell','newTopFullSub','newBotFullSub','newPreAddFull',...
    'cellMask','cropImTop','cropImBot','fittedCellNoConv','rotCMask');

%Ivars =load(char(varsList(iiFile)),'Irot','cropImBot','cropImTop','th',...
%    'saveFileNameTraj','iiCell','newTopFullSub','newBotFullSub','newPreAddFull',...
%    'cellMask','cropImTop','cropImBot','fittedCellNoConv','rotCMask');
%%
figure;imshow(Ivars.rotCMask,[])
se1 = strel('line',1,0)
se2 = strel('line',1,90)
SE=strel('disk',1);
cellMask = imdilate(Ivars.rotCMask,[se1 se2],'full')
cellMask = imdilate(cellMask,SE)%[se1 se2],'full')
cellMask = imdilate(cellMask,[se1 se2],'full')
cellMask = imdilate(cellMask,SE)%[se1 se2],'full')
figure;imshow(cellMask,[])


%%
clear maskedCellFull

for iiFrame = 1:1000
    iiFrame
    
    tempTop = Ivars.cropImTop(:,:,iiFrame);
    tempBot = Ivars.cropImBot(:,:,iiFrame);
    imageInNoConv= tempTop+tempBot;

    resMatNoConv = imageInNoConv -0.6*max(imageInNoConv(:))*Ivars.fittedCellNoConv; 
    resMatNoConvRot = imrotate(resMatNoConv,-Ivars.th,'bicubic');
    maskedCellNoConv =resMatNoConvRot.*cellMask;
    maskedCellFull(:,:,iiFrame)=maskedCellNoConv;
    %imshow(maskedCellFull(:,:,iiFrame),[])
    %drawnow
    %maskedCellTop(:,:,iiFrame)=tempTop.*cellMask;
    %maskedCellBot(:,:,iiFrame)=tempBot.*cellMask;

%Irot = imrotate(maskedCellFull,-Ivars.th);
end
%%

Irot = maskedCellFull;
 for iiFrame =1:1000
     iiFrame
     
    temp = traj (traj(:,6)==iiFrame,:);
    if isempty(temp)
        temp = newTraj (newTraj(:,6)==iiFrame-1,:);
    end
        
    filteredImage = Irot(:,:,iiFrame);
         
    [tempPeaks,~]= centfinderDVM (filteredImage, temp(1:2),9);
    newTraj(iiFrame,:)=[tempPeaks,temp(3:end-1),iiFrame];



 end
 %%
 saveNewTraj =newTraj;
 trajOg = newTraj(1:end-1,1:2);
 trajShift = circshift(newTraj(:,1:2),-1);
 trajShift = trajShift(1:end-1,:);
 stepSz = trajShift-trajOg
 
 dist = sqrt(stepSz(:,1).^2+stepSz(:,2).^2)
 
 frames =1:1000;
 thresh= quantile(dist,.95);
 missFrame = frames(dist>thresh)+1
 
 
 while (~isempty(missFrame))
    iiFrame = missFrame(1)
    temp = newTraj (newTraj(:,6)==iiFrame-1,:);
    
    filteredImage = Irot(:,:,iiFrame);
         
    [tempPeaks,~]= centfinderDVM (filteredImage, temp(1:2),9);
   
    if iiFrame <999
        test = sqrt((trajOg(iiFrame+1,1)-tempPeaks(1))^2+ (trajOg(iiFrame+1,2)-tempPeaks(2))^2)
    else
        test = 0
    end
        
    if test>=thresh
        1
        
        if ~sum(missFrame==missFrame(1)+1)
            missFrame(1) = missFrame(1)+1;
        else
            missFrame(1)=[];
        end
    else 
        2
        missFrame(1)=[];
    end
    
    newTraj(iiFrame,:)=[tempPeaks,newTraj(iiFrame,3:end)];
    
 end
    

 

  
 
%%

%traj = csvread(strcat(char(positionFiles(iiFile)),'.csv')); 

varList = uipickfiles();
positionFiles = uipickfiles();
%%
for iiFile =9%1:length(positionFiles)
traj = csvread(char(positionFiles(iiFile))); 
    
%traj = csvread(strcat(char(positionFiles(iiFile))));%,'.csv')); 
%varList = varsList;

Ivars =load(char(varList(iiFile)),'Irot','cropImBot','cropImTop','th',...
    'saveFileNameTraj','iiCell','newTopFullSub','newBotFullSub','newPreAddFull',...
    'cellMask','cropImTop','cropImBot','fittedCellNoConv','rotCMask','minBot','maxBot',...
    'maskedCellFull');

minBot = Ivars.minBot;
maxBot = Ivars.maxBot;
nucDirPath=fullfile (strcat('Z:\dsmendez\diana\', fileparts(Ivars.saveFileNameTraj)),'nucleoid0') 
memDirPath=fullfile (strcat('Z:\dsmendez\diana\', fileparts(Ivars.saveFileNameTraj)),'membrane0') 

 matFilesNuc = dir(fullfile(nucDirPath,'*.tif'));
matFilesMem = dir(fullfile(memDirPath,'*.tif'));
    for iiSlice =1:length(matFilesNuc)
        iiSlice
        nucFile=fullfile(matFilesNuc(iiSlice).folder,matFilesNuc(iiSlice).name);
        memFile=fullfile(matFilesMem(iiSlice).folder,matFilesMem(iiSlice).name);
        
        Inuc = imread(char(nucFile))-PEDESTAL;
        Imem = imread(char(memFile))-PEDESTAL;
        
        for iiRow=1:7
            Imem(iiRow,:) = Imem(8,:);
        end
        for iiCol =1:3
            Imem(:,iiCol) = Imem(:,4);
        end
                
        Imem = imtranslate(Imem,[3, 7]); %empirically seen translation to match displacement of blue channel. 
        intMem(iiSlice) = sum(sum(Imem));
        intNuc(iiSlice) = sum(sum(Inuc));
        
%         
%         cropInuc= Inuc(minBot(2):maxBot(2),minBot(1):maxBot(1));
%         cropInuc = imrotate(cropInuc,-Ivars.th,'bicubic');
%         %figure;imshow(cropInuc,[])
%         cropInuc = double(cropInuc) .* Ivars.rotCMask;
%         %figure;imshow(cropInuc,[])
%         cropImem= Imem(minBot(2):maxBot(2),minBot(1):maxBot(1));
%         cropImem = imrotate(cropImem,-Ivars.th,'bicubic');
%         %figure;imshow(cropImem,[])
%         cropImem = double(cropImem) .* Ivars.rotCMask;
%         %figure;imshow(cropImem,[])
% 
% 
%         intMem(iiSlice) = sum(sum(cropImem));
%         intNuc(iiSlice) = sum(sum(cropInuc));
    end
 

intMem = intMem./max(intMem);
intNuc = intNuc./max(intNuc) ;

[intCorrMem,fitInt]= bleachCorrection(intMem');
[intCorrNuc,fitInt]= bleachCorrection(intNuc(1:end)');
%figure; plot(intCorrMem)
%figure; plot(intCorrNuc)

intCorrNucMin = intCorrNuc - min(intCorrNuc)
[fitresult, gof,xData,yData]=createGaussianFit1D(1:size(intCorrNucMin),intCorrNucMin)
fitresult.b1


% 
% % %Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'cropSliver vs. cropZidx', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel cropZidx
% ylabel cropSliver
% grid on
% drawnow

intCorrNucMin = intCorrMem - min(intCorrMem)
[fitresult, gof,xData,yData]=createGaussianFit1D(1:size(intCorrNucMin),intCorrNucMin)
fitresult.b1

% % %Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'cropSliver vs. cropZidx', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel cropZidx
% ylabel cropSliver
% grid on
% drawnow

end
%%
%posA = positionFiles
%varA = varList 
%%
figure
varList = dataVars;
positionFiles =goodTracks;

%%
figure
varList = varA;
positionFiles =posA;
%%
%positionFiles([4,14,18])=[]
%varList([4,14,18])=[]
badIdx=[];
for iiFile =1:length(positionFiles)
traj = csvread(char(positionFiles(iiFile))); 
traj=traj(1:500,:);
%traj = csvread(strcat(char(positionFiles(iiFile))));%,'.csv')); 
%varList = varsList;

Ivars =load(char(varList(iiFile)),'Irot','cropImBot','cropImTop','th',...
    'saveFileNameTraj','iiCell','newTopFullSub','newBotFullSub','newPreAddFull',...
    'cellMask','cropImTop','cropImBot','fittedCellNoConv','rotCMask','minBot','maxBot',...
        'maskedCellFull');

minBot = Ivars.minBot;
maxBot = Ivars.maxBot;
nucDirPath=fullfile (strcat('Z:\dsmendez\diana\', fileparts(Ivars.saveFileNameTraj)),'nucleoid0') ;
memDirPath=fullfile (strcat('Z:\dsmendez\diana\', fileparts(Ivars.saveFileNameTraj)),'membrane0') ;

Itop = Ivars.cropImTop;
Itop = imrotate(Itop,-Ivars.th,'bicubic');
Ibot =Ivars.cropImBot;
Ibot = imrotate(Ibot,-Ivars.th,'bicubic');
topFullSub =Ivars.newTopFullSub;
botFullSub =Ivars.newBotFullSub;

newTraj = traj;

figure; imshow(Ivars.maskedCellFull(:,:,1),[])
for iiFrame =1:size(Ivars.maskedCellFull,3)
    tempTop=Itop(:,:,iiFrame).*Ivars.rotCMask;
    tempTop(tempTop==0)=nan;
    avgIntTop(iiFrame )=nanmean(tempTop(:)); 
    tempBot=Ibot(:,:,iiFrame).*Ivars.rotCMask;
    tempBot(tempBot==0)=nan;
    avgIntBot(iiFrame) =nanmean(tempBot(:)); 
end

avgTop = mean(avgIntTop);
avgBot = mean(avgIntBot);


avgTopF(iiFile)=avgTop;
avgBotF(iiFile)=avgBot;

meanTop(iiFile) = mean(traj(:,4));
meanBot(iiFile) = mean(traj(:,5));
maxTop(iiFile) = max(traj(:,4)); 
maxBot(iiFile) = max(traj(:,5)); 



SNRmaxTop(iiFile) = maxTop(iiFile)/avgTop;
SNRmaxBot(iiFile) = maxBot(iiFile)/avgBot;
SNRmeanTop(iiFile) = meanTop(iiFile)/avgTop;
SNRmeanBot(iiFile) = meanBot(iiFile)/avgBot;



try
    hello=topFullSub(:,:,:,500)
catch
    for ii=size(topFullSub,4):500
        topFullSub(:,:,:,ii)=topFullSub(:,:,:,size(topFullSub,4)-1);
        botFullSub(:,:,:,ii)=botFullSub(:,:,:,size(botFullSub,4)-1);
   
    end
    
end

[posXYZrot,mult,logDiv,mald,mild]= getZFunction(newTraj,Itop, Ibot,topFullSub,botFullSub,Ivars.rotCMask,size(Ivars.maskedCellFull,3) );
%[posXYZrot,mult,logDiv,mald,mild]= getZFunctionOld(newTraj,Itop, Ibot,topFullSub,botFullSub,Ivars.rotCMask,size(Ivars.maskedCellFull,3),Ivars.newPreAddFull );
saveVars =posXYZrot;
try
    traj = posXYZrot ;
    traj= getRidOfNans(traj);
    posXYZrot=traj;
catch
    badIdx = [badIdx,iiFile];
end
   



mild ;
mald ;
maxLD(iiFile)=max(logDiv);
minLD(iiFile)=min(logDiv);
%posXYZrot = traj;
trajXVector =posXYZrot(:,1) *80;
trajYVector = posXYZrot(:,2)* 80;
trajZVector =posXYZrot(:,3)*100 ;
zeroX = (max(trajXVector)+min(trajXVector))/2;
zeroY = (max(trajYVector)+min(trajYVector))/2;
zeroZ = (max(trajZVector)+min(trajZVector))/2;
zeroZ=25.75*100;%(fitresult.b1+7.25)*100;
Xprime = trajXVector - zeroX; 
Yprime = trajYVector - zeroY;
Zprime = (trajZVector - zeroZ)/0.59;

%figure; plot(Xprime,Yprime); axis equal
figure; plot(Yprime,Zprime,'.r'); axis equal
title(strcat('logDiv max:',num2str(maxLD(iiFile)),'logDiv min',num2str(minLD(iiFile))))
%figure; plot(Xprime,Zprime); axis equal

%if iiFile==4
 %   saveFileGrainRot = fullfile(strcat(Ivars.saveFileNameTraj,num2str(Ivars.iiCell),'_fixedZFull3_2','.csv'));

%else
    saveFileGrainRot = fullfile(strcat(Ivars.saveFileNameTraj,num2str(Ivars.iiCell),'_fixedZFull3','.csv'));
%end
    dlmwrite(saveFileGrainRot,posXYZrot,'delimiter',',','precision',32);

% hold on 
% plot (logDiv)
% drawnow



end 
%%
%%
fileStruct =uipickfiles()
load(char(fileStruct))
%%
PEDESTAL =99;

clear intMem
clear intNuc
for iiNucImages = 6
    nucDirPath =char(fileStruct.nucleoidParent(iiNucImages));
    memDirPath =char(fileStruct.membraneParent(iiNucImages));
    matFilesNuc = dir(fullfile(nucDirPath,'*.tif'));
    matFilesMem = dir(fullfile(memDirPath,'*.tif'));
    for iiSlice =1:length(matFilesNuc)
        iiSlice
        nucFile=fullfile(matFilesNuc(iiSlice).folder,matFilesNuc(iiSlice).name)
        memFile=fullfile(matFilesMem(iiSlice).folder,matFilesMem(iiSlice).name)
        
        Inuc = imread(char(nucFile))-PEDESTAL;
        Imem = imread(char(memFile))-PEDESTAL;

        Imem = imtranslate(Imem,[3, 7]); %empirically seen translation to match displacement of blue channel. 
        for iiRow=1:7
            Imem(iiRow,:) = Imem(8,:);
        end
        for iiCol =1:3
            Imem(:,iiCol) = Imem(:,4);
        end

        intMem(iiSlice) = sum(sum(Imem));
        intNuc(iiSlice) = sum(sum(Inuc));
    end
end 

[intCorrMem,fitInt]= bleachCorrection(intMem');
[intCorrNuc,fitInt]= bleachCorrection(intNuc(1:end)');
figure; plot(intCorrMem)

figure; plot(intCorrNuc)

intCorrNucMin = intCorrNuc - min(intCorrNuc)
createGaussianFit1D(1:size(intCorrNucMin),intCorrNucMin)


%%
posXYZrot=newTraj;% csvread(strcat(saveFileGrainRot));
figure
for iiFrame = 1:200
    iiFrame

    imshow(maskedCellFull(:,:,iiFrame),[])
    hold on 
    plot(posXYZrot(iiFrame,1),posXYZrot(iiFrame,2),'.r')
    drawnow

    %maskedCellTop(:,:,iiFrame)=tempTop.*cellMask;
    %maskedCellBot(:,:,iiFrame)=tempBot.*cellMask;

%Irot = imrotate(maskedCellFull,-Ivars.th);
end


%%
for iiFile =1:length(goodTracks )
%traj = csvread(char(positionFiles(iiFile))); 
    
traj = csvread(strcat(char(goodTracks(iiFile))));%,'.csv')); 

figure
plot(traj(:,2),traj(:,3))
end 




