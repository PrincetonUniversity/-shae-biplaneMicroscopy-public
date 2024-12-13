function []=recalculatePositionsConvProb(files2Fix,trajOrProb,iiFile)
%  positionFiles = uipickfiles();
%  varList =uipickfiles();
load (files2Fix)
if (trajOrProb)
    positionFiles = trajFilesList;
else 
    positionFiles = trajFilesListProb;
end

traj = csvread(strcat(char(positionFiles(iiFile)),'.csv')); 
traj = traj(:,1:6);

%for iiFile =1:length(positionFiles)
%traj = csvread(char(positionFiles(iiFile))); 
    


Ivars =load(char(varsList(iiFile)),'Irot','cropImBot','cropImTop','th',...
    'saveFileNameTraj','iiCell','newTopFullSub','newBotFullSub','newPreAddFull',...
    'cellMask','cropImTop','cropImBot','fittedCellNoConv','rotCMask');
Ivars

figure;imshow(Ivars.rotCMask,[])
se1 = strel('line',1,0);
se2 = strel('line',1,90);
SE=strel('disk',1);
cellMask = imdilate(Ivars.rotCMask,[se1 se2],'full');
cellMask = imdilate(cellMask,SE);%[se1 se2],'full')
cellMask = imdilate(cellMask,[se1 se2],'full');
cellMask = imdilate(cellMask,SE);%[se1 se2],'full')
figure;imshow(cellMask,[])


%%
clear maskedCellFull

for iiFrame = 1:size(Ivars.cropImTop,3)
    iiFrame
    
    tempTop = Ivars.cropImTop(:,:,iiFrame);
    tempBot = Ivars.cropImBot(:,:,iiFrame);
    imageInNoConv= tempTop+tempBot;

    resMatNoConv = imageInNoConv -0.3*max(imageInNoConv(:))*Ivars.fittedCellNoConv; 
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
 for iiFrame =1:size(Ivars.cropImTop,3)
     iiFrame
     
    temp = traj (traj(:,6)==iiFrame,:)
    stFrame = iiFrame-1;
    while isempty(temp)
        if stFrame ==0
            newTraj=[0,0,0,0,0,1];
            temp =  traj (1,:)
            break;
        end
        stFrame
        temp = newTraj (newTraj(:,6)==stFrame,:)
        stFrame = stFrame -1;
    end
        
    filteredImage = Irot(:,:,iiFrame);
    
    temp
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
 
 frames =1:size(Ivars.cropImTop,3);
 thresh= quantile(dist,.95);
 missFrame = frames(dist>thresh)+1
 
 
 while (~isempty(missFrame))
    iiFrame = missFrame(1)
    temp = newTraj (newTraj(:,6)==iiFrame-1,:)
    
    filteredImage = Irot(:,:,iiFrame);
         
    [tempPeaks,~]= centfinderDVM (filteredImage, temp(1:2),9);
   
    if iiFrame <size(Ivars.cropImTop,3)-1
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

Itop = Ivars.cropImTop;
Itop = imrotate(Itop,-Ivars.th,'bicubic');
Ibot =Ivars.cropImBot;
Ibot = imrotate(Ibot,-Ivars.th,'bicubic');
topFullSub =Ivars.newTopFullSub;
botFullSub =Ivars.newBotFullSub;


posXYZrot= getZFunction(newTraj,Itop, Ibot,topFullSub,botFullSub,Ivars.rotCMask,size(Itop,3)  );


saveFileGrainRot = fullfile(strcat(Ivars.saveFileNameTraj,num2str(Ivars.iiCell),'_fixedZFullProb','.csv'))
dlmwrite(saveFileGrainRot,posXYZrot,'delimiter',',','precision',32);

%end 
 
end
 
% %%
% posXYZrot=newTraj;% csvread(strcat(saveFileGrainRot));
% figure
% for iiFrame = 1:200
%     iiFrame
% 
%     imshow(maskedCellFull(:,:,iiFrame),[])
%     hold on 
%     plot(posXYZrot(iiFrame,1),posXYZrot(iiFrame,2),'.r')
%     drawnow
% 
%     %maskedCellTop(:,:,iiFrame)=tempTop.*cellMask;
%     %maskedCellBot(:,:,iiFrame)=tempBot.*cellMask;
% 
% %Irot = imrotate(maskedCellFull,-Ivars.th);
% end