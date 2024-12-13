function [posXYZrot,multFull,logDivFull,maxLD,minLD]= getZFunction(traj,Itop, Ibot,topFullSub,botFullSub,rotCMask,sizeMask )
%%
index=1;
clear posXYZ

for iiFrame =1:sizeMask
    tempTop=Itop(:,:,iiFrame).*rotCMask;
    tempTop(tempTop==0)=nan;
    avgIntTop(iiFrame )=nanmean(tempTop(:)); 
    tempBot=Ibot(:,:,iiFrame).*rotCMask;
    tempBot(tempBot==0)=nan;
    avgIntBot(iiFrame) =nanmean(tempBot(:)); 
end

avgTop = mean(avgIntTop);
avgBot = mean(avgIntBot);


maxTop = max(traj(:,4)); 
maxBot = max(traj(:,5)); 
SNRmaxTop = maxTop/avgTop;
SNRmaxBot = maxBot/avgBot;

    
%maxLD= 0.0085*SNRmaxBot-0.49; %Empirically found Vuldi
%minLD= -0.008*SNRmaxTop +0.99; %Empirically found Vuldi
    
%maxLD = 0.65%A
%minLD= -0.0071*SNRmaxTop +0.97; %Empirically found A

% maxLD = 0.0092*SNRmaxBot-0.62; %Empirically found PFV
% minLD =-0.0031*SNRmaxTop +0.3;%-0.1;

% maxLD =0.9% -18
% minLD =0.0011*SNRmaxTop -0.14;%-18

% maxLD =0.0056*SNRmaxBot-0.075% +15
% minLD =-0.0035*SNRmaxTop +0.18;%+15

 maxLD =0.0048*SNRmaxBot-0.1;% vuldi exp
 minLD =-0.0076*SNRmaxTop +0.84;%vuldi exp
%%
for iiFrame =1:size(traj,1)
    size(traj,1)
    iiFrame

    topIm = Itop(:,:,iiFrame);
    botIm = Ibot(:,:,iiFrame);
    

    
%  try
tempTop = topFullSub(:,:,:,iiFrame);
tempBot = botFullSub(:,:,:,iiFrame);
%tempPreAdd = preAddFull(:,:,:,iiFrame);

roughXY =[traj(:,1:5),(1:size(traj,1))'];
    posX = (roughXY(roughXY(:,6)==iiFrame,1));
    posY = (roughXY(roughXY(:,6)==iiFrame,2));
    rndposX = round(roughXY(roughXY(:,6)==iiFrame,1));
    time= round(roughXY(roughXY(:,6)==iiFrame,6))
    
    lastFrame =iiFrame;
    while rndposX==0
        lastFrame =lastFrame -1;
         posX = (roughXY(roughXY(:,6)==lastFrame,1));
        posY = (roughXY(roughXY(:,6)==lastFrame,2));
        rndposX = round(roughXY(roughXY(:,6)==lastFrame,1));
        time= lastFrame;
    end

    sizeBox =9;
    intTop=calcInt(topIm,[posX,posY],sizeBox);
    intBot=calcInt(botIm,[posX,posY],sizeBox);
    
    
    %voll5Top = calcVol(topIm/maxTop,[posX,posY])/intTop;
    %voll5Bot = calcVol(botIm/maxBot,[posX,posY])/intBot;
    voll5Top = calcVol(topIm,[posX,posY])/intTop;
    voll5Bot = calcVol(botIm,[posX,posY])/intBot;
    logDiv= log(voll5Bot/voll5Top)
    logDivFull(iiFrame) = logDiv;%/(maxTop+maxBot)*2
    

    m=(0.8/(maxLD-minLD));
    b=0.9-m*maxLD;
    mult = m*logDiv+b;
    multFull(iiFrame)=mult;    
    
    
%     %mult = (1/1.2) * logDiv +(1/12);% (-0.1 maps to 0 and 1.1 maps to 1)PFV!
%     %mult = (1/1.6)* logDiv+0.3125;vuldi
%     %mult = 1.25 * logDiv +0.125% (-0.1 maps to 0 and 1.1 maps to 1) AqLS!
%     mult = 1.43*logDiv +0.143 %aqls!
%     mult = 1/1.05*logDiv +8/21;
    % m=0%(1/0.95)
     %b=1-m*0.65;
     %mult = m*logDiv+b
%     multFull(iiFrame)=mult;

counterBig =0;
counterSmall=0;

    if mult > 1
        counterBig = counterBig+1;
        mult =1;
    end
    if mult < 0
        mult =0;
        counterSmall=counterSmall+1;
    end

    clear zCom
    clear varZcom
    clear posZ
    
    tempAdd = (1-mult)*tempTop + mult*tempBot;

    
    for iiElement=1%:size(rndposX)
        try
        iiElement
        tempC = tempAdd(:,:,iiElement);  
        tempC= double(tempC>0).*tempC;
        tempT= tempTop(:,:,iiElement); 
        tempB= tempBot(:,:,iiElement); 
        sliverX = tempC(:,rndposX(iiElement));

         sliverT = tempT(:,rndposX(iiElement));
         sliverB = tempB(:,rndposX(iiElement));
        [~,midx]=max(sliverX);
         [~,midT]=max(sliverT);
         [~,midB]=max(sliverB);
        [~,midx]=max(sliverX);

        
        cropSliver = sliverX(midx-5:midx+5);
        zIdxPos = 1:size(tempC,1);
        cropZidx = zIdxPos(midx-5:midx+5);
        cropSliverT = sliverT(midT-5:midT+5);
        zIdxPosT = 1:size(tempT,1);
        cropZidxT = zIdxPosT(midT-5:midT+5);
        cropSliverB = sliverB(midB-5:midB+5);
        zIdxPosB = 1:size(tempB,1);
        cropZidxB = zIdxPosB(midB-5:midB+5);
%     
        
        %[fitG,gofG,x1,y1]=createGaussianFit1D(cropZidx, cropSliver);
        [fitGT,gofGT,x1T,y1T]=createGaussianFit1D(cropZidxT, cropSliverT);
        [fitGB,gofGB,x1B,y1B]=createGaussianFit1D(cropZidxB, cropSliverB);
        [fitG,gofG,x1,y1]=createGaussianFit1D(cropZidx, cropSliver);
        
         a=fitGT.b1*(1-mult);
         b=fitGB.b1*mult;
         %if fitG.b1 <20 | fitG.b1>30
          %  posZ(iiElement) =NaN;
         %else
             posZ(iiElement)= a+b
         %end
         posT(iiElement)=fitGT.b1
         posB(iiElement)=fitGB.b1;
         posG(iiElement)=fitG.b1;
       if iiFrame ==1

         %figure
         subplot (1,3,1)
         imshow(tempT*(1-mult),[0,max(tempC(:))])
         hold on
         plot(rndposX(iiElement),fitGT.b1,'g*')
         title(strcat('Top',num2str(fitGT.b1)))

         subplot (1,3,2)
         imshow(tempB*mult,[0,max(tempC(:))])
         hold on
         plot(rndposX(iiElement),fitGB.b1,'b*')
         title(strcat('Bot',num2str(fitGB.b1)))


         subplot (1,3,3)
         imshow(tempC,[0,max(tempC(:))])
         title(strcat('w_sum',num2str(mult),',',num2str(a+b)))

         hold on
         plot(rndposX(iiElement),posZ(iiElement),'r*')
         plot(rndposX(iiElement),fitGT.b1,'g*')
         plot(rndposX(iiElement),fitGB.b1,'b*')
         
         drawnow
       end
        
         catch
             warning(strcat('frame',num2str(iiFrame) ,' did not work'))
             posZ(iiElement)=NaN;
             posT=NaN;
             posB=NaN;
             posG=NaN;
         end
    end


        posXYZrot(index:index+size(rndposX,1)-1,:)=[posX,posY,posZ',intTop,intBot,time,logDiv,posT',posB',posG,mult];

        index = index + size(rndposX,1);

end

posXYZrot
