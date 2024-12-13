 function stitched = stitchingTrajProbConvBleachCorr(traj)
    maxDisp = 300;
   %filename = strcat(saveFileGrainRot,'.csv');
   % traj = csvread(filename);
   counter =0;
    %%
    %for trouble shooting:

    size(traj);
%     figure
%     plot(traj(:,6),traj(:,7),'b.')
%     hold on
%     plot(traj(:,6),traj(:,8),'g.')
%     figure 
%     histogram(traj(:,7))
%     plot(traj(:,6),traj(:,7)+traj(:,8),'r.')
    sumTraj = traj(:,4)+traj(:,5);
    [intCorr,fitInt]= bleachCorrection(sumTraj);
    %size(resVals)
    size(traj(:,6))
   % [resCorr,fitRes]= bleachCorrection(resVals,traj(:,6));
    
    %ratioRes =resVals./sumTraj;
    %figure; plot(traj(:,6),ratioRes,'r.')
clear minRes
    for iiTimeStep = 1:max(traj(:,6))
        
        try
        temp = intCorr(traj(:,6)==iiTimeStep);
        tempIndx = find (traj(:,6)==iiTimeStep);
        
        tempMax=tempIndx(temp==max(temp),1);
        %size(temp)
        %min(temp)
        %max(temp)
        maxInt(iiTimeStep,:) = [iiTimeStep,max(temp),tempMax(1)];
        %minRes(iiTimeStep)=ratioRes(tempMax);
%         
%           tempR = ratioRes(traj(:,6)==iiTimeStep);
%           tempIndx = find (traj(:,6)==iiTimeStep);
%           tempMin=tempIndx(tempR==min(tempR),1);
          
          %maxInt(iiTimeStep,:) = [iiTimeStep,min(tempR),tempMin(1)];
        catch
            maxInt(iiTimeStep,:) = [iiTimeStep,NaN,NaN];
        end

    end
    %plot(maxInt(:,1), maxInt(:,2), 'r.')
    %meanRes = nanmean(ratioRes);
    %stdRes = nanstd(ratioRes);
    %rangeRes = [meanRes-1.5*stdRes,meanRes + 1.5*stdRes]
    
    perc =quantile(maxInt(:,2),0.8);
    indMax = maxInt(:,2)>perc;
    %maxInt(:,2)
    %perc =quantile(maxInt(:,2),0.05);
    %indMax = maxInt(:,2)<perc;
%     figure
%     plot (maxInt(indMax,1), maxInt(indMax,2) ,'*g')
%     figure
%     plot(traj(:,6),ratioRes,'.b')
%     hold on
%     plot(maxInt(indMax,1),minRes(indMax),'rp')

    startPts = traj(maxInt(indMax, 3), :);
    size(startPts);
    
    preStitched =zeros(traj(end,6),size(traj,2));
    
    for iiStep = 1:traj(end,6)
        if(any(startPts(:,6)==iiStep))
            preStitched(iiStep,:)=startPts((startPts(:,6)==iiStep),:);
        end
    end
 
    
    
    trajTemp=traj;
    %here we found the most likely peaks. 
    %Then we are gonna build trajectories from each of this points. 
    
    
    %Now the question is how to build the trajectories? 
    %for 
    %stitched = p
    %%
 
        
        topInt = traj(:,4);
        botInt = traj(:,5);
        stitched = zeros(trajTemp(end,6),size(traj,2));
        totalDist =0;
        totalElements =0;      
        


        if (startPts(1,6)~= 1)
            % stitch backwards
            %%
            iiFirst =0;
            while (startPts(1,6)-iiFirst~=1 )
                iiFirst=iiFirst+1;
                
                if (iiFirst == 1)
                    stitched(startPts(1,6),:)=startPts(1,:);      
                end            
                iiFrNum = startPts(1,6)-iiFirst;
                dataPts = trajTem p(trajTemp(:,6)==iiFrNum,:);
                if isempty(dataPts)
                    warning('no data point for frame %d', iiFrNum)
                else
                    counter=1;
                    while ~any(stitched(iiFrNum+counter,:))
                        counter =counter +1;
                    end
                    prevPts =stitched(iiFrNum+counter,:);    


                    if iiFrNum == trajTemp(1,6)
                        nextPts =dataPts;                  
                    else
                        counter =1;
                        nextPts = trajTemp(trajTemp(:,6)==iiFrNum-counter,:);

                        while (isempty(nextPts))        
                            counter = counter +1;
                            nextPts = trajTemp(trajTemp(:,6)==iiFrNum-counter,:);
                            if (iiFrNum -counter ==1)
                                nextPts =dataPts;
                            end
                        end
                    end

                    nextDs = eucDistance(dataPts(:,1:3),nextPts(:,1:3));
                    prevDs = eucDistance(prevPts(:,1:3),dataPts(:,1:3));

                    if size(dataPts,1)==1
                        if (prevDs+min(nextDs))<maxDisp
                            stitched(iiFrNum,:)=dataPts;
                            totalDist=totalDist+prevDs+min(nextDs);
                            totalElements = totalElements +1;
                        else
                             warning('Distance is %d, in frame %d', prevDs+min(nextDs), iiFrNum)
                        end
                    else

                        minGlobal = 10000;
                        for iiPoints = 1:size(dataPts,1)
                            minD= min(prevDs(iiPoints)+nextDs(:,iiPoints));
                            minGlobal =min(minGlobal, minD);
                            if minD == minGlobal
                                minIndex = iiPoints;
                            end
                        end
                        minGlobal;
                        if (minGlobal<maxDisp)
                            stitched(iiFrNum,:)=dataPts(minIndex,:);
                            totalDist=totalDist+minGlobal;
                            totalElements = totalElements +1;
                        else 

                            warning('minGlobal is %d, in frame %d', minGlobal, iiFrNum)
                        end
                    end                
                end
                 
            end
            
                
            firstPoint = stitched(1,:);
        else
            firstPoint =startPts (1,:);
        end
        
        

          
        
       stitched(1,:)= firstPoint;
       

        for iiFrame = startPts(1,6)+1:trajTemp(end,6)   
        
            if (any(preStitched(:,6)==iiFrame))
                5;
                stitched (iiFrame,:) = preStitched(preStitched(:,6)==iiFrame,:);
            else
                 dataPts = trajTemp(trajTemp(:,6)==iiFrame,:);
                if isempty(dataPts)
                    warning('no data point for frame %d', iiFrame)
                else
                    prevPts =stitched(iiFrame-1,:);
                    counter=2;

                    while prevPts==zeros(1,size(traj,2))
                        prevPts = stitched(iiFrame-counter,:);
                        counter = counter +1;
                        6;
                    end
                    
                    if iiFrame == trajTemp(end,6)
                        nextPts =dataPts;
                    elseif any(preStitched(:,6)==iiFrame+1)
                        nextPts=preStitched(preStitched(:,6)==iiFrame+1,:) ;
                    else
                        nextPts = trajTemp(trajTemp(:,6)==iiFrame+1,:);
                    end
                    counter=2;
                    while isempty(nextPts)
                        nextPts = trajTemp(trajTemp(:,6)==iiFrame+counter,:);
                        counter = counter +1;
                        7;
                    end
                    
                    nextDs = eucDistance(dataPts(:,1:3),nextPts(:,1:3));
                    prevDs = eucDistance(prevPts(:,1:3),dataPts(:,1:3));
                    
                     if size(dataPts,1)==1
                        if (prevDs+min(nextDs))<maxDisp
                            stitched(iiFrame,:)=dataPts;
                            totalDist=totalDist+prevDs+min(nextDs);
                            totalElements = totalElements +1;
                        else
                             warning('Single Distance is %d, in frame %d', prevDs+min(nextDs), iiFrame)
                        end
                    else

                        minGlobal = 10000;
                        for iiPoints = 1:size(dataPts,1)
                            minD= min(prevDs(iiPoints)+nextDs(:,iiPoints));
                            minGlobal =min(minGlobal, minD);
                            if minD == minGlobal
                                minIndex = iiPoints;
                            end
                        end
                        minGlobal;
                        if (minGlobal<maxDisp)
                            stitched(iiFrame,:)=dataPts(minIndex,:);
                            totalDist=totalDist+minGlobal;
                            totalElements = totalElements +1;
                        else 
                            maxDisp2 = maxDisp +300;
                            if (minGlobal<maxDisp2)
                                stitched(iiFrame,:)=dataPts(minIndex,:);
                                totalDist=totalDist+minGlobal;
                                totalElements = totalElements +1;
                            else 
%                                 maxDisp3 = maxDisp2 +300;
%                                  if (minGlobal<maxDisp3)
%                                     stitched(iiFrame,:)=dataPts(minIndex,:);
%                                     totalDist=totalDist+minGlobal;
%                                     totalElements = totalElements +1;
%                                  else 
                                      warning('minGlobal is %d, in frame %d', minGlobal, iiFrame)
                                      counter =counter +1;
                                 %end
                            end
                            
                        end
                     end                 
                end
            end
            
            %index =totalDist./totalElements==min(totalDist./totalElements);
            %stitched = stitched(:,:,index);
        end
         

counter 
    
    %%
    %for iiCell = traj(1,5):traj(end,5)
%         trajTemp = traj(traj(:,5)==iiCell,:) ;
%         firstPoint =trajTemp(trajTemp(:,6)==1,:)
%         topInt = traj(:,7);
%         botInt = traj(:,8);
%         
%         
%         
%         if(isempty(firstPoint))
%             firstPoint =trajTemp(trajTemp(:,6)==trajTemp(1,6),:);
%         end
%         stitched = zeros(trajTemp(end,6),size(traj,2),size(firstPoint,1));
%         totalDist =zeros(size(firstPoint,1),1);
%         totalElements =zeros(size(firstPoint,1),1);
%        
%             stitched(1,:,iiStart)= firstPoint(iiStart,:);
%          
% 
%             for iiFrame = 2:trajTemp(end,6)   
%                 dataPts = trajTemp(trajTemp(:,6)==iiFrame,:);
%                 if isempty(dataPts)
%                     warning('no data point for frame %d', iiFrame)
%                 else
%                     prevPts =stitched(iiFrame-1,:,iiStart);
%                     counter=2;
%                     while prevPts==zeros(1,size(traj,2))
%                         prevPts = stitched(iiFrame-counter,:,iiStart);
%                         counter = counter +1;
%                     end
%                     if iiFrame == trajTemp(end,6)
%                         nextPts =dataPts;
%                     else
%                         nextPts = trajTemp(trajTemp(:,6)==iiFrame+1,:);
%                     end
%                     counter=2;
%                     while isempty(nextPts)
%                         nextPts = trajTemp(trajTemp(:,6)==iiFrame+counter,:);
%                         counter = counter +1;
%                     end
%                     nextDs = eucDistance(dataPts(:,1:3),nextPts(:,1:3));
%                     prevDs = eucDistance(prevPts(:,1:3),dataPts(:,1:3));
% 
%                     if size(dataPts,1)==1
%                         if (prevDs+min(nextDs))<maxDisp
%                             stitched(iiFrame,:,iiStart)=dataPts;
%                             totalDist(iiStart)=totalDist(iiStart)+prevDs+min(nextDs);
%                             totalElements(iiStart) = totalElements(iiStart) +1;
%                         else
%                              warning('Distance is %d, in frame %d', prevDs+min(nextDs), iiFrame)
%                         end
%                     else
% 
%                         minGlobal = 10000;
%                         for iiPoints = 1:size(dataPts,1)
%                             minD= min(prevDs(iiPoints)+nextDs(:,iiPoints));
%                             minGlobal =min(minGlobal, minD);
%                             if minD == minGlobal
%                                 minIndex = iiPoints;
%                             end
%                         end
%                         minGlobal;
%                         if (minGlobal<maxDisp)
%                             stitched(iiFrame,:,iiStart)=dataPts(minIndex,:);
%                             totalDist(iiStart)=totalDist(iiStart)+minGlobal;
%                             totalElements(iiStart) = totalElements(iiStart) +1;
%                         else 
% 
%                             warning('minGlobal is %d, in frame %d', minGlobal, iiFrame)
%                         end
%                     end
% 
%                 end
%             
%         
%        end
%         index =totalDist./totalElements==min(totalDist./totalElements);
%         stitched = stitched(:,:,index);
%     end

end 