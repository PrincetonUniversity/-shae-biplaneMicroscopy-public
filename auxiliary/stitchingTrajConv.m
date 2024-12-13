 function stitched = stitchingTrajConv(traj)
    maxDisp = 500;
   %filename = strcat(saveFileGrainRot,'.csv');
   % traj = csvread(filename);
    %%
    size(traj)
    for iiCell =1% traj(1,3):traj(end,3)
        trajTemp = traj;
        firstPoint =trajTemp(trajTemp(:,6)==1,:)
        
        
        if(isempty(firstPoint))
            firstPoint =trajTemp(trajTemp(:,6)==trajTemp(1,6),:);
        end
        stitched = zeros(trajTemp(end,6),size(traj,2),size(firstPoint,1));
        totalDist =zeros(size(firstPoint,1),1);
        totalElements =zeros(size(firstPoint,1),1);
       for iiStart = 1:size(firstPoint,1) 
            stitched(1,:,iiStart)= firstPoint(iiStart,:);
         

            for iiFrame = 2:trajTemp(end,6)   
                %iiFrame
                dataPts = trajTemp(trajTemp(:,6)==iiFrame,:);
                if isempty(dataPts)
                    warning('no data point for frame %d', iiFrame)
                else
                    prevPts =stitched(iiFrame-1,:,iiStart);
                    counter=2;
                    while prevPts==zeros(1,size(traj,2))
                        prevPts = stitched(iiFrame-counter,:,iiStart);
                        counter = counter +1;
                    end
                    if iiFrame == trajTemp(end,6)
                        nextPts =dataPts;
                    else
                        nextPts = trajTemp(trajTemp(:,6)==iiFrame+1,:);
                    end
                    counter=2;
                    while isempty(nextPts)
                        nextPts = trajTemp(trajTemp(:,6)==iiFrame+counter,:);
                        counter = counter +1;
                    end
                    nextDs = eucDistance(dataPts(:,1:3),nextPts(:,1:3));
                    prevDs = eucDistance(prevPts(:,1:3),dataPts(:,1:3));

                    %nextDs = eucDistance(dataPts(:,1:2),nextPts(:,1:2));
                    %prevDs = eucDistance(prevPts(:,1:2),dataPts(:,1:2));
                    if size(dataPts,1)==1
                        if (prevDs+min(nextDs))<maxDisp
                            stitched(iiFrame,:,iiStart)=dataPts;
                            totalDist(iiStart)=totalDist(iiStart)+prevDs+min(nextDs);
                            totalElements(iiStart) = totalElements(iiStart) +1;
                        else
                             warning('Distance is %d, in frame %d', prevDs+min(nextDs), iiFrame)
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
                            stitched(iiFrame,:,iiStart)=dataPts(minIndex,:);
                            totalDist(iiStart)=totalDist(iiStart)+minGlobal;
                            totalElements(iiStart) = totalElements(iiStart) +1;
                        else 
                            if (minGlobal<maxDisp+500)
                                stitched(iiFrame,:,iiStart)=dataPts(minIndex,:);
                                totalDist(iiStart)=totalDist(iiStart)+minGlobal;
                                totalElements(iiStart) = totalElements(iiStart) +1;
                            end
                            warning('minGlobal is %d, in frame %d', minGlobal, iiFrame)
                        end
                    end

                end
            end
        
       end
        index =totalDist./totalElements==min(totalDist./totalElements);
        stitched = stitched(:,:,index);
    end

end 