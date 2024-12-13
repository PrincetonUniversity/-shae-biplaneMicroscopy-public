function [peaks,int]=pkfinderDVM(imageIn, threshold, sz)

%Locate peaks above threshold
%%
% imageIn =temp
% threshold =max(temp)*0.5
% sz= 15
imageIn(imageIn>threshold);
tempPeaks= imageIn>threshold;
intVals =imageIn(boolean(tempPeaks));
intVals=sort(intVals,'descend');

iiVal = 1;
clear rV;clear cV

while (~isempty(intVals))
        iiVal;

        [tempR,tempC]=find(imageIn.*tempPeaks==intVals(1));
        
        for iiCoords=1:length(tempC)

            rV(iiVal)= tempR(iiCoords);
            cV(iiVal)=tempC(iiCoords);
            mask=ones(size(imageIn));
            if (rV(iiVal)-floor(sz/2)<1)
                minRv = 1;
            else 
                minRv = rV(iiVal)-floor(sz/2);
            end
            if (cV(iiVal)-floor(sz/2)<1)
                minCv =1;
            else
                minCv = cV(iiVal)-floor(sz/2);
            end

            if (rV(iiVal)+floor(sz/2)>size(imageIn,1))
                maxRv = size(imageIn,1);
            else 
                maxRv = rV(iiVal)+floor(sz/2);
            end
            if (cV(iiVal)+floor(sz/2)>size(imageIn,2))
                maxCv = size(imageIn,2);
            else 
                maxCv = cV(iiVal)+floor(sz/2);
            end

            mask(minRv:maxRv,minCv:maxCv)=0;

            tempPeaks= tempPeaks.*mask();
            intVals =imageIn(boolean(tempPeaks));
            intVals=sort(intVals,'descend');
            iiVal =iiVal+1;
        end
    end

    peaks =[cV(:),rV(:)];

    for iiPeak =1:length(cV)
        int(iiPeak) = imageIn(rV(iiPeak),cV(iiPeak));
    end

