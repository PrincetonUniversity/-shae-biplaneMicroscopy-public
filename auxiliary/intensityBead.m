function [intTop,intBot,beadTop,beadBot]= intensityBead(imTop,imBot,peaksIn,sz)

for iiElement = 1:length (peaksIn(:,1))
   % try
    %figure; imshow(imTop,[])
    imTop = (imTop>0).*imTop;
    min(imTop);
    %figure; imshow(imTop,[])
    %iiElement
    rV = round(peaksIn(iiElement,2));
    cV =round( peaksIn(iiElement,1));
    if (rV-floor(sz/2)<1)
        minRv = 1;
    else 
        minRv = rV-floor(sz/2);
    end
    if (cV-floor(sz/2)<1)
        minCv =1;
    else
        minCv = cV-floor(sz/2);
    end

    if (rV+floor(sz/2)>size(imTop,1))
        maxRv = size(imTop,1);
    else 
        maxRv = rV+floor(sz/2);
    end
    if (cV+floor(sz/2)>size(imTop,2))
        maxCv = size(imTop,2);
    else 
        maxCv = cV+floor(sz/2);
    end
    tmpBeadTop = imTop (minRv:maxRv, minCv:maxCv);
    tmpBeadBot = imBot (minRv:maxRv, minCv:maxCv);

    intTop(iiElement)=sum(tmpBeadTop(:));
    intBot(iiElement)=sum(tmpBeadBot(:));

    beadTop(1:size(tmpBeadTop,1),1:size(tmpBeadTop,2),iiElement)=tmpBeadTop;
    beadBot(1:size(tmpBeadTop,1),1:size(tmpBeadTop,2),iiElement)=tmpBeadBot;
    
end