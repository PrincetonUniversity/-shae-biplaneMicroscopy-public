function [peaksOut,resPeaks] = centfinderDVM (imageIn, peaksIn,sz)

%%
peaksIn
for iiElement = 1:length (peaksIn(:,1))
   % try
    %figure; imshow(imageIn,[])
    imageIn = (imageIn>0).*imageIn;
    min(imageIn);
    %figure; imshow(imageIn,[])
    %iiElement
    rV =round( peaksIn(iiElement,2));
    cV = round(peaksIn(iiElement,1));
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

    if (rV+floor(sz/2)>size(imageIn,1))
        maxRv = size(imageIn,1);
    else 
        maxRv = rV+floor(sz/2);
    end
    if (cV+floor(sz/2)>size(imageIn,2))
        maxCv = size(imageIn,2);
    else 
        maxCv = cV+floor(sz/2);
    end

    
%         figure
%     imshow(imageIn,[])

    tmpBead = imageIn (minRv:maxRv, minCv:maxCv);
    tmpBead = tmpBead/max(tmpBead(:));
%     figure
%     imshow(tmpBead,[])

    
    xCoordMat =zeros(size(tmpBead,1),size(tmpBead,2))+[1:size(tmpBead,2)];                                                                                                                                                                                                                                                                                                     
    yCoordMat =zeros(size(tmpBead,1),size(tmpBead,2))+[1:size(tmpBead,1)]';
    xavg(iiElement)=sum(sum(tmpBead.*xCoordMat))/sum(tmpBead(:));
    yavg(iiElement)=sum(sum(tmpBead.*yCoordMat))/sum(tmpBead(:));
    posX = xCoordMat-xavg(iiElement);
    posY=yCoordMat-yavg(iiElement);
    wX = posX .*tmpBead;
    wY =posY.*tmpBead;
    rg = (sum(sum((posX.^2+posY.^2).*tmpBead))/sum(tmpBead(:)));     

    % using the centroid information, fit to a gauss2d
    b = [min(tmpBead,[],2)',min(tmpBead,[],1)];
    a0(1) = mean(b(:));
    a0(2) = tmpBead(floor(yavg(iiElement)),floor(xavg(iiElement)));
    a0(3) = xavg(iiElement) ;  % maybe use the center of the box 'r' here
    a0(4) = yavg(iiElement);   % maybe use the center of the box 'r' here
    a0(5) = sqrt(rg);
    a0(6) = sqrt(rg);
    a0(7) = 0;
a0
% Do Fitting
options = optimset('Display','none','Jacobian','off',...
'MaxFunEvals',3000,'MaxIter',1000);
[Param,resnorm,residual,exitflag,~]=lsqcurvefit(@gauss2D,a0,tmpBead,tmpBead,[],[],options);
  Param;
exitflag;

if exitflag < 1
    % if fitting fails output centroid method centers
    3;
    Param = zeros(7,1);
    xoffset = xavg(iiElement);
    yoffset = yavg(iiElement) ;
    resnorm=NaN;
    residual=NaN;

    else
    4;
    xoffset = Param(3);
    yoffset = Param(4);
end

% figure 
% imshow(tmpBead,[])
% hold on
% plot(xoffset,yoffset,'*r')

xoffset = xoffset - ceil(size(tmpBead,1)/2)+peaksIn(iiElement,1);
yoffset = yoffset - ceil(size(tmpBead,2)/2)+peaksIn(iiElement,2);

if (xoffset<minCv | yoffset< minRv | xoffset> maxCv | yoffset>maxRv)
    
    peaksOut(iiElement,:) = [xavg(iiElement)- ceil(size(tmpBead,1)/2)+peaksIn(iiElement,1),...
        yavg(iiElement)- ceil(size(tmpBead,2)/2)+peaksIn(iiElement,2)];
else
    peaksOut(iiElement,:) = [xoffset,yoffset];
end
size(peaksOut);
size(resnorm);
size(residual);
resPeaks (iiElement,:)=resnorm;%,residual]

% figure 
% imshow(imageIn,[])
% hold on
% plot(xoffset,yoffset,'*r')
% plot(xavg- ceil(size(tmpBead,1)/2)+peaksIn(iiElement,1),yavg- ceil(size(tmpBead,2)/2)+peaksIn(iiElement,2),'b*')


   % catch
     %   warning('something went wrong')
   % end
end








