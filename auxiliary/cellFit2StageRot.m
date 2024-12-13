function [cellMat2,a1] = cellFit2StageRot (imageIn,a0)

% a0(1) = 37; %3um
% a0(2) = 6; %1um
% a0(3) = size(imageIn,2)/2;
% a0(4) = size(imageIn,1)/2;
% a0(5) = 0.9;
    
cellParams = fitCellScriptRot (a0,imageIn) ;
cellMat =fitCellonlyRot(cellParams,imageIn);
figure;imshow(cellMat,[])

a1(1)=cellParams(1)-2*cellParams(2)-2;
a1(2)=cellParams(2);
a1(3)=cellParams(3);
a1(4)=cellParams(4);
a1(5)=cellParams(5);
a1(6)=90;
a1(7)=cellParams(6);
theta =cellParams(6);

cellMat2 = fitCellonlyEndCapsRot(a1,imageIn);
figure;imshow(cellMat2,[])

C = imfuse(cellMat2,imageIn/max(imageIn (:)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
figure; imshow(C)


%figure;imshow(imfuse(cellMat,cellMat2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]))

end








