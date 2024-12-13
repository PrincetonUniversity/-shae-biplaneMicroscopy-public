function [cellMat2,theta,x0,y0,cellFit]= cellFit2StageNoConvRot (imageIn,a0)

% a0(1) = 37; %3um
% a0(2) = 6; %1um
% a0(3) = size(imageIn,2)/2;
% a0(4) = size(imageIn,1)/2;
% a0(5) = 0.9;
 figure;imshow(fitCellonlyNoConvRot(a0(1:5),{imageIn,a0(6)}))   
cellParams = fitCellScriptNoConvRot (a0,imageIn) ;

cellMat =fitCellonlyNoConvRot(cellParams(1:5),{imageIn,cellParams(6)});
C = imfuse(cellMat,imageIn/max(imageIn (:)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
figure; imshow(C)

figure;imshow(cellMat,[])
x0 =cellParams(3);5
a1(1)=cellParams(1)%-2*cellParams(2);
a1(2)=cellParams(2);
a1(3)=cellParams(3);
a1(4)=cellParams(4);
a1(5)=cellParams(5);
a1(6)=90;
a1(7)=cellParams(6);
theta =cellParams(6);


y0 =cellParams(4);
[cellMat2,cellFit] = fitCellonlyEndCapsNoConvRot(a1,imageIn);
%figure;imshow(cellMat2,[])
6
C = imfuse(cellMat2,imageIn/max(imageIn (:)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
figure; imshow(C)
7
end








