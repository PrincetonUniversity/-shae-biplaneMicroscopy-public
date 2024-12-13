function [cellMat,Param,fitcell1,fitcell2]= cellFitZConv (imageIn,a0)

%figure; imshow(imageIn,[0,max(imageIn(:))])
%title('Imagein')
%figure; imshow(fitZConvData(a0,imageIn),[])
%title('Initial guess')
%

% Do Fitting
%,'Algorithm','levenberg-marquardt',
options = optimset('Display','none','Jacobian','off',...
'MaxFunEvals',3000,'MaxIter',1000);
%options = optimoptions('lsqcurvefit','Display','none',...
%'MaxFunEvals',3000,'MaxIter',1000);
[Param,resnorm,residual,exitflag,~]=lsqcurvefit(@fitZConvData,a0,imageIn,imageIn/max(imageIn(:)),[],[],options);

[fitCell,fitcell1,fitcell2] =fitZConvData(Param,imageIn);
%figure; imshow(abs(residual),[0,1])
%title('Residuals')
%figure; imshow(fitCell,[0,1])
%title('fit')

cellMat =fitCell;
  
end








