function [fitresult, gof,xData,yData] = createGaussianFit1D(cropZidx, cropSliver)
%CREATEFIT(CROPZIDX,CROPSLIVER)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : cropZidx
%      Y Output: cropSliver
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 20-Feb-2022 02:09:03


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( cropZidx, cropSliver );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [1.25668573286826 28 4.11039315057867];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

%Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'cropSliver vs. cropZidx', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel cropZidx
% ylabel cropSliver
% grid on
% drawnow

