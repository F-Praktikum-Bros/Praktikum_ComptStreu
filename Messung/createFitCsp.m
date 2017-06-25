function [fitresult, gof] = createFitCsp(EpVec, Csp, errCsp)
%CREATEFIT(EPVEC,CSP,ERRCSP)
%  Create a fit.
%
%  Data for 'fitCsp' fit:
%      X Input : EpVec
%      Y Output: Csp
%      Weights : errCsp
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 25-May-2017 12:39:01


%% Fit: 'fitCsp'.
[xData, yData, weights] = prepareCurveData( EpVec, Csp, errCsp );

% Set up fittype and options.
ft = fittype( 'sqrt(2/pi)*a/w *exp(-((x-x0)/w)^2)+b*x+c+d/x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = [0.980214181846737 0.701608546221584 0.557442985873489 0.795199901137063 0.186872604554379 0.489764395788231];
opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'fitCsp' );
% h = plot( fitresult, xData, yData );
% legend( h, 'Csp vs. EpVec with errCsp', 'fitCsp', 'Location', 'NorthEast' );
% % Label axes
% xlabel('Energie [MeV]')
% ylabel('Anzahl Ereignisse pro Sekunde [1/s]')
% grid off