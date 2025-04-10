function [DE, C, estPatt, refPatt] = directivityerror(sourceCoeff, estimatedCoeff, sourceOr, ...
    thetaAx)
% function DIRECTIVITYERROR
% This function computes the directivity error and the correlation between 
% the estimated directivity pattern and the actual pattern of a sound 
% source.
% Parameters:                                                            
%   - sourceCoeff: the coefficients of the actual pattern
%   - estimatedCoeff: the coefficients estimated
%   - sourceOr: the looking angle of the source
%   - thetaAx: the angle ax in radians
%                                              
% Returns:
%   - DE: directivity error as the mean in frequency of the MSE of the
%   estimated radiance pattern
%   - C: the correlation matrix or vector between the reference pattern and
%   the estimated one

orderRef = size(sourceCoeff, 2)-1;
% Compute the reference pattern
[~, refPatt] = ch2pol(360, sourceCoeff, sourceOr);
refPatt = refPatt.^2;
% The order of the estimated pattern
orderEst = size(estimatedCoeff, 2)/2 - 1;
estPatt = (getchbasis(thetaAx, orderEst) * squeeze(estimatedCoeff).').^2;

% Directivity error
DE = mean(mean((repmat(refPatt,1,size(estPatt,2)) - estPatt).^2));
yDiff = refPatt - estPatt;
% Correlation 
C = 1 - std(yDiff);%corr(estPatt, refPatt, 'type', 'Kendall');

end