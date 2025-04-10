function [doaWeight, doa, meanPseudoSpectrum] = doaestimator(signal, H, thetaAx, mask)
%% DOAESTIMATOR
% This function compute the DOA for a circular microphone array based on
% the peak evaluation of the pseudospectrum obtained as the absolute value
% of the output of a beamformer
%   Params:
%       - signal: the microphone signals in frequency
%       - H: the beamformer filter dim [micN dirN freqN]
%       - thetaAx: the direction axis
%       - mask: binary mask of active frequncy bins
%   Returns:
%       - doaWeight: the doa value
%       - doa: the direction of arrival
%       - meanPseudoSpectrum: averaged pseudospectrum

if isempty(mask) == true
    mask = ones(size(signal));
end

% Apply the mask
signal(~mask) = 0;  

% Compute the pseudo-spectrum
aux2 = repelem(signal, 1,1,length(thetaAx));        % Repeat signal for 360
aux2 = permute(aux2,[2,3,1]);                       % Move indexes
pseudoSpectrum = abs(squeeze(sum(H .* aux2, 1)));   % OUTPUT
% Normalize the pseudo spectrum
pseudoSpectrum = pseudoSpectrum ./ max(pseudoSpectrum(:));

% Apply the mask
% pseudoSpectrum(:,~mask) = 0;

% Aveaging the pseudospectrum in frequency
% meanPseudoSpectrum = geomean(pseudoSpectrum, 2);
meanPseudoSpectrum = mean(pseudoSpectrum, 2);

repMeanPseudoSpectrum = repelem(meanPseudoSpectrum, 1,3);
doaAx = repelem(thetaAx', 1, 3);
% DOA estimation and pruning
% [doaWeight, doa] = findpeaks(meanPseudoSpectrum, thetaAx, 'SortStr', ...
%     'descend');

% Normalize
% repMeanPseudoSpectrum = repMeanPseudoSpectrum ./ max(repMeanPseudoSpectrum);
[doaWeight, doaIdx] = findpeaks(repMeanPseudoSpectrum(:), 'SortStr', ...
    'descend', 'MinPeakDistance', 5);
doa = doaAx(doaIdx);
doa = unique(doa, 'stable');
doaWeight = unique(doaWeight, 'stable');

end