function [completeEstimate, directEstimate] = parametricvirtualmiking(array, source, cptPts, sphParams, params, macro)
%% parametricvirtualmiking
% This function performs virtual miking using a parametric technique.
% Parameters:
%   array: struct containing information on the distributed arrays. 
%   source: struct containing information on the sources. 
%   cptPts: struct containing information on the VMs.
%   sphParams: 
% 
% Mirco Pezzoli, 2021
% 22/06/2021
% v. 0.1
%
% Silvio Attoli
% 23/02/2025
% v. 0.2
% Adds Habets algorithm for generating vm diffus component, and DRR weights

fprintf("parametric vm method...\n");
%% Localize the sources
sourcePos = cell2mat(source.position.');
if macro.LOCALIZATION_PRS == true
    localizationParams.dereverb = true;
    localizationParams.secPerPosEstimation = -1;
    localizationParams.stepAngle = 1;
    localizationParams.plotPRS = false;
    localizationParams.localizationTest = 1000;
    [medianPosition, bestPosition, medianError, bestError] = ...
        sourcelocalization(array, source, localizationParams, params, macro);
    
else
    bestPosition = sourcePos(:,1:2);
    medianPosition = bestPosition;
    bestError = 0;
    medianError = 0;
end
fprintf('Best estimated source position(s): \n')
disp(bestPosition)
fprintf('Reference source position(s): \n')
disp(sourcePos)
fprintf('Best localization Error: \n')
disp(bestError)
fprintf('Median estimated source position(s): \n')
disp(medianPosition)
fprintf('Median localization Error: \n')
disp(medianError)

source.bestPosition = bestPosition;
source.medianPosition = medianPosition;


% Remove reverberation from the microphone signals
dereverbParams.sourcePosition = 'median';
dereverbParams.cdrMicN = array.micN;
array.cdrMicN = dereverbParams.cdrMicN;

[meanDerev, meanDerevSTFT, meanDiffuse, meanDiffuseSTFT] = dereverbarraywithdoa(array, source, ...
    dereverbParams, params, macro);

array.meanDerev = meanDerev;
array.meanDerevSTFT = meanDerevSTFT;
array.meanDiffuse = meanDiffuse;
array.meanDiffuseSTFT = meanDiffuseSTFT;
 
%% % estimate DRR of microphone arrays for the best weight assignment
mean3 = @(x) sqrt(mean(abs(x).^2,3));
totalDirectSTFT = cellfun(mean3, array.meanDerevSTFT, ...
    'UniformOutput', false); % mean for each array
totalDiffuseSTFT = cellfun(mean3, array.meanDiffuseSTFT, ...
    'UniformOutput', false);

DRRs_array = zeros(1, array.N);
for aa = 1:array.N
    time_direct_array = my_istft(totalDirectSTFT{aa}, params.analysisWin, ...
               params.synthesisWin, params.hop, params.Nfft, params.Fs);
    time_diffuse_array = my_istft(totalDiffuseSTFT{aa}, params.analysisWin, ...
               params.synthesisWin, params.hop, params.Nfft, params.Fs);
    
    power_direct_array = mean(abs(fft(time_direct_array)).^2, 2);
    power_diffuse_array = mean(abs(fft(time_diffuse_array)).^2, 2);
    %disp(size(power_direct_array));
    %disp(size(power_diffuse_array));
    DRRs_array(1, aa) = power_direct_array./power_diffuse_array;
end

% DRR of arrays
fig = figure('Visible', 'off');
plot(db(DRRs_array))
xlabel('VM index');
ylabel('[dB]');
title('DRR ESTIMATE of arrays')
grid on;
saveas(fig, 'results/DRR_of_arrays_estimate.png');
close(fig);

%% Spherical harmonics expansion estimation
if isempty(sphParams)
    sphParams.arraySignal = 'estimate';
    sphParams.sourcePosition = 'median';
    sphParams.maxOrder = 1;
    sphParams.cdrMicN = array.cdrMicN;
    sphParams.regParam.method = 'tikhonov';
    sphParams.regParam.nCond = 35;
    sphParams.type = 2;
    
end
if source.N > 1
    sphParams.maxOrder = sphParams.maxOrder*ones(1:source.N);  % Spherical harmonics max order
end

hCoeff = sphericalharmonicsestimation(array, source, sphParams, ...
    params, macro);


%% Estimate the direct signal using the sph expansion

[directEstimateSTFT, directSourceEstimate, arrayEstimateSTFT] = estimatedirectsignal(cptPts, hCoeff, array, source, ...
    sphParams, params, macro);



for mm = 1:cptPts.N
    directEstimate(:,mm) = my_istft(directEstimateSTFT(:,:,mm), ...
        params.analysisWin, params.synthesisWin, params.hop, ...
        params.Nfft, params.Fs);
end

for mm = 1:array.N*array.cdrMicN
    arrayEstimate(:,mm) = my_istft(arrayEstimateSTFT(:,:,mm), ...
        params.analysisWin, params.synthesisWin, params.hop, ...
        params.Nfft, params.Fs);
end
minLength = min([size(directEstimate,1),size(array.arraySignal{1}(:,1),1)]);

nmse = @(x) mean(sum(abs(x * arrayEstimate(1:minLength,:) - ...
    cat(2,array.meanDerev{:})).^2) ./ sum(abs(cat(2,array.meanDerev{:}).^2)));
bestScale = fminbnd(nmse, 0, 10);
% Scale the estimate for honest evaluation
directEstimate = bestScale * directEstimate;
cptPts.directEstimate = directEstimate;

%% Estimate full signal with the diffuse component

[completeEstimate, diffuseContribution] = estimatecompletesignal(cptPts,...
    array, params);

for mm = 1:cptPts.N
    [directEstimate(:,mm), completeEstimate(:,mm)] = ...
        alignsignals(directEstimate(:,mm), completeEstimate(:,mm), [], ...
        'truncate');
end
cptPts.directEstimate = directEstimate;
cptPts.completeEstimate = completeEstimate;
end
