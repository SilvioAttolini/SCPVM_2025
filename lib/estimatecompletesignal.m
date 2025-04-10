function [completeEstimate, diffuseContribution] = estimatecompletesignal(cptPts, array, params)
%% ESTIMATECOMPLETESIGNAL
% This function estimate the diffuse components and computes the complete
% signal (direct + diffuse).

addpath(genpath('habets'));

fprintf('Estimate the full signal (direct + diffuse)...\n');
fLen = params.fLen;
tLen = params.tLen;
arrayCenter = cell2mat(array.center);
arrayCenter = arrayCenter(:,1:2);

distArrayCptPts = pdist2(cptPts.position,  arrayCenter);

[~, closestArray]= min(distArrayCptPts, [], 2);
diffuseContribution = zeros(size(cptPts.directEstimate));

mean3 = @(x) sqrt(mean(abs(x).^2,3));
totalDiffuseSTFT = cellfun(mean3, array.meanDiffuseSTFT, ...
    'UniformOutput', false);

% Compute the diffuse contribution for each control point
prevArrayIdx = 0; % for the first iteration
%mic_of_couple = 1;
for mm = 1:cptPts.N
    distFactor = (1 ./ distArrayCptPts(mm,:).');
    distFactor = distFactor ./ sum(distFactor);

    arrayIdx = closestArray(mm);
    % fprintf("closest array: ");
    % display(arrayIdx);

    diffuseSum = zeros(fLen, tLen);
    
    % the amplitude is the weighted sum
    for dd = 1:array.N
        diffuseSum = diffuseSum + (distFactor(dd) * totalDiffuseSTFT{dd});
    end

    whichMic = 1;
    if arrayIdx == prevArrayIdx
        whichMic = whichMic + 1;
    end

    diffuseSum = diffuseSum .* ...
        exp(1j*angle(array.meanDiffuseSTFT{arrayIdx}(:,:,whichMic)));
       
    diffuseContribution(:,mm) = my_istft(diffuseSum, params.analysisWin, ...
           params.synthesisWin, params.hop, params.Nfft, params.Fs);

    prevArrayIdx = arrayIdx;
end

% applay habets' sc here
num_of_couples = cptPts.N / 2;
vm = 1;
for cc = 1:num_of_couples
    couple_diffuse = [diffuseContribution(:,vm).', ...
                      diffuseContribution(:,vm+1).'].'; 

    couple_diff = habets(params, cptPts, couple_diffuse);

    [adj_diff_a, adj_diff_b] = adjust_restrict_diffuses(...
        diffuseContribution(:, vm),diffuseContribution(:, vm+1),...
        couple_diff(:, 1), couple_diff(:, 2));

    diffuseContribution(:, vm) = adj_diff_a;
    diffuseContribution(:, vm+1) = adj_diff_b;
    
    vm = vm+2;  % align to next couple
end

% The complete sound field is given by direct + diffuse
completeEstimate = cptPts.directEstimate + diffuseContribution;


end
