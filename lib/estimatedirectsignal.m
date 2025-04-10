function [directEstimateSTFT, directSourceEstimate, arrayEstimateSTFT] = estimatedirectsignal(cptPts, hCoef, array, source, sphParams, params, macro)
%% ESTIMATEDIRECTSIGNAL
% This function estimate the direct signal of some control points from
% using the spherical harmonics propagation
%
fprintf('Estimate the direct signal...\n');

fLen = params.fLen;
tLen = params.tLen;

if ~isfield(sphParams,'cdrMicN')
    micPos = cell2mat(array.position);
    micPos = micPos(:,1:2);
    arrayEstimateSTFT = zeros(fLen, tLen, array.N*array.micN);
else
    micPos = [];
    for aa = 1:array.N
        micL = array.position{aa}(1:sphParams.cdrMicN,1:2);
        micPos = [micPos; micL];
    end
    arrayEstimateSTFT = zeros(fLen, tLen, array.N*array.cdrMicN);
end

frequency = params.frequency;
kvec = 2*pi*frequency/params.c;
type = sphParams.type;
if strcmp(sphParams.sourcePosition, 'median')
    sourcePos = source.medianPosition;
elseif strcmp(dereverbParams.sourcePosition, 'best')
    sourcePos = source.bestPosition;
else
    sourcePos = source.position;
end

sourceSTFT = source.sourceSTFT;

directEstimateSTFT = zeros(fLen, tLen, cptPts.N);
% testDirectSTFT = directEstimateSTFT;
% testDirectSource = cell(1,source.N);
% testDirectSourceSTFT = testDirectSource;


% Figure for plots
% figure();
for ss = 1:source.N
    fprintf(['Estimate direct signal for source ', num2str(ss), '\n']);
    p_idx(:,ss) = pdist2(cptPts.position,sourcePos(ss,1:2))>0.25;
        
    N = sphParams.maxOrder(ss); % Spherical harmonics max order
    
    if macro.modelType == 1
        start = (2*N+1)*(ss-1)+1;
        stop = (2*N+1)*ss;
        if N == 0
            tmpCoeff = squeeze(hCoef(:,start:stop,:)).';
        else
            tmpCoeff = squeeze(hCoef(:,start:stop,:));
        end
    elseif macro.modelType == 2
        start = (N+1)^2*(ss-1)+1;
        stop = (N+1)^2*ss;
        if N == 0
            tmpCoeff = squeeze(hCoef(:,start:stop,:));
        else
            tmpCoeff = squeeze(hCoef(:,start:stop,:));
        end
    end
    
    [outSignalCpts,outSignalCptsSTFT] = propagate(tmpCoeff,cptPts.position,sourcePos(ss,1:2), ...
        N, type, kvec, params, macro.modelType);
    
    tt_ax = 0:1/params.Fs:1/params.Fs*size(outSignalCpts,2);
    directSourceEstimate{ss} = outSignalCpts.';
    directEstimateSTFT = directEstimateSTFT + outSignalCptsSTFT;
    %     plotWaveField(outSignalCpts,size(XX),p_idx,tt_ax(1:end-1))% Simulate signal at the microphones
    [~,currentEstimate] = propagate(tmpCoeff,micPos,sourcePos(ss,1:2), ...
        N, type, kvec, params, macro.modelType);
    arrayEstimateSTFT = arrayEstimateSTFT + currentEstimate;
end


end