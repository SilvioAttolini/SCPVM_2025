function [medianPosition, bestPosition, medianError, bestError] = ...
    sourcelocalization(array, source, localizationParams, params, macro)
%% SOURCELOCALIZATION
% This function localize the N acoustic sources using the Projective Ray
% Space and RANSAC algorithm.

fprintf('Localize the sources...\n')
arraySTFT = array.arraySTFT;
frequency = params.frequency;
tAx = params.tAx;
fLen = params.fLen;

if localizationParams.dereverb == true
    [dereverbSTFT, ~] = dereverbarraynodoa(array, params, macro);
else
    dereverbSTFT = arraySTFT;
end

secPerPosEstimation = localizationParams.secPerPosEstimation;

[nominalMicPos(:,1),nominalMicPos(:,2)] = pol2cart(0:2*pi/array.micN ...
    :2*pi*(1-1/array.micN), array.radius);

stepAngle = localizationParams.stepAngle;
thetaAx = 0:stepAngle:360-stepAngle;
thetaAx = deg2rad(thetaAx);


% Normalize the center corrdinates translating to the center
arrayCenter = cell2mat(array.center);
arrayCenter = arrayCenter(:,1:2);
[radiusROI, centerROI] = ApproxMinBoundSphereND(arrayCenter(:,1:2));
arrayCenter = arrayCenter - centerROI;
arrayCenter_ = cell(array.N, 1);
pS = cell(array.N,1);
arrayCenterStruct = cell(1);
srcPos = cell2mat(source.position');
srcPos = srcPos(:,1:2);


fprintf('Localization based on beamforming...\n')

sdFilter = superdirectivefilter(nominalMicPos, thetaAx, frequency, ...
    params.c);

if secPerPosEstimation == -1
    timeFrameStep = 1;
else
    [~,timeFrameStep] = min(abs(tAx-secPerPosEstimation));
end

timeAxis = 1:timeFrameStep:length(tAx)-timeFrameStep;
allDoa = zeros(source.N, array.N, length(timeAxis));
allWeight = allDoa;
doaHist = cell(size(timeAxis));
doaVal = cell(size(timeAxis));
for aa = 1:array.N            % Averaging signal on time
    fprintf(['\nArray ', num2str(aa), '...\n'])
    dix = 1;
    for idx = 1:length(timeAxis)
        tt = timeAxis(idx);
        if mod(dix,10) == 0, fprintf('.'); end               % Feedback
        if mod(dix,39) == 0, fprintf('\n'); end
        dix = dix+1;
        
        stopFrame = min(timeFrameStep, size(dereverbSTFT{aa},2));
        micSignal = dereverbSTFT{aa}(:, tt+stopFrame, :);
        mask = [];
        micSignal = squeeze(mean(micSignal, 2));
        [doaWeight, doa, pS{aa}] = doaestimator(micSignal, sdFilter, ...
            thetaAx, mask);

        % Prune external DOA
        idxPrev = mod((aa-1)-1, array.N) + 1;
        idxNext = mod((aa-1)+1, array.N) + 1;
        tmpPrev = array.center{idxPrev} - array.center{aa};
        tmpSucc = array.center{idxNext} - array.center{aa};
        angPrev = atan2(tmpPrev(2), tmpPrev(1));
        angSucc = atan2(tmpSucc(2), tmpSucc(1));
        [~, idxDoa] = removeExternalDOA(doa, angSucc, angPrev);
        doaNew = doa(idxDoa);
        doaNew = (doaNew(1:min(length(doaNew), source.N)));
        doaN = length(doaNew);
        
        [~, doaIdx] = intersect(wrapTo2Pi(doa), doaNew, 'stable'); % indices for doa values
        doaVal = doaWeight(doaIdx);
        
        if isempty(doaNew)             % NO DOA FOUND
            fprintf('\bs')
            doaNew = NaN * ones(1,source.N);
            doaVal = doaNew;
            doaN = source.N;
            %                 continue;               % Skip!
        end
        
        allDoa(:,aa,idx) = doaNew';
        allWeight(:,aa,idx) = doaVal;
    end
end
fprintf('\n')                    % Feedback

arrayCenter = cell2mat(array.center);
arrayCenter = arrayCenter(:,1:2);

% dist = sourcePos(:,1:2) - arrayCenter;
% allDoa = atan2(dist(:,2), dist(:,1));
doaRaySpace = [];
weightRaySpace = [];
for aa = 1:array.N
    arrayAllDoa = squeeze(allDoa(:,aa,:));
    arrayAllWeight = squeeze(allWeight(:,aa,:));
    % Mapping doas Projective Ray Space
    arrayAllDoa = rmmissing(arrayAllDoa(:));
    arrayAllWeight = rmmissing(arrayAllWeight(:));
    alpha = 1;
    l1s = alpha' .* sin(arrayAllDoa) ;
    l2s = -alpha' .* cos(arrayAllDoa);
    l3s = alpha' .* ((arrayCenter(aa,2) * cos(arrayAllDoa)) ...
        - (arrayCenter(aa,1) * sin(arrayAllDoa)));
    doaRaySpace = [doaRaySpace; [l1s l2s l3s]];
    weightRaySpace=[weightRaySpace; arrayAllWeight];
end

if localizationParams.plotPRS == true
    figure
    scatter3(doaRaySpace(:,1), doaRaySpace(:,2), doaRaySpace(:,3));
    hold on
    [l1, l2] = ndgrid(-1:.1:1, -1:.1:1);
    sourcePos = cell2mat(source.position.');
    for ss = 1:source.N
        d = -[0,0,0]*sourcePos(ss,:)';
        l3 = (-sourcePos(ss,1)*l1 - sourcePos(ss,2)*l2);
        surf(l1,l2,l3);
    end
    title('Point in the Projective Ray Space');
end

localizationTest = localizationParams.localizationTest;
sourcePos = cell2mat(source.position.');
sourcePos = sourcePos(:,1:2);
bestError = inf*ones(source.N,1);
bestPosition = zeros(size(sourcePos));
allEstimatedPosition = zeros(source.N, 2, localizationTest);

fprintf(['Estimate source position with RANSAC in ', ...
    num2str(localizationTest), ' tests.\n'])
for tst = 1:localizationTest
    % Localize for each time frame
    estimatedPosition = (cluster_and_estimate(doaRaySpace.', weightRaySpace.', ...
        source.N, 10^-2))';
    if source.N > 1 % Two sources are present in the scene
        v = zeros(source.N,1);
        srcIdx = zeros(source.N,1);
        for ss = 1:source.N
            % Find the closest estimated position
            distance = (pdist2(sourcePos, estimatedPosition(ss,:)));
            [v(ss), idxM] = min(distance);
            if idxM ~= ss, srcIdx(ss) = idxM; else, srcIdx(ss) = ss; end
        end
        if srcIdx(1) == srcIdx(2)   % Swap indexes if necessary
            [~,idMin] = min(v);
            if idMin == 1
                srcIdx(2) = mod(srcIdx(2),2) + 1;
            else
                srcIdx(1) = mod(srcIdx(2),2) + 1;
            end
        end
        estimatedPosition = estimatedPosition(srcIdx,:);    % Reorder!!
    end
    localizationError = (localizationerror(estimatedPosition, sourcePos));
    if localizationError < bestError
        bestError = localizationError;
        bestPosition = estimatedPosition;
    end
    allEstimatedPosition(:,:,tst) = estimatedPosition;
end

medianPosition = median(allEstimatedPosition, 3);
medianError = localizationerror(medianPosition, sourcePos);
end