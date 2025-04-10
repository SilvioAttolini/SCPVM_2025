function hCoeff = sphericalharmonicsestimation(array, source, sphParams, params, macro)
%% SPHERICALHARMONICSESTIMATION
% This function estimates the coefficients of the spherical harmonics
% expansion

fprintf('Esimate direct signal spherical harmonics expansion...\n');
if ~isfield(sphParams,'cdrMicN')
    micPos = cell2mat(array.position);
    micPos = micPos(:,1:2);

else
    micPos = [];
    for aa = 1:array.N
        micL = array.position{aa}(1:sphParams.cdrMicN,1:2);
        micPos = [micPos; micL];
    end
end

% sphCoef = cell(source.N, 1);
frequency = params.frequency;
kvec = 2*pi*frequency/params.c;
L = params.tLen;
regParam = sphParams.regParam;
type = sphParams.type;
if strcmp(sphParams.arraySignal, 'direct')
    meanDereverbSTFT = arrayDirectSTFT;
else
    meanDereverbSTFT = array.meanDerevSTFT;
end

if strcmp(sphParams.sourcePosition, 'median')
    sourcePos = source.medianPosition;
elseif strcmp(dereverbParams.sourcePosition, 'best')
    sourcePos = source.bestPosition;
else
    sourcePos = source.position;
end
    

micSTFT = cat(3, meanDereverbSTFT{:});

maxOrder = sphParams.maxOrder;               % Spherical harmonics max order

if macro.modelType == 1   % 2D circular harmonics expansion
    % Coefficients of the harmonics expansion
    hCoeff = zeros(length(frequency),(2*maxOrder(1)+1)*source.N,L);
    for ff = 2:length(frequency)
        tempSignal = squeeze(micSTFT(ff,:,:)).';
        k = kvec(ff);
        basisF = circularbasis(micPos, sourcePos(:,1:2), k, maxOrder, ...
            type);
        basisFInv = pinv(basisF);
        hCoeff(ff,:,:) = (pinv(basisF) * tempSignal);
    end
else   % 3D spherical harmonics expansion
    % Coefficients of the harmonics expansion
    hCoeff = zeros(length(frequency),((maxOrder(1)+1)^2)*source.N,L);
    indexNotZero = correctindexes(maxOrder(1));
    if source.N == 2
        indexNotZero = [indexNotZero, ((maxOrder(1)+1)^2)+ indexNotZero];
    end
    for ff = 2:length(frequency)
        tempSignal = squeeze(micSTFT(ff,:,:)).';
        k = kvec(ff);
        basisF = sphericalbasis(micPos, sourcePos(:,1:2), k, ...
            maxOrder, type);
        
%         indexZero = mean(abs(basisF),1) < 1e-14; % zeros(1,size(basisF,2));
%         condNumber(ff) = cond(basisF(:,~indexZero));
        cleanBasisF = basisF(:,indexNotZero);
%         auxBasisF = pinv(basisF(:,~indexZero));
        auxBasisF = svd_inverse_matrix(cleanBasisF, regParam);
        aux = (auxBasisF * tempSignal);
        hCoeff(ff,indexNotZero,:) = aux;
    end
end



end