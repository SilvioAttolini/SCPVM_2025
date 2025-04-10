function [outSignalCpts,outSignalCptsSTFT] = propagate(hCoef,cptPts,sourcePos,maxOrder,type,kvec,params,modelType,varargin)

L = size(hCoef,3);
outSignalCptsSTFT = zeros(length(kvec),params.tLen,size(cptPts,1));
for ff = 2:length(kvec)
    if maxOrder == 0
        tempHcoeffs = squeeze(hCoef(ff,:,:));
    else
        tempHcoeffs = squeeze(hCoef(ff,:,:));
%         tempHcoeffs = tempHcoeffs(:);
    end
    k = kvec(ff);
    if modelType == 1
        basisF = circularbasis(cptPts, sourcePos, k, maxOrder, type);
    elseif modelType == 2
        basisF = sphericalbasis(cptPts, sourcePos, k, maxOrder, type);
    end
    outSignalCptsSTFT(ff,:,:) = (basisF*tempHcoeffs).';
end

if ~isempty(varargin)
orientation = varargin{1};
cptPts_type = varargin{2};

angAttenuation = computeAttenuation(cptPts,sourcePos,orientation,cptPts_type);
outSignalCptsSTFT = outSignalCptsSTFT.*...
    repmat(reshape(angAttenuation,1,1,length(angAttenuation))...
    ,size(outSignalCptsSTFT,1), size(outSignalCptsSTFT,2));
end
    
xlen = params.winLength + (params.tLen-1)*params.hop;
outSignalCpts = zeros(size(cptPts,1),xlen);
for ii = 1:size(cptPts,1)
    outSignalCpts(ii,:) = my_istft(squeeze(outSignalCptsSTFT(:,:,ii)), params.analysisWin,...
        params.synthesisWin, params.hop, params.Nfft, params.Fs);
end

end

function angAttenuation = computeAttenuation(cptPts,sourcePos,orientation,cptPts_type)
    
Ncpts = size(cptPts,1);
vec1 = repmat(sourcePos,Ncpts,1)-cptPts;
vec1 = normr(vec1);

vec2 = [cos(orientation(:)),sin(orientation(:))];
vec2 = normr(vec2);

ang = acos(sum(vec1.*vec2,2));

rho = ones(Ncpts,1);
rho(cptPts_type == 'b') = 0;
rho(cptPts_type == 'h') = 0.25;
rho(cptPts_type == 'c') = 0.5;
rho(cptPts_type == 's') = 0.75;

angAttenuation = rho +(1-rho).*cos(ang);

end

