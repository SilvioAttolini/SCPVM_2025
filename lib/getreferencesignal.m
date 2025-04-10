function [testDirectSTFT, testCompleteSTFT, testDirectSourceSTFT] = getreferencesignal(cptPts, source, room, sphParams, params)
%% GETREFERENCESIGNAL
% This function computes the direct and the complete signal at the control
% points (virtual microphones).

tLen = params.tLen;
fLen = params.fLen;
testDirectSTFT = zeros(fLen, tLen, cptPts.N);
testDirectSourceSTFT = cell(1,source.N);
testCompleteSTFT = testDirectSTFT;

sourceSTFT = source.sourceSTFT;
if strcmp(sphParams.sourcePosition, 'median')
    sourcePos = source.medianPosition;
elseif strcmp(sphParams.sourcePosition, 'best')
    sourcePos = source.bestPosition;
else
    sourcePos = cell2mat(source.position);
end
fprintf('Compute the VMs reference signals...\n');

for ss = 1:source.N
%     p_idx(:,ss) = pdist2(cptPts.position,sourcePos(ss,1:2))>0.25;
    
    for mm = 1:cptPts.N
        if mod(mm, 39) == 0, fprintf('\n'); else, fprintf('.'); end
        [~, h] = rir(params.c, params.Fs, source.position{ss},...
            [cptPts.position(mm,:),room.z/2], room.dim, 0, params.Nfft,...
            source.type{ss}, 0, 3, source.orientation{ss}, false);
        hFrame = repmat(h.', 1,tLen);
        
        % Compute the direct signal at the control points
        current = sourceSTFT{ss} .* hFrame;
        testDirectSTFT(:,:,mm) = testDirectSTFT(:,:,mm) + current;
        testDirectSourceSTFT{ss}(:,:,mm) = current;
    
        % Compute the complete sound field at the control points
        
        [~, h] = rir(params.c, params.Fs, source.position{ss},...
            [cptPts.position(mm,:),room.z/2], room.dim, room.T60, params.Nfft,...
            source.type{ss}, room.reflectionOrder, 3,...
            source.orientation{ss}, false, room.diffuseTime);
        hFrame = repmat(h.', 1,tLen);
        current = sourceSTFT{ss} .* hFrame;
        testCompleteSTFT(:,:,mm) = testCompleteSTFT(:,:,mm) + current;
    end
    fprintf('\n')
end

end