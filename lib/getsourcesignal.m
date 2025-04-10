function [sourceSignal, sourceSTFT, source] = getsourcesignal(source, params)
%% getsourcesignal
% This function computes the source signals

fprintf('Source signal computation...\n');

if source.N == 2
    source.signalType{1} = 'speech';
    source.signalType{2} = 'speech';
else
    source.signalType{1} = 'speech';
end


sourceSignal = cell(source.N, 1);
sourceSTFT = sourceSignal;
nameNum = 49;
for iSrc = 1:source.N
    if strcmp(source.signalType{iSrc}, 'white')         % white noise
        tmp = randn(1,params.Fs*source.signalLength);
        tmp = [zeros(1, 2*params.winLength) tmp];
        sourceSignal{iSrc} = tmp;
    elseif strcmp(source.signalType{iSrc}, 'chirp')     % sinusoidal sweep
        tmp = chirp(0:1/params.Fs:(source.signalLength-1/params.Fs), 600, ...
            source.signalLength, 3500);
        tmp = [zeros(1, 2*params.winLength) tmp];
        sourceSignal{iSrc} = tmp;
        if source.N == 2
            if iSrc == 1
                tmp = chirp(0:1/params.Fs:(2-1/params.Fs), 10, 2, 3500);
                sourceSignal{iSrc} = tmp;
            else
                tmp = chirp(0:1/params.Fs:(2-1/params.Fs), 4000, 2, 7990);
                sourceSignal{iSrc} = tmp;
            end
        end
    else                                                % speech signal
        [tmp, oFs] = audioread([source.filePath, '/' num2str(nameNum), '.flac']);
        start = find(tmp > 0.5 * var(tmp), 1);          % signal is present
        stop = start + oFs * source.signalLength;
        tmp = tmp(start:stop-1, 1);
        tmp = resample(tmp, params.Fs, oFs);
        tmp = [zeros(2*params.winLength,1); tmp];
        sourceSignal{iSrc} = tukeywin(length(tmp), 0.99)' .* tmp';
        nameNum = nameNum +1;
    end
    sourceSignal{iSrc} = normalize(sourceSignal{iSrc}');
    sourceSignal{iSrc} = [zeros(size(params.winLength,1)); sourceSignal{iSrc}; zeros(size(params.winLength,1))]; 
%     source.signalLength = length(sourceSignal{iSrc});
    [sourceSTFT{iSrc}, ~, source.tAx] = my_stft(sourceSignal{iSrc}, ...
        params.analysisWin, params.hop, params.Nfft, params.Fs);
end
source.signalLength = length(sourceSignal{1}) / params.Fs;


end