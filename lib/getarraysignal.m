function [arraySignal, arraySTFT] = getarraysignal(array, source, room, params)
%% GETARRAYSIGNAL
% This function computes the signals propagating the source signals to the
% array microphones in the specified room.

arraySignal = cell(array.micN, 1);
arraySTFT = cell(array.micN, 1);
arrayDirectSTFT = arraySTFT;
arrayDirectSignal = arrayDirectSTFT;
% arraySDR = arrayDirectSTFT;
% arraySDRFrequency = arraySDR;
diffuseTime =  room.diffuseTime;
fLen = params.fLen;
tLen = params.tLen;

sourcePos = cell2mat(source.position');

for iSrc = 1:source.N
    for aa = 1:array.N
        if iSrc == 1
            arraySTFT{aa} = zeros(fLen, tLen, array.micN);
            arrayDirectSTFT{aa} = zeros(fLen, tLen, array.micN);
        end
        arrayPos = cell2mat(array.position(aa));         % Mics coordinates
        for mm = 1:array.micN
            [~, impResp] = rir(params.c, params.Fs, ...
                sourcePos(iSrc,:), arrayPos(mm,:), room.dim, ...
                room.T60, params.Nfft, source.type{iSrc}, ...
                room.reflectionOrder, 3, source.orientation{iSrc}, false, ...
                diffuseTime);
            impRespFrame = repmat(impResp.', 1, tLen);
            current = source.sourceSTFT{iSrc} .* impRespFrame;
            arraySTFT{aa}(:,:,mm) = arraySTFT{aa}(:,:,mm) + current;
            
            [~, h] = rir(params.c, params.Fs, sourcePos(iSrc,:),...
                arrayPos(mm,:), room.dim, 0, params.Nfft,...
                source.type{iSrc}, 0, 3, source.orientation{iSrc}, false);
            hFrame = repmat(h.', 1,tLen);
            current = source.sourceSTFT{iSrc} .* hFrame;
            arrayDirectSTFT{aa}(:,:,mm) = arrayDirectSTFT{aa}(:,:,mm) + current;
        end
    end
end

dix = 1;
for aa = 1:array.N
    for mm = 1:array.micN
        if mod(dix, 39) == 0, fprintf('\n'); else, fprintf('.'); end
        arraySignal{aa}(:,mm) = my_istft(arraySTFT{aa}(:,:,mm), params.analysisWin, ...
            params.synthesisWin, params.hop, params.Nfft, params.Fs);
        arrayDirectSignal{aa}(:,mm) = my_istft(arrayDirectSTFT{aa}(:,:,mm),...
            params.analysisWin, params.synthesisWin, params.hop, ...
                params.Nfft, params.Fs);
        varS = var(arraySignal{aa}(:,mm));          % Signal energy
        varN = varS / (10^(params.SNR/10));         % Noise energy
        % Generate microphones noise
        noise = sqrt(varN) * randn(1,  size(arraySignal{aa},1));
        % Noisy microphone signal
        arraySignal{aa}(:,mm) = arraySignal{aa}(:,mm) + noise';
        % STFT the microphone signal
       noiseSTFT = my_stft(noise, params.analysisWin, params.hop, ...
           params.Nfft, params.Fs);
       arraySTFT{aa}(:,:,mm) = arraySTFT{aa}(:,:,mm) + noiseSTFT;
       dix = dix + 1;
%        arraySDR{aa}(mm) = bss_eval_sources(arraySignal{aa}(:,mm).', arrayDirectSignal{aa}(:,mm).');
%        [arraySDRFrequency{aa}(:,mm), F0] = sdrfrequency(arraySignal{aa}(:,mm), arrayDirectSignal{aa}(:,mm), '1/3 octave', params.Fs);
    end
end
fprintf('\n')

%% % Theoretical DRR of microphone arrays for the best weight assignment
mean3 = @(x) sqrt(mean(abs(x).^2,3));
totalDirectSTFT = cellfun(mean3, arrayDirectSTFT, ...
    'UniformOutput', false); % mean for each array
array_diffuse = arraySignal;
array_diffuse_stft = arraySTFT;
for aa = 1:array.N
    for mm = 1:array.micN
        array_diffuse{aa}(:,mm) = arraySignal{aa}(:,mm) - arrayDirectSignal{aa}(:,mm);
        array_diffuse_stft{aa}(:,:,mm) = my_stft(array_diffuse{aa}(:,mm), params.analysisWin, params.hop, params.Nfft, params.Fs);
    end
end
totalDiffuseSTFT = cellfun(mean3, array_diffuse_stft, ...
    'UniformOutput', false);

DRRs_array = zeros(1, array.N);
for aa = 1:array.N
time_direct_array = my_istft(totalDirectSTFT{aa}, params.analysisWin, ...
           params.synthesisWin, params.hop, params.Nfft, params.Fs);
time_diffuse_array = my_istft(totalDiffuseSTFT{aa}, params.analysisWin, ...
           params.synthesisWin, params.hop, params.Nfft, params.Fs);
% disp(size(time_direct_array));
% disp(size(time_diffuse_array));

power_direct_array = mean(abs(fft(time_direct_array)).^2, 2);
power_diffuse_array = mean(abs(fft(time_diffuse_array)).^2, 2);
DRRs_array(1, aa) = power_direct_array./power_diffuse_array;
end

% DRR of arrays
fig = figure('Visible', 'off');
plot(db(DRRs_array))
xlabel('VM index');
ylabel('[dB]');
title('DRR of arrays')
grid on;
saveas(fig, 'results/DRR_of_arrays.png');
close(fig);

end
