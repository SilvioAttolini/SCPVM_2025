function [meanDerev, meanDereverbSTFT, meanDiffuse, meanDiffuseSTFT] = dereverbarraywithdoa(array, source, dereverbParams, params, macro)
%% DEREVERBARRAYWITHDOA
% This function remove the reverberate part from the microphones signals

if strcmp(dereverbParams.sourcePosition, 'median')
    sourcePos = source.medianPosition;
elseif strcmp(dereverbParams.sourcePosition, 'best')
    sourcePos = source.bestPosition;
else
    sourcePos = source.position;
end
arraySTFT = array.arraySTFT;

fprintf('Compute the dereverberation with DOA...\n')
cfg.estimator = @estimate_cdr_robust_unbiased;    % unbiased, "robust" estimator (CDRprop2)

if ~isfield(dereverbParams,'cdrMicN')
    micN = array.micN;
else
    micN = dereverbParams.cdrMicN;
end
dereverbSTFT = cell(array.N, source.N);
dereverbSignal = cell(array.N, source.N);

for iSrc = 1:source.N
    for aa = 1:array.N
        micPosition = cell2mat(array.position(aa));
        micPosition = micPosition(:,1:2);
        for mm = 1:micN
            % TEST
            if micN == 1
                nextMic =  mod(mm,micN+1) + 1;
            else
                nextMic =  mod(mm,micN) + 1;
            end
                        
            tstSignal = my_istft(arraySTFT{aa}(:,:,mm), params.analysisWin,...
                params.synthesisWin, params.hop, params.Nfft, params.Fs);
            tstSignalNextMic = my_istft(arraySTFT{aa}(:,:,nextMic), params.analysisWin,...
                params.synthesisWin, params.hop, params.Nfft, params.Fs);
            
            winLength = 512;
            analysisWin = hamming(winLength,'periodic');
            synthesisWin = hamming(winLength,'periodic');
            hop = winLength / 4;
            Nfft = 4*2^nextpow2(winLength);
            tstFrequency = linspace(0,params.Fs/2, Nfft/2+1)'; % frequency axis
            testSTF{aa}(:,:,mm) = my_stft(tstSignal, analysisWin, hop, Nfft, params.Fs);
            testSTF{aa}(:,:,nextMic) = my_stft(tstSignalNextMic, analysisWin, hop, Nfft, params.Fs);
            PSD(:,:,mm) = estimate_psd(testSTF{aa}(:,:,mm), params.lambda);
            PSD(:,:,nextMic) = estimate_psd(testSTF{aa}(:,:,nextMic), params.lambda);
            centerMic = (micPosition(mm,:) + micPosition(nextMic,:)) / 2;
            micDist = pdist2(micPosition(mm,:), micPosition(nextMic,:));
            % Local x axis of the pair of microphones
            xAxis = (micPosition(mm,:) - micPosition(nextMic,:));
            xAxis = xAxis / norm(xAxis); % Normalization

            % Vector in the direction of the source
            dd = (sourcePos(iSrc,:) - centerMic);
            dd = dd / norm(dd);         % Normalization
            
            DOA =  acos(xAxis * dd');
            TDOA = micDist*cos(DOA)/params.c;
            % define coherence models
            Css = exp(1i * 2 * pi * tstFrequency * (TDOA));  % target signal coherence; not required for estimate_cdr_nodoa
            Cnn = sinc(2 * tstFrequency * micDist/params.c); % diffuse noise coherence; not required for estimate_cdr_nodiffuse
            
            crossPSD = estimate_cpsd(testSTF{aa}(:,:,mm), ...
                testSTF{aa}(:,:,nextMic), params.lambda) ./ ...
                sqrt(PSD(:,:,mm).*PSD(:,:,nextMic));
            
            CDR = cfg.estimator(crossPSD, Cnn, Css);
            CDR = max(real(CDR),0);
            
            % weights = spectral_subtraction(CDR, params.alpha, params.beta, ...
            %     params.mu, params.floor);
            % version 6
            curr_mu = params.mu(aa);
            weights = spectral_subtraction(CDR, params.alpha, params.beta, ...
                curr_mu);
%             weights = max(weights, params.floor);
%             weights = min(weights, 1);
%             if mm == 1
%                 figure
%                 imagesc(db(weights))
%                 title(['Source ', num2str(iSrc)])
%             end
            directFilter{aa}(:,:,mm,iSrc) = (weights);
            inputSignal = cat(3, testSTF{aa}(:,:,mm), ...
                testSTF{aa}(:,:,nextMic));
            postFilter = sqrt(mean(abs(inputSignal).^2,3)) .* ...
                exp(1j*angle(testSTF{aa}(:,:,mm)));
            inputDereverb{aa}(:,:,mm) = postFilter;
            dereverbSTFTLocal = sqrt(weights) .* postFilter;
            
            dereverbSignal{aa,iSrc}(:,mm) = ...
                my_istft(dereverbSTFTLocal, analysisWin, ...
                synthesisWin, hop, Nfft, params.Fs);
            if macro.PRINT_WIENER == true
                figure()
                subplot(2,1,1)
                imagesc(db(CDR))
                set(gca,'YDir','normal')
                caxis([-15 15])
                colorbar
                title('Estimated CDR (=SNR) [dB]')
                xlabel('frame index')
                ylabel('subband index')
                subplot(2,1,2)
                imagesc(weights)
                set(gca,'YDir','normal')
                caxis([0 1])
                colorbar
                title(['Filter gain mic: ', num2str(mm)])
                xlabel('frame index')
                ylabel('subband index')
                
                maxC = max(db(abs(postFilter(:))));
                minC = min(db(abs(postFilter(:))));
                
                figure()
                subplot(2,1,1)
                imagesc(db(abs(postFilter)), [minC, maxC])
                set(gca,'YDir','normal')
                colorbar
                title('Input of the Wiener filter [dB]')
                
                subplot(2,1,2)
                imagesc(db(abs(dereverbSTFT{aa,iSrc}(:,:,mm))), ...
                    [minC, maxC])
                set(gca,'YDir','normal')
                colorbar
                title('Output of the Wiener filter [dB]')
            end
        end
    end
end

for aa = 1:array.N
   meanDirectFilter{aa} = mean(directFilter{aa}, 4);
   diffuseFilter{aa} = sqrt(1 - meanDirectFilter{aa}.^2); % Diffuse filter
   meanDereverbSTFTLocal = meanDirectFilter{aa} .* inputDereverb{aa};
   meanDiffuseSTFTLocal{aa} = diffuseFilter{aa} .*  inputDereverb{aa};
   for mm = 1:micN
       meanDerev{aa}(:,mm) = my_istft(meanDereverbSTFTLocal(:,:,mm), analysisWin, ...
           synthesisWin, hop, Nfft, params.Fs);
       meanDereverbSTFT{aa}(:,:,mm) = my_stft(meanDerev{aa}(:,mm), params.analysisWin, params.hop,...
           params.Nfft, params.Fs);
       meanDiffuse{aa}(:,mm) = my_istft(meanDiffuseSTFTLocal{aa}(:,:,mm), analysisWin, ...
           synthesisWin, hop, Nfft, params.Fs);
       meanDiffuseSTFT{aa}(:,:,mm) = my_stft(meanDiffuse{aa}(:,mm), params.analysisWin, params.hop,...
           params.Nfft, params.Fs);
   end
end


end