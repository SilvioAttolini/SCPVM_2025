%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPLY LETTER SCRIPT
% Evaluate the performance including the early reflections in the signals.%
%
% v. 0.1
% Mirco Pezzoli
% 05/06/2020
%
% v. 0.2
% Silvio Attolini
% 23/02/2025
% adds Spatial Coherence consraints to VM diffuse component synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

fprintf('Test script for function parametricvirtualmiking.\n');

%% Setup
fprintf('Setting up...\n');
addpath(genpath('lib'));
addpath(genpath('harmonicY'));
addpath(genpath('audio'));

macro.modelType = 2; % 1: 2D, 2: 3D
macro.LOCALIZATION_PRS = false; %true
macro.PRINT_WIENER = false;
macro.PRINT_SETUP = false;

% STFT parameters
params.c = 342;
params.Fs = 16000;
params.winLength = 4096;
params.analysisWin = hamming(params.winLength,'periodic');
params.synthesisWin = hamming(params.winLength,'periodic');
params.hop = params.winLength / 8;
params.Nfft = 2*2^nextpow2(params.winLength);
params.SNR = 60;
params.SDR = inf;

params.lambda = 0.68; % smoothing factor for PSD estimation
params.alpha = 1;
params.beta = 1; % magnitude subtraction
%params.mu = 1.6; % Weight of wiener filter. was 1.3, ver 5

% wiener filter weights for current scenario
h = 1.8; % high_weight -> regions where direct >> diffuse
m = 2.3; % medium-weight -> intermediate regions
l = 5; %low-weight -> regions where direct << diffuse -> behind or far
params.mu = [h h m h m l l m h];
weight_h = h;
weight_m = m;
weight_l = l;
params.floor = 10^(-30/20);

frequency = linspace(0,params.Fs/2,params.Nfft/2+1)'; % frequency axis
fLen = length(frequency);
params.fLen = fLen;
params.frequency = frequency;
% Source settings
source.N = 2; % 2;
source.signalType = cell(source.N, 1);
source.signalLength = 5;                    % Signal length in [s]
source.filePath = './audio';
% Room settings
room.x = 5;                             % x length
room.y = 4;                             % y length
room.z = 3;                             % z length
room.dim = [room.x, room.y, room.z];    % room dimensions
room.volume = room.x * room.y * room.z;
room.surface = room.x*room.z*2 + room.y*room.z*2 + room.x*room.y*2;
room.T60 = 0.4;                         % room T60
room.reflectionOrder = 20;
room.diffuseTime = [];                  % Consider also the early
% reflections

%% Load the source signals and define their location
[sourceSignal, sourceSTFT, source] = getsourcesignal(source, params);
tLen = length(source.tAx);
params.tLen = tLen;
tAx = source.tAx;
params.tAx = tAx;

if source.N == 2
    source.position{1} = [1.75, 2, (room.z/2)];          % source positions
    source.position{2} = [3.25, 2.75 (room.z/2)];
    source.orientation{1} = [pi/4 0];             % Source looking angle
    source.orientation{2} = [-pi/2 0];

    source.type{1} = 'c';                   % Source directivity
    source.type{2} = 'c';

    if source.type{1} == 'c'                % Source pattern coefficients
        source.coefficient{1} = 0.5;
    else
        source.coefficient{1} = 0;
    end

    if source.type{2} == 'c'
        source.coefficient{2} = 0.5;
    else
        source.coefficient{2} = 0;
    end
else
    source.position{1} = [1.75, 2, (room.z/2)];          % source positions
    source.orientation{1} = [pi/4 0];              % Source looking angle
    source.type{1} = 'c';

    if source.type{1} == 'c'
        source.coefficient{1} = 0.5;
    else
        source.coefficient{1} = 0;
    end

end

%% Place the array in the scene
array.N = 9;                % Number of arrays
array.micN = 4;             % Number of microphone per array
array.radius = 0.04;       % Radius of the array

array = placearray(array, room);
fprintf(['Placed ', num2str(array.N), ' arrays in the room...\n']);


%% Place the virtual mics
th_ax = linspace(0,2*pi,37); %th_ax = linspace(0,2*pi,5);
th_ax = th_ax(1:(37-1));
indices_to_keep = [1, 2, 8, 9, 15, 16, 22, 23, 29, 30];
filtered_th_ax = th_ax(indices_to_keep);
th_ax = filtered_th_ax(1:end) + deg2rad(1.5);
distCpts = 1;
[xx,yy] = pol2cart(th_ax,distCpts);
cptPts.position = [];
for ss = 1:source.N
    cptPts.position = [cptPts.position;[xx(:)+source.position{ss}(1),yy(:)+source.position{ss}(2)]];
end
cptPts.N = size(cptPts.position,1);
d = norm(cptPts.position(1, :)-cptPts.position(2, :));
cptPts.distance = d;
fprintf(['Placed ', num2str(cptPts.N), ' VMs in the room...\n']);

%% Plot the setup

if macro.PRINT_SETUP == true
    fprintf('Plotting setup...\n');
    %figure(1)

    fig = figure('Visible', 'off');
    % Draw room
    rectangle('Position', [0,0,room.dim(1:2)], 'LineWidth', 1);
    hold on
    h = zeros(3, 1);
    
    for aa = 1:array.N
        % Draw mics
        tmp = cell2mat(array.position(aa));
        h(1) = scatter(tmp(:,1),tmp(:,2),30,'bo','filled');
    end
    % Draw source
    sourcePos = cell2mat(source.position');
    h(2) = scatter(sourcePos(:,1), sourcePos(:,2), 50, 'red', ...
        'diamond', 'filled');
    for iSrc = 1:source.N
        [tt, rr] = ch2pol(180, source.coefficient{iSrc}', ...
            source.orientation{iSrc}(1));
        [xPatt, yPatt] = pol2cart(tt, abs(rr*0.5));
        plot(xPatt + sourcePos(iSrc,1), yPatt + sourcePos(iSrc,2), ...
            'red', 'linewidth', 1.5)
    end
    % Draw virtual mics
    scatter(cptPts.position(:,1), cptPts.position(:,2), 50, 'k', 'square', 'filled');
    %grid on;
    %axis equal
    %     legend(h,'Mics','Sources','Location', 'best');
    %xlabel('x [m]')
    %ylabel('y [m]')
    %title('Geometric setup')
    grid on;
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    title('Geometric setup')
    saveas(fig, 'scene.png');
    fprintf("Setup saved as scene.png\n");

    close(fig);
end

%% Compute the array microphone signals
fprintf('Computing the microphone signals...\n');
source.sourceSTFT = sourceSTFT;

[arraySignal, arraySTFT] = getarraysignal(array, source, room, params);

array.arraySTFT = arraySTFT;
array.arraySignal = arraySignal;

%% Compute the reference signals at the contol points
sphParams.sourcePosition = "";
sphParams.arraySignal = 'estimate';
sphParams.maxOrder = 1;
sphParams.cdrMicN = array.micN;
sphParams.regParam.method = 'tikhonov';
sphParams.regParam.nCond = 35;
sphParams.type = 2;

[directReferenceSTFT, completeReferenceSTFT, directSourceReferenceSTFT] = ...
    getreferencesignal(cptPts, source, room, sphParams, params);
time_l = size(arraySignal{1},1);
directReference = zeros(time_l, cptPts.N);
completeReference = zeros(time_l, cptPts.N);
for mm = 1:cptPts.N
    directReference(:,mm) = my_istft(directReferenceSTFT(:,:,mm), params.analysisWin,...
        params.synthesisWin, params.hop, params.Nfft, params.Fs);
    completeReference(:,mm) = my_istft(completeReferenceSTFT(:,:,mm), params.analysisWin,...
        params.synthesisWin, params.hop, params.Nfft, params.Fs);
end

%% Export the GT audio files, stereo for each couple
Fs = params.Fs;
i = 1;
while i < cptPts.N
    audio_of_curr_couple = [completeReference(:, i), completeReference(:, i+1);]; 
    
    % Normalize the audio to prevent clipping
    max_val = max(abs(audio_of_curr_couple), [], 'all'); % Find the maximum absolute value
    if max_val > 1
        audio_of_curr_couple = audio_of_curr_couple / max_val; % Scale down the data
    end
    
    audiowrite(['audio_out/GT_couple_', num2str(i), '_', num2str(i+1), '.wav'],audio_of_curr_couple,Fs);
    i = i + 2;
end


%% Compute the VM signals!
[completeEstimate, directEstimate] = parametricvirtualmiking(array, ...
    source, cptPts, [], params, macro);

%% Export the synthesized audio files, stereo for each couple
Fs = params.Fs;
i = 1;
while i < cptPts.N
    audio_of_curr_couple = [completeEstimate(:, i), completeEstimate(:, i+1);]; 
    
    % Normalize the audio to prevent clipping
    max_val = max(abs(audio_of_curr_couple), [], 'all'); % Find the maximum absolute value
    if max_val > 1
        audio_of_curr_couple = audio_of_curr_couple / max_val; % Scale down the data
    end
    
    audiowrite(['audio_out/couple_', num2str(i), '_', num2str(i+1), '.wav'],audio_of_curr_couple,Fs);
    i = i + 2;
end


%% Spatial Coherence evaluation
fprintf('Spatial Coherence evaluation...\n')
K_eval = 256;

% evaluate the SC of the Reference vm pairs
diffuseReference = completeReference - directReference;
% Calculalte STFT and PSD of all ground truth signals
X_ref = stft(diffuseReference,'Window',hanning(K_eval),'OverlapLength',0.75*K_eval,'FFTLength',K_eval,'Centered',false);
X_ref = X_ref(1:K_eval/2+1,:,:);
phi_x_ref = mean(abs(X_ref).^2,2);

% evaluate the SC of the Estimated vm pairs signals
diffuseEstimate = completeEstimate - directEstimate;
% Calculalte STFT and PSD of all output signals
X_est = stft(diffuseEstimate,'Window',hanning(K_eval),'OverlapLength',0.75*K_eval,'FFTLength',K_eval,'Centered',false);
X_est = X_est(1:K_eval/2+1,:,:);
phi_x_est = mean(abs(X_est).^2,2);

% evaluate the Theoretical SC
d = norm(cptPts.position(1, :)-cptPts.position(2, :)); 
ww = 2*pi*Fs*(0:K_eval/2)/K_eval;
sc_theory = sinc(ww*d/(params.c*pi));

% Plot spatial coherence of each pair
% Theoretical vs Estimated vs Reference
% nmse_of_SCs = zeros(cptPts.N, 1);
K = 256; % = K_eval
sc_ref_bank = zeros(K/2+1, cptPts.N/2);
sc_est_bank = zeros(K/2+1, cptPts.N/2);  % num pairs = num mics / 2
vm = 1;
pair = 1;
while vm < cptPts.N
    % REFERENCE 
    % Compute cross-PSD of x_1 and x_(m+1)
    psi_x_ref =  mean(X_ref(:,:,vm) .* conj(X_ref(:,:,vm+1)),2);
    % Compute real-part of complex coherence between x_1 and x_(m+1)
    sc_ref = real(psi_x_ref ./ sqrt(phi_x_ref(:,1,vm) .* phi_x_ref(:,1,vm+1))).';
    sc_ref_bank(:,pair) = sc_ref;


    % ESTIMATED
    % Compute cross-PSD of x_1 and x_(m+1)
    psi_x_est =  mean(X_est(:,:,vm) .* conj(X_est(:,:,vm+1)),2);
    % Compute real-part of complex coherence between x_1 and x_(m+1)
    sc_est = real(psi_x_est ./ sqrt(phi_x_est(:,1,vm) .* phi_x_est(:,1,vm+1))).';
    sc_est_bank(:,pair) = sc_est;

    % Calculate NMSE Est vs THEORETICAL
    NMSE_est_vs_theor = db(new_nmse(sc_theory, sc_est));
    
    Freqs=0:(Fs/2)/(K/2):Fs/2;
    fig = figure('Visible', 'off');
    plot(Freqs/1000,sc_theory,'-k','LineWidth',1.5)
    hold on;
    plot(Freqs/1000,sc_est,'-.b','LineWidth',1.5)
    hold on;
    plot(Freqs/1000,sc_ref,'-.r','LineWidth',1.5)
    hold off;
    xlabel('Frequency [kHz]');
    ylabel('Real(Spatial Coherence)');
    title(sprintf('Inter sensor distance %1.2f m',d));
    legend('Theory',sprintf('Proposed Method (NMSE = %2.1f dB)',NMSE_est_vs_theor),...
        'Reference');
    grid on;
    saveas(fig, ['sc_out/sc_couple_', num2str(vm), '_', num2str(vm+1), '.png']);
    close(fig);

    vm = vm+2;  % skip the second mic of each pair, since
                % we consider two mics at a time
    pair = pair + 1;  % store the results for each pair
end

%% Compute the metrics on the full signal
%
%   METRICS
%
fprintf('Compute the Metrics...\n')


powerDirect = mean(directReference.^2, 1);
powerDirectEstimate = mean(directEstimate.^2, 1);

powerDiffuse = mean((completeReference - directReference).^2,1) ;
powerDiffuseEstimate = mean((completeEstimate - directEstimate).^2,1);

signalDiffuseRatio = powerDirect ./ powerDiffuse;
signalDiffuseRatioEstimate = powerDirectEstimate ./ powerDiffuseEstimate;


% DRR
fig = figure('Visible', 'off');
plot(db(signalDiffuseRatio))
hold on;
plot(db(signalDiffuseRatioEstimate))
hold off;
xlabel('VM index');
ylabel('[dB]');
legend('GT', 'Estimate')
title('DRR')
grid on;
saveas(fig, 'results/DRR.png');
close(fig);

% power of directs
fig = figure('Visible', 'off');
plot(db(powerDirect))
hold on;
plot(db(powerDirectEstimate))
hold off;
xlabel('VM index');
ylabel('[dB]');
legend('GT', 'Estimate')
title('Power of directs')
grid on;
saveas(fig, 'results/pow_directs.png');
close(fig);

% power of diffuses
fig = figure('Visible', 'off');
plot(db(powerDiffuse))
hold on;
plot(db(powerDiffuseEstimate))
hold off;
xlabel('VM index');
ylabel('[dB]');
legend('GT', 'Estimate')
title('Power of diffuses')
grid on;
saveas(fig, 'results/pow_diffuses.png');
close(fig);

% NMSE between Estimated and Reference Signals, for each vm
nmse_of_ffts = zeros(cptPts.N,1);
for vm = 1:cptPts.N
    nmse_of_ffts(vm,1) = new_nmse(abs(fft(completeReference(:, vm))), abs(fft(completeEstimate(:, vm))));
end

fig = figure('Visible', 'off');
plot(db(nmse_of_ffts))
xlabel('VM index');
ylabel('NMSE [dB]');
title('Power of Est vs GT')
grid on;
saveas(fig, 'results/NMSE_power_of_complete.png');
close(fig);

% SC nmse for each pair
% NMSE of Ref vs Theor SC
% NMSE of Est vs Ref SC
nmse_of_SCs_est_vs_theor = zeros(cptPts.N/2,1);
nmse_of_SCs_est_vs_ref = zeros(cptPts.N/2,1);
for pair = 1:cptPts.N/2
    nmse_of_SCs_est_vs_theor(pair) = new_nmse(sc_theory.', sc_est_bank(:,pair));
    nmse_of_SCs_est_vs_ref(pair) = new_nmse(sc_ref_bank(:,pair), sc_est_bank(:,pair));
end

fig = figure('Visible', 'off');
plot(db(nmse_of_SCs_est_vs_theor)) %,'-k','LineWidth',1.5)
hold on;
plot(db(nmse_of_SCs_est_vs_ref)) % ,'-.r','LineWidth',1.5)
hold off;
xlabel('Pair number')
ylabel('NMSE [dB]')
legend('vs Theor', 'vs Ref')
title('NMSE of Est vs Theor and Ref')
xlim([0, 11]);
grid on;
saveas(fig, 'results/NMSE_SC.png');
close(fig);

fprintf("Done.");
