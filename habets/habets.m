function couple_signals = habets(params, cptPts, couple_noise)

    Fs = params.Fs;
    c = params.c;
    K = 256;  % optimized parameter
    M = 2; % calculate the SC between one couple at a time
    d = cptPts.distance;
    % type_nf = 'spherical'; % Type of noise field always assumed to be 3D
    L = size(cptPts.directEstimate);
    L = L(1);

    data = couple_noise;
    data = data - mean(data);  % removes continuous contribution (?)

    babble = zeros(L,M);
    for m=1:M
        babble(:,m) = data((m-1)*L+1:m*L);
    end

    %% Generate matrix with desired spatial coherence
    ww = 2*pi*Fs*(0:K/2)/K;
    DC = zeros(M,M,K/2+1);
    for p = 1:M
        for q = 1:M
            if p == q
                DC(p,q,:) = ones(1,1,K/2+1);
            else
                DC(p,q,:) = sinc(ww*abs(p-q)*d/(c*pi));  % case 'spherical'
            end
        end
    end

    %% Generate sensor signals with desired spatial coherence
    couple_signals = mix_signals(babble,DC,'eigen');  % 'cholesky' / 'eigen'
end
