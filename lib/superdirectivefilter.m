function H = superdirectivefilter(micPosition, thetaAx, fAx, c)
%% PSEUDOSPECTRUM
% This function computes the filter for an array over 360
% degrees with superdirective beamforming algorithm
%
% Params:
%   micPosition: the relative microphone position
%   c: sound speed
% Returns:
%   H: the filter in frequency
fLen = length(fAx);
micsN = size(micPosition, 1);

dl = 10^(-20/10);                   % Diagonal loading
% thetaAx = -180:stepAngle:180-stepAngle;
% thetaAx = deg2rad(thetaAx);

% Filter definition
% eigBeamforming = zeros(length(thetaAx)*(fLen), micsN);
% idx = 1;
gSteeringAll = zeros(micsN, length(thetaAx), fLen);     % Steering Vector
k = [cos(thetaAx)', sin(thetaAx)'];
for mm = 1:micsN
    for aa = 1:length(thetaAx)  
        gSteeringAll(mm, aa, :) = exp(1i*2*pi*fAx/c .* (k(aa,:) * ...
            micPosition(mm,:)')).';
    end
end

dist = squareform(pdist(micPosition));                  % Microphone dist
Gamma2 = zeros(micsN*micsN, fLen);                      % Cross-correlation

for mm = 1:micsN*micsN
    Gamma2(mm,:) = sinc(2*pi*fAx * dist(mm) / (pi*c));
end
GammaRS = reshape(Gamma2, micsN, micsN, fLen);
GammaRS = GammaRS + dl * max(max(GammaRS)) .* repmat(eye(micsN), 1,1, fLen);


for ff = 1:fLen
   gSt = gSteeringAll(:,:,ff);                      % Steeering vector
   gamma = GammaRS(:,:,ff);                         % Cross-correlation
   DD = diag(gSt' * (gamma \ gSt));                 % Denominator
   h(:,:,ff) =  (gamma \ gSt) ./ DD';               % Filter!
end
H = conj(h);          

% DIRECT IMPLEMENTATION
% idx = 1;
% for ff = 1:fLen
%         Gamma = squareform(pdist(micPosition));
%         Gamma = sinc(2*pi*fAx(ff) * Gamma / (pi*c));
%         Gamma = Gamma + dl * max(max(Gamma))*eye((micsN));
% %         GammaInv = inv(Gamma);
%     for aa = 1:length(thetaAx)
%         k = [cos(thetaAx(aa)), sin(thetaAx(aa))];
%         gSteering = exp(1i*2*pi*fAx(ff) / c*(k*micPosition')).';
%         h2 = (Gamma \ gSteering) / (gSteering' * (Gamma\gSteering));
% %         %         h = (GammaInv * gSteering) / (gSteering' * GammaInv * ...
% %         %             gSteering);
%         eigBeamforming(idx,:) = h2';
%         idx = idx+1;
%     end
% end
% 
% aux2 = repelem(signal, 1,1,length(thetaAx));        % Repeat signal for 360
% aux2 = permute(aux2,[2,3,1]);                       % Move indexes
% P2 = squeeze(sum(aux2 .* H, 1));      % OUTPUT
% 
% % DIRECT IMPLEMENTATION
% % aux = repelem(signal, length(thetaAx), 1);
% % P = reshape(sum(aux.*eigBeamforming, 2), ...
% %     length(thetaAx), (fLen));
% P = abs(P2);
end