function H = sphericalbasis(micPosition, ...
    sourcePosition, waveNumber, maxOrder, type)
%SPHERICALBASIS Summary of this function goes here
%   Detailed explanation goes here

nMic = size(micPosition, 1);
nSource = size(sourcePosition, 1);
H = cell(1,nSource);

for ii = 1:nSource
    radius = pdist2(micPosition, sourcePosition(ii, :));
%     disp(radius)
    dd = sourcePosition(ii,:) - micPosition;
    phi = (atan2(dd(:,2), dd(:,1)));
    phi = phi(:);
%     n = -maxOrder(ii):maxOrder(ii);
%     expo = exp(1i*n'*phi(:)').';
    
%     h = zeros(length(radius), length(n));
    H{ii} = zeros(nMic, (maxOrder(ii)+1)^2, length(waveNumber));
    id = 1;
    for n = 0:maxOrder(ii)
        h = sphericalhankel(n, type, radius(:)*waveNumber(:).');
        
        for mm = -n:n
            sphHarmonic = harmonicY(n, mm, pi/2*ones(size(phi)), phi);
            H{ii}(:,id,:) = repmat(sphHarmonic, 1,length(waveNumber)) .* h;
            id = id+1;
        end
    end
    
%     H{ii} = h .* expo;
end

H = cell2mat(H);

end

