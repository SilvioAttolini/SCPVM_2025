function LE = localizationerror(estimatedPosition, sourcePosition)
% function LOCALIZATIONERROR
% This function computes the localization error between two sources.
% Parameters:
%   - estimatedPosition: [x, y] coordinates of the estimated location
%   - sourcePosition: [x, y] coordinates of the actual location
% Returns:
%   - LE: localization error as the distances of the estimated sources with
%   respect to the actual source position.

    dist = pdist2(estimatedPosition, sourcePosition);
    LE = min(dist, [], 1);
end