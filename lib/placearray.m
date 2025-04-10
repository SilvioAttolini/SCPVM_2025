function array = placearray(array, room)
%%PLACEARRAY
% This function computes the array locations according to the number of
% arrays

array.center = cell(array.N, 1);
array.position = cell(array.N, 1);
[centerCoordinates(:,1), centerCoordinates(:,2)] = pol2cart(0:2*pi/array.N:2*pi*(1-1/array.N), 1.7);
centerCoordinates = centerCoordinates + room.dim(1:2)/2';

[micCoordinate(:,1),micCoordinate(:,2)] = ...
    pol2cart(0:2*pi/array.micN:2*pi*(1-1/array.micN), array.radius);

% Place the microphones in the environment.
for aa = 1:array.N
    array.center{aa} = [centerCoordinates(aa,:),(room.z/2)];         % Arrays center
    array.position{aa} = [micCoordinate,zeros(array.micN,1)]+ repmat(array.center{aa}, ...
        array.micN, 1);
end


end