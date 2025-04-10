function  h = sphericalhankel(nu,K,z)
%SPHERICALBASIS Summary of this function goes here
%   Detailed explanation goes here

h = sqrt(pi./(2*z)) .* besselh(nu+0.5,K,z);

end

