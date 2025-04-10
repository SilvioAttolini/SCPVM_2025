function [ new_doa,idx ] = removeExternalDOA(doa,low_ang,up_ang)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

doa=wrapTo2Pi(doa);
low_ang=wrapTo2Pi(low_ang);
up_ang=wrapTo2Pi(up_ang);
if low_ang<=up_ang
   new_doa=doa(doa>= low_ang & doa<= up_ang);
   idx=doa>= low_ang & doa<= up_ang;
else
    doa(doa<=up_ang)=doa(doa<=up_ang)+2*pi;
    new_doa=doa(doa>= low_ang & doa<= up_ang+2*pi);
    idx=doa>= low_ang & doa<= up_ang+2*pi;
end
end

