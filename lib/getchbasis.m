function [ B ] = getchbasis( th_ax,order,angle )
%GETCHBASIS Summary of this function goes here
%   Detailed explanation goes here
if nargin==3
    B=zeros(length(th_ax),order+1);
    for n=0:order
        B(:,n+1)=cos(n*(th_ax-angle));
    end
else
    Acos=zeros(length(th_ax),order+1);
    Asin=zeros(length(th_ax),order+1);
    for n=0:order
       Acos(:,n+1)=cos(n*(th_ax));
       Asin(:,n+1)=sin(n*(th_ax));
    end 
    B=[Acos,Asin];
end
end

