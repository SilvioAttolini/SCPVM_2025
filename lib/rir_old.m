function [IR, RTF,res ] = rir(c, Fs, mic_pos,s, L, beta, NFFT, mtype, nOrder,dim,orientation,doPlot,varargin)

if length(beta)==1
    if beta(1)~=0
        V = L(1)*L(2)*L(3);
        S = 2*(L(1)*L(2))+2*(L(1)*L(3))+2*(L(2)*L(3));
        alfa = 24*V*log(10.0)/(c*S*beta(1));
        if (alfa > 1)
            error('Error: The reflection coefficients cannot be calculated using the current room parameters, i.e. room size and reverberation time Please specify the reflection coefficients or change the room parameters.');
        else
            beta= zeros(1,6);
            for i=1:6
                beta(i) = sqrt(1-alfa);
            end
        end
    else
        beta= zeros(1,6);
    end
end

s = s';
mic_pos = mic_pos';
r = zeros(3,1);
res =[];
Rp_plus_Rm =zeros(3,1);
Rm = zeros(3,1);
for mx = -nOrder:nOrder
    Rm(1) = 2*mx*L(1);
    
    for my = -nOrder:nOrder
        Rm(2) = 2*my*L(2);
        
        for mz = -nOrder:nOrder
            Rm(3) = 2*mz*L(3);
            
            for q = 0:1
                Rp_plus_Rm(1) = (1-2*q)*s(1) - r(1) + Rm(1);
                refl(1) = beta(1)^(abs(mx-q)) * beta(2)^(abs(mx));
                
                for jj = 0:1
                    Rp_plus_Rm(2) = (1-2*jj)*s(2) - r(2) + Rm(2);
                    refl(2) = beta(3)^(abs(my-jj)) * beta(4)^(abs(my));
                    
                    for k = 0:1
                        Rp_plus_Rm(3) = (1-2*k)*s(3) - r(3) + Rm(3);
                        refl(3) = beta(5)^(abs(mz-k)) * beta(6)^(abs(mz));
                        
                        %                         dist = sqrt(Rp_plus_Rm(1)^2 + Rp_plus_Rm(2)^2 + Rp_plus_Rm(3)^2);
                        
                        if (abs(2*mx-q)+abs(2*my-jj)+abs(2*mz-k) <= nOrder || nOrder == -1)
                            
                            res =[res,[Rp_plus_Rm;refl(1)*refl(2)*refl(3)]];
                            
                        end
                    end
                end
            end
        end
    end
end

f_ax =linspace(0,Fs,NFFT+1)';
f_ax = f_ax(1:NFFT/2+1);
k = 2*pi*f_ax/c;
RTF = zeros(size(mic_pos,2),length(k));
if (~isempty(varargin)) && (~isempty(varargin{1}))
    earlyDistance = varargin{1}*c;
    directDistance = pdist2(res(1:3,:).', mic_pos.');
    [minVal, directIdx] = min(directDistance);
    lateDistance = minVal + earlyDistance;
    lateIdx = directDistance > lateDistance;
    lateIdx(directIdx) =  1;
    res = res(:,lateIdx);
end
    
for mic_indx = 1:size(mic_pos,2)
    gain = zeros(1,size(res,2));
    for src_indx = 1:size(res,2)
        x = res(1,src_indx)-mic_pos(1,mic_indx);
        y = res(2,src_indx)-mic_pos(2,mic_indx);
        z = res(3,src_indx)-mic_pos(3,mic_indx);
        
        gain(src_indx) = res(4,src_indx)*sim_microphone(x, y, z, orientation, mtype)/(4*pi*pdist2(res(1:3,src_indx)',mic_pos(:,mic_indx)'));
    end
    
    dist = pdist2(res(1:3,:)',mic_pos(:,mic_indx)')';
    aux = exp(-1i*k*dist);
    aux = aux.*repmat(gain,length(k),1);
    RTF(mic_indx,:) = sum(aux,2);
end

IR = ifft([RTF,conj(RTF(:,end-1:-1:2))],[],2,'symmetric');

if doPlot==1
    figure
    plotcube(L,[0,0,0],.4,[0 1 0]); hold on
    scatter3(s(1),s(2),s(3),50,'filled'); 
    scatter3(mic_pos(1,:),mic_pos(2,:),mic_pos(3,:),50,'filled');
    scatter3(res(1,:),res(2,:),res(3,:));
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    grid on
    axis equal
    
end

end

function gain = sim_microphone(x, y, z, angle, mtype)

if (strcmp(mtype,'b') || strcmp(mtype,'c')|| strcmp(mtype,'s') || strcmp(mtype,'h'))
    
    switch(mtype)
        case 'b'
            rho = 0;
        case 'h'
            rho = 0.25;
        case 'c'
            rho = 0.5;
        case 's'
            rho = 0.75;
    end
    
    vartheta = acos(z/sqrt(x^2+y^2+z^2));
    varphi = atan2(y,x);
    
    gain = sin(pi/2-angle(2)) * sin(vartheta) * cos(angle(1)-varphi) + cos(pi/2-angle(2)) * cos(vartheta);
    gain = rho + (1-rho) * gain;
    
else
    gain = 1;
end
end

