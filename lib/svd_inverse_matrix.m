function [G_INV, outParam] = svd_inverse_matrix(G,regParams)

method = regParams.method;
q      = regParams.nCond;

if strcmp(method,'truncated_svd')
    
    [U,E,V] = svd(G);
    N = size(G,1);
    M = size(G,2);
    
    
    sing_values = diag(E(1:M,1:M));
    max_sv = sing_values(1);
    k = length(find(sing_values >= max_sv/q));
    E2 = diag(1./sing_values(1:k));
    E2_inv = zeros(M,N);
    E2_inv(1:k,1:k) = E2;
    G_INV  = (V * E2_inv * U');
    numsv = k;
    outParam = numsv; % number of singular values
    
elseif strcmp(method,'tikhonov')
    if size(G,1) < size(G,2)
        G_INV = pinv(G);
        outParam = 0;
    else
    
    [~,E,~] = svd(G);
    M = size(G,2);
    
    sing_values = diag(E(1:M,1:M));
    max_sv = sing_values(1);
    min_sv=sing_values(end);
    tik=(max_sv^2-min_sv^2*q^2)/(q^2-1);

    if max_sv/min_sv<=q
        G_INV  = (G'*G)\G';  
        outParam = 0;
    else    
        G_INV  = (G'*G+tik*eye(M))\G';
        outParam = tik; % regularization parameter
    end
    end
end
