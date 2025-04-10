function [est_pos,idx] = cluster_and_estimate(peaks_pos,peaks_value,num,th)
addpath(genpath('CodeRANSAC'));
dim = 2;
idx = ransac_subspaces(peaks_pos,dim,num,th);

est_pos = zeros(2,num);
for j = 1:num
    L = diag(peaks_value(idx == j))*(peaks_pos(:,idx == j)');%(peaks_pos(1:2,:)');
%     c = -peaks_pos(3,:)';
    [~,~,V] = svd(L'*L);
    est_pos(:,j) = V(1:2,end)/V(3,end);%pinv(L)*c;
end

end