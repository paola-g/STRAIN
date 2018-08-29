% cls: clustering matrix n_points x n_clusterings
% k: number of cluster in the final solution
% linkageMethod: as in matlab linkage function ('single', 'complete', 'average', ...)
%
function [cl,Z,M] = consensus_cl(cls,k,linkageMethod) 
[m,n] = size(cls);

M = zeros(m,m);
for i=1:n
    idx = cls(:,i);

    v = unique(idx);
    nclust = length(v);
    for l=1:nclust
        index = idx==v(l);
        M(index,index) = M(index,index) + 1;
    end
end
M = M./n;
M(logical(eye(size(M)))) = 1;
M = 1 - M;
M = squareform(M);
Z = linkage(M, linkageMethod);
cl = cluster(Z,'maxclust',k);
M = squareform(M);
end