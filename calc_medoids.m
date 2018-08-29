% corrmat: full correlation matrix (input of kmeans)
% c: clustering labels (npoints x 1)  
% n: number of clusters
% method: either 'mean' or 'sum'
% mds: vector containing medoids indices
% meancorr: vector containing mean correlation (or sum of correlations) of medoids with cluster points


function [mds, meancorr] = calc_medoids(corrmat,c,method)
uniq = unique(c);
nn = length(uniq);
meancorr = zeros(nn,1);
mds = zeros(nn,1);

for i=1:nn
    idxc = c == uniq(i);
    cluster = corrmat(idxc, idxc);
    switch(method)
        case 'sum'
            [meancorr(i),center] = max(sum(cluster));
        case 'mean'
            [meancorr(i),center] = max(mean(cluster));
        otherwise
    end
    tmp = find(idxc);
    mds(i) = tmp(center);
end
end
