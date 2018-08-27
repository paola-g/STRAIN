% data matrix n_points x n_features
% c clustering labels (npoints x 1) 
% meanmds vector containing clusters means


function [meanmds] = calc_means(data,c)
uniq = unique(c);
[~,d] = size(data);
nn = length(uniq);
meanmds = zeros(nn,d);

for i=1:nn
    idxc = c == uniq(i);
    cluster = data(idxc, :);
    meanmds(i,:) = mean(cluster);
        
end
meanmds = meanmds';
end
