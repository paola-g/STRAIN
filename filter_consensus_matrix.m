% M: consensus matrix (dissimlarity matrix as returned by consensus_cl)
% th: 1 - th -> proportion of times two objects are clustered together
% nneigh: -> proportion of neighbours exceeding threshold
% filter_idx: index of filtered objects
function [filter_idx] = filter_consensus_matrix(M,th,nneigh)
n = size(M,1);
filter_idx = [];
for i=1:n
    flag = 0;
    for j=1:n
        if(j>i)
            if M(i,j) <= th
                flag = flag +1;
            end
        end
    end
    if(flag>=nneigh)
        filter_idx = [filter_idx i];
    end
end