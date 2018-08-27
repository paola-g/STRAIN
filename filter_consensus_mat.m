% M: consensus matrix
% sigma: real value between 0 and 1, how often two features have to be clustered together to be considered a stable pair.
% nu: a feature is selected if it is included in a number of stable pairs greater than nu
% filter_idx: index of filtered features
function [filter_idx] = filter_consensus_mat(M,sigma,nu)
n = size(M,1);
filter_idx = [];
for i=1:n
    flag = 0;
    for j=1:n
        if(j>i)
            if M(i,j) <= sigma
                flag = flag +1;
            end
        end
    end
    if(flag>=nu)
        filter_idx = [filter_idx i];
    end
end