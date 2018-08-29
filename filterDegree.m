% deg punteggio per ogni cluster
% n numero di punti
% thresh soglia
% idxs indici logici punti dei cluster filtrati
function [idxs, cls] = filterDegree(c, deg, thresh)
cls = find(deg>thresh);
u = unique(c);
idxs = ismember(c, u(cls));
end