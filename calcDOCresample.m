% clusters: cluster labels
% clusTH: minimum number of objects in cluster
% idxFiltered: index of objects to be selected prior to computing DOC,
%	e.g. as returned by filter_consensus_matrix, or 1:no_objects if no filter is applied
function [deg] = calcDOCresample(clusters, clusTH, idxFiltered)
nCenters = length(unique(clusters));

if(~exist('idxFiltered','var'))
    idxFiltered = 1:length(clusters);
end

disp('Loading BVQX objects.');
mask = xff('wholebrain.msk');
mask4mm = resampleVolume(mask.Mask, 0.5, 3);
maskIdx = find(mask4mm>0.9);
[x y z]=ind2sub(size(mask4mm), maskIdx);
zmap=zeros(size(x));
disp('Calculating degree of clustering.');
zmap(idxFiltered)=clusters;
% nCenters is the number of clusters
deg=zeros(1,nCenters);
u = unique(clusters);
for i=1:nCenters
    [l,~]=anatomicalClusterDistance(zmap, [x y z], u(i), clusTH);
    deg(i)=sum(abs(l>0))/sum(abs(zmap==u(i)));
end
deg=deg';
end