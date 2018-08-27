load multi_disease_DMN_ext_fix DMNdata
code = '76';
load(sprintf('md_split_%s', code), 'train_idx');
train_idx = train_idx + 1;
fulldata = DMNdata';
fulldata = fulldata(:,train_idx);
[nVoxel,nSubj] = size(fulldata);
clear DMNdata
%%
nBoot = 100;
linkageMethod= 'average';
propSamp = 0.67;
nSample = round(nSubj*propSamp);
nBootstrap = 100;
sigma = 0.5;
nu = 3;
%%
for nCenters=100:100:500
    cls = zeros(nVoxel,nBoot);
    cls_consensus = zeros(nVoxel,nBootstrap);
    cls_consensus_filter = zeros(nVoxel,nBootstrap);
    for j=1:nBootstrap
        smpl = sort(randsample(nSubj,nSample,true));
        for i=1:nBoot
            cls(:,i) = kmeans(fulldata(:,smpl), nCenters, 'Distance', 'correlation', 'EmptyAction', 'singleton','MaxIter',200);
        end
        [cls_consensus(:,j),~,M] = consensus_f(cls,nCenters,linkageMethod);
        filter_idx = filter_consensus_matrix(M,sigma,nu);
        cls_consensus_filter(filter_idx,j) = cls_consensus(filter_idx,j); 
        save(sprintf('stability_filter_%d.mat',nCenters), 'cls_consensus','cls_consensus_filter', 'j', ...
            'nCenters','sigma','nu','propSamp','nBoot','nBootstrap')
    end
end






