load multi_disease_DMN DMNdata
load md_split_train
train_idx = train_idx + 1;
fulldata = DMNdata';
fulldata = fulldata(:,train_idx);
[nVoxel,nSubj] = size(fulldata);
clear DMNdata
%%
nCenters = 500;
nBoot = 100;
linkageMethod= 'average';
propSamp = 0.67;
nSample = round(nSubj*propSamp);
nBootstrap = 30;
%%
cls = zeros(nVoxel,nBoot);
cls_consensus = zeros(nVoxel,nBootstrap);
M = zeros(nVoxel,nVoxel);

for j=1:nBootstrap
    smpl = sort(randsample(nSubj,nSample,true));
    for i=1:nBoot
        cls(:,i) = kmeans(fulldata(:,smpl), nCenters, 'Distance', 'correlation', 'EmptyAction', 'singleton');
    end
    [cls_consensus(:,j),~,~] = consensus_f(cls,nCenters,linkageMethod);
end







