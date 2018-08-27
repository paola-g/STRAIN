load multi_disease_DMN_ext_fix 
code = '76';
load(sprintf('md_split_%s', code), 'train_idx');
train_idx = train_idx + 1;
fulldata = DMNdata';
fulldata = fulldata(:,train_idx);
[nVoxel,nSubj] = size(fulldata);
clear DMNdata
control_idx = find(DMNlabels(train_idx)==1);
sla_idx = find(DMNlabels(train_idx)==2);
park_idx = find(DMNlabels(train_idx)==3);
%%
nBoot = 100;
linkageMethod= 'average';
propSamp = 0.67;
nSample = round(nSubj*propSamp);
nBootstrap = 100;
sigma = 0.5;
nu = 3;
%%
nCenters=500;



for j=1:nBootstrap
    cls = zeros(nVoxel,nBoot);
    cls_consensus_filter = zeros(nVoxel,1);
    control_smp = randsample(control_idx,round(size(control_idx,2)*propSamp),true);
    sla_smp = randsample(sla_idx,round(size(sla_idx,2)*propSamp),true);
    park_smp = randsample(park_idx,round(size(park_idx,2)*propSamp),true);
    smpl = sort([control_smp, sla_smp, park_smp]);
    mydata = fulldata(:,smpl);
    parfor i=1:nBoot;
        cls(:,i) = kmeans(mydata, nCenters, 'Distance', 'correlation', 'EmptyAction', 'singleton','MaxIter',200);
    end
    [cls_consensus,~,M] = consensus_f(cls,nCenters,linkageMethod);
    filter_idx = filter_consensus_matrix(M,sigma,nu);
    cls_consensus_filter(filter_idx) = cls_consensus(filter_idx); 
    mydata=fulldata';
    [~,ranks_controls,~] = perSubjectRankAbs(mydata(control_smp,filter_idx),cls_consensus(filter_idx),'median');
    [~,ranks_sla,~] = perSubjectRankAbs(mydata(sla_smp,filter_idx),cls_consensus(filter_idx),'median');
    [~,ranks_park,~] = perSubjectRankAbs(mydata(park_smp,filter_idx),cls_consensus(filter_idx),'median');
    save(sprintf('subsamples_filter_boot%d.mat',j), 'cls_consensus','cls_consensus_filter', 'j', ...
        'nCenters','sigma','nu','propSamp','nBoot','nBootstrap', 'filter_idx','ranks_controls', 'ranks_sla', 'ranks_park')
end






