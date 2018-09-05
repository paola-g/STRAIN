n_cls = 100;
n_k = 5;
scores = zeros(n_cls*(n_cls-1)/2,n_k);
scores_filter = zeros(n_cls*(n_cls-1)/2,n_k);
i_k=0;
for i=100:100:500
    load(sprintf('stability_filter_%d.mat', i)) % output of md_kmeans_filter_stability.m
    idx=0;
    i_k=i_k+1;
    for j=1:n_cls
       for k=1:n_cls
          if j<k 
              idx = idx+1;            
              cls_consensus(cls_consensus(:,j)==0,j) = grp2idx(1000:999+length(find(cls_consensus(:,j)==0)));
              cls_consensus_filter(cls_consensus_filter(:,j)==0,j) = grp2idx(1000:999+length(find(cls_consensus_filter(:,j)==0)));
              cls_consensus(cls_consensus(:,k)==0,k) = grp2idx(1000:999+length(find(cls_consensus(:,k)==0)));
              cls_consensus_filter(cls_consensus_filter(:,k)==0,k) = grp2idx(1000:999+length(find(cls_consensus_filter(:,k)==0)));
              scores(idx, i_k) = ami(cls_consensus(:,j), cls_consensus(:,k));
              scores_filter(idx, i_k) = ami(cls_consensus_filter(:,j), cls_consensus_filter(:,k));

          end
       end
    end
    fprintf('Mean NMI score for %d clusters: %.2f (+/- %.2f) (with filter %.2f (+/- %.2f))\n', i, mean(scores(:,i_k)), std(scores(:,i_k)),...
        mean(scores_filter(:,i_k)), std(scores_filter(:,i_k)))
end


%% 

for i=1:6
fprintf('Mean NMI score for %d clusters: %.2f (+/- %.2f) (with filter [on %.    f%% voxels] %.2f (+/- %.2f))\n', ...
    n_cls(i), mean(scores(:,i)),std(scores(:,i)), 100*mean(nvoxels(:,i))/13423, mean(scores_filter(:,i)), std(scores_filter(:,i)))
end

%%
load ami_results
n_cls = 100;
n_k = 6;
nfiltered = zeros(n_cls,n_k);
i_k=0;
for i=[90,100:100:500]
    i_k=i_k+1;
    load(sprintf('stability_filter_%d.mat', i))
    nfiltered(:,i_k) = sum(cls_consensus_filter~=0);
    fprintf('Mean NMI score for %d clusters: %.2f (+/- %.2f) (with filter %.2f (+/- %.2f)) [filtered %.f%% voxels intersection %.f%%] \n', ...
    i, mean(scores(:,i_k)),std(scores(:,i_k)), mean(scores_filter(:,i_k)), std(scores_filter(:,i_k)), 100*mean(nfiltered(:,i_k))/13423, 100*mean(nvoxels(:,i_k))/13423)
end

%%
load ami_results
n_cls = 100;
n_k = 6;
norm = zeros(n_cls,n_k);
i_k=0;
idx=0;
for i=[90,100:100:500]
    idx=0;
    i_k=i_k+1;
    load(sprintf('stability_filter_%d.mat', i))
    for j=1:n_cls
       for k=1:n_cls
          if j<k 
              idx = idx+1;
              norm(idx, i_k) = mean([sum(cls_consensus_filter(:,j)~=0),sum(cls_consensus_filter(:,k)~=0)]);
          end
       end
    end
end

%%
load ami_results
n_cls = 100;
n_k = 6;
nfiltered = zeros(n_cls,n_k);
i_k=0;
mean_nmi = zeros(n_k,2);
std_nmi = zeros(n_k,2);
perc_vox = zeros(n_k,2);
std_perc = zeros(n_k,2);
for i=[90,100:100:500]
    i_k=i_k+1;
    load(sprintf('stability_filter_%d.mat', i))
    nfiltered(:,i_k) = sum(cls_consensus_filter~=0);
    normed_perc = nvoxels./norm;
    fprintf('Mean NMI score for %d clusters: %.2f (+/- %.2f) (with filter %.2f (+/- %.2f)) [filtered %.f%% (+/- %.2f) voxels intersection %.f%% (+/- %.2f)] \n', ...
    i, mean(scores(:,i_k)),std(scores(:,i_k)), mean(scores_filter(:,i_k)), std(scores_filter(:,i_k)), ...
    100*mean(nfiltered(:,i_k))/13423, std(nfiltered(:,i_k))/13423, 100*mean(normed_perc(:,i_k)), std(normed_perc(:,i_k)))
    mean_nmi(i_k,1) =  mean(scores(:,i_k));
    mean_nmi(i_k,2) =  mean(scores_filter(:,i_k));
    std_nmi(i_k,1) =  std(scores(:,i_k));
    std_nmi(i_k,2) =  std(scores_filter(:,i_k));
    perc_vox(i_k,1) = mean(nfiltered(:,i_k))/13423;
    perc_vox(i_k,2) = mean(normed_perc(:,i_k));
    std_perc(i_k,1) = std(nfiltered(:,i_k))/13423;
    std_perc(i_k,2) = std(normed_perc(:,i_k));
end

%%
addpath(genpath('C:\\Users\\Paola\\Documents\\neuronelab\\NeuroElf_v09c'))
addpath(genpath('C:\\Users\\Paola\\Documents\\neuronelab'))
n_cls = 100;
n_k = 6;
i_k=0;
idx=0;

for i=[90,100:100:500]
    i_k=i_k+1;
    load(sprintf('stability_filter_%d.mat', i))
    cls_doc = zeros(13423,100);
    for j=1:n_cls
              deg_j = calcDOCresample(cls_consensus_filter(cls_consensus_filter(:,j)~=0,j), 10, cls_consensus_filter(:,j)~=0);
              [idxs_j, ~] = filterDegree(cls_consensus_filter(:,j), deg_j, 0.5);
              cls_doc(idxs_j, j) = cls_consensus_filter(idxs_j,j);
    end
    save(sprintf('stability_doc_%d.mat', i), 'cls_doc')
end

%%
n_cls = 100;
n_k = 6;
norm = zeros(n_cls,n_k);
norm_doc = zeros(n_cls,n_k);
i_k=0;
idx=0;
for i=[90,100:100:500]
    idx=0;
    i_k=i_k+1;
    load(sprintf('stability_filter_%d.mat', i))
    load(sprintf('stability_doc_%d.mat', i))
    for j=1:n_cls
       for k=1:n_cls
          if j<k 
              idx = idx+1;
              norm(idx, i_k) = mean([sum(cls_consensus_filter(:,j)~=0),sum(cls_consensus_filter(:,k)~=0)]);
              norm_doc(idx, i_k) = mean([sum(cls_doc(:,j)~=0),sum(cls_doc(:,k)~=0)]);
          end
       end
    end
end
%%
load ami_results
load ami_doc_3
n_cls = 100;
n_k = 6;
nfiltered = zeros(n_cls,n_k);
nfiltered_doc = zeros(n_cls,n_k);
i_k=0;
for i=[90,100:100:500]
    i_k=i_k+1;
    load(sprintf('stability_filter_%d.mat', i))
    load(sprintf('stability_doc_%d.mat', i))
    nfiltered(:,i_k) = sum(cls_consensus_filter~=0);
    nfiltered_doc(:,i_k) = sum(cls_doc~=0);
    normed_perc = nvoxels./norm;
    normed_perc_doc = nvoxels_doc./norm_doc;
    fprintf('Mean NMI score for %d clusters: %.2f (+/- %.2f)\nWith filter %.2f (+/- %.2f) [filtered %.f%% voxels intersection %.f%%] \n', ...
    i, mean(scores(:,i_k)),std(scores(:,i_k)), mean(scores_filter(:,i_k)), std(scores_filter(:,i_k)), 100*mean(nfiltered(:,i_k))/13423, 100*mean(normed_perc(:,i_k)))
    fprintf('With DoC %.2f (+/- %.2f) [filtered %.f%% voxels intersection %.f%%]\n\n', mean(scores_doc(:,i_k)),std(scores_doc(:,i_k)),100*mean(nfiltered_doc(:,i_k))/13423, 100*mean(normed_perc_doc(:,i_k)))
end

%%
n_cls = 100;
n_k = 6;
nclusters = zeros(n_cls,n_k);
nclusters_doc = zeros(n_cls,n_k);
i_k=0;
for i=[90,100:100:500]
    idx=0;
    i_k=i_k+1;
    load(sprintf('stability_filter_%d.mat', i))
    load(sprintf('stability_doc_%d.mat', i))
    for j=1:n_cls
        idx = idx+1;
        nclusters(idx,i_k) = length(unique(cls_consensus_filter(:,j)))-1;
        nclusters_doc(idx,i_k) = length(unique(cls_doc(:,j)))-1;
    end
    
end
