load md_ext76_kmeans_filters_fix filters
load multi_disease_DMN_ext_fix
load md_split_76
train_idx = train_idx+1;
test_idx=test_idx+1;
load md_ext76_kmeans_sol_consensus_fix
filt = filters{1,3};
cl_all = cl;
cl_filt = cl(filt);
[deg] = calcDOCresample(cl_filt, 10,filt);
[idxs, ~] = filterDegree(grp2idx(cl_filt), deg, .5);
DMNtest = DMNdata(test_idx,:);
test_labels = DMNlabels(test_idx);
control_idx = find(test_labels==1);
sla_idx = find(test_labels==2);
park_idx = find(test_labels==3);
%% permutation testing
n_perm = 2000;
perm_idx = zeros(size(test_labels,2), n_perm);
ranks_controls = zeros(size(control_idx,2), 405, n_perm);
ranks_sla = zeros(size(sla_idx,2), 405, n_perm);
ranks_park = zeros(size(park_idx,2), 405, n_perm);

for i=1:n_perm
    perm_idx(:,i) = randperm(size(test_labels,2));
    perm_labels = test_labels(perm_idx(:,i));
    control_idx = find(perm_labels==1);
    sla_idx = find(perm_labels==2);
    park_idx = find(perm_labels==3);
    [~,ranks_controls(:,:,i),~] = perSubjectRankAbs(DMNtest(control_idx,idxs),cl_filt(idxs),'median');
    [~,ranks_sla(:,:,i),~] = perSubjectRankAbs(DMNtest(sla_idx,idxs),cl_filt(idxs),'median');
    [~,ranks_park(:,:,i),~] = perSubjectRankAbs(DMNtest(park_idx,idxs),cl_filt(idxs),'median');
end
save md_ext76_kmeans_f05s3_doc_permranks2000 ranks_controls ranks_sla ranks_park perm_idx

%% input topklists anat 
load multi_disease_DMN_ext_fix DMNlabels
load medianvois
control_idx_all = find(DMNlabels==1);
sla_idx_all = find(DMNlabels==2);
park_idx_all = find(DMNlabels==3);
n_perm = 2000;
perm_idx = zeros(size(DMNlabels,2), n_perm);
ranks_ctrl = zeros(size(control_idx_all,2), 90, n_perm);
ranks_sla = zeros(size(sla_idx_all,2), 90, n_perm);
ranks_park = zeros(size(park_idx_all,2), 90, n_perm);
for i=1:n_perm
    perm_idx(:,i) = randperm(size(DMNlabels,2));
    perm_labels = DMNlabels(perm_idx(:,i));
    control_idx_all = find(perm_labels==1);
    sla_idx_all = find(perm_labels==2);
    park_idx_all = find(perm_labels==3);  
    medians_ctrl = medianvois(control_idx_all, :);
    medians_sla = medianvois(sla_idx_all, :);
    medians_park = medianvois(park_idx_all, :);
    [~,ranks_sla(:,:,i)] = sort(abs(medians_sla),2,'descend');
    [~,ranks_park(:,:,i)] = sort(abs(medians_park),2,'descend');
    [~,ranks_ctrl(:,:,i)] = sort(abs(medians_ctrl),2,'descend');
end
save md_anatvoi_permranks2000 ranks_sla ranks_park ranks_ctrl perm_idx
