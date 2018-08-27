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
%% same experiments on random subsamples
n_samples = 1000;
prop_sample = .67;
ranks_controls = zeros(round(size(control_idx,2)*prop_sample), 405, n_samples);
ranks_sla = zeros(round(size(sla_idx,2)*prop_sample), 405, n_samples);
ranks_park = zeros(round(size(park_idx,2)*prop_sample), 405, n_samples);

for i=1:n_samples
    control_smp = randsample(control_idx,round(size(control_idx,2)*prop_sample),true);
    sla_smp = randsample(sla_idx,round(size(sla_idx,2)*prop_sample),true);
    park_smp = randsample(park_idx,round(size(park_idx,2)*prop_sample),true);
    [~,ranks_controls(:,:,i),~] = perSubjectRankAbs(DMNtest(control_smp,idxs),cl_filt(idxs),'median');
    [~,ranks_sla(:,:,i),~] = perSubjectRankAbs(DMNtest(sla_smp,idxs),cl_filt(idxs),'median');
    [~,ranks_park(:,:,i),~] = perSubjectRankAbs(DMNtest(park_smp,idxs),cl_filt(idxs),'median');
end
save md_ext76_kmeans_f05s3_doc_subsranks ranks_controls ranks_sla ranks_park

%% input topklists anat 
load multi_disease_DMN_ext_fix DMNlabels
load medianvois
control_idx_all = find(DMNlabels==1);
sla_idx_all = find(DMNlabels==2);
park_idx_all = find(DMNlabels==3);
n_samples = 1000;
prop_sample = .67;

ranks_ctrl = zeros(round(size(control_idx_all,2)*prop_sample), 90, n_samples);
ranks_sla = zeros(round(size(sla_idx_all,2)*prop_sample), 90, n_samples);
ranks_park = zeros(round(size(park_idx_all,2)*prop_sample), 90, n_samples);
for i=1:n_samples
    control_smp = randsample(control_idx_all,round(size(control_idx_all,2)*prop_sample),false);
    sla_smp = randsample(sla_idx_all,round(size(sla_idx_all,2)*prop_sample),false);
    park_smp = randsample(park_idx_all,round(size(park_idx_all,2)*prop_sample),false);
    medians_ctrl = medianvois(control_smp, :);
    medians_sla = medianvois(sla_smp, :);
    medians_park = medianvois(park_smp, :);
    [~,ranks_sla(:,:,i)] = sort(abs(medians_sla),2,'descend');
    [~,ranks_park(:,:,i)] = sort(abs(medians_park),2,'descend');
    [~,ranks_ctrl(:,:,i)] = sort(abs(medians_ctrl),2,'descend');
end
save md_anatvoi_subranks_norep ranks_sla ranks_park ranks_ctrl 

%% input topklists consensus
matFiles = dir('C:/Users/Paola/Documents/neuronelab/subsamples_filter_boot*.mat');
n_samples = 1000;
prop_sample = .67;
ranks_controls = zeros(round(size(control_idx,2)*prop_sample), 405, n_samples);
ranks_sla = zeros(round(size(sla_idx,2)*prop_sample), 405, n_samples);
ranks_park = zeros(round(size(park_idx,2)*prop_sample), 405, n_samples);
for k = 1:length(matFiles)
  baseFileName = matFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  storedStructure = load(fullFileName);
  % Now do whatever you need to do with the structure you recalled.
end
