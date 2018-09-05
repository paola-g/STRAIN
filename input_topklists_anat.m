%% input topklists
load medianvois % median values computed from AAL parcellation
medians_ctrl = medianvois(control_idx, :);
medians_sla = medianvois(sla_idx, :);
medians_park = medianvois(park_idx, :);
[Y,rank_sla] = sort(abs(medians_sla),2,'descend');
[Y,rank_park] = sort(abs(medians_park),2,'descend');
[Y,rank_ctrl] = sort(abs(medians_ctrl),2,'descend');
save md_anatvoi_abs_median_ranks rank_sla rank_park rank_ctrl

%% on test set
%% input topklists
load medianvois % median values computed from AAL parcellation
load multi_disease_DMN_ext_fix DMNlabels
load md_split_76
train_idx = train_idx+1;
test_idx=test_idx+1;
test_labels = DMNlabels(test_idx);
control_idx = find(test_labels==1);
sla_idx = find(test_labels==2);
park_idx = find(test_labels==3);
medians_ctrl = medianvois(control_idx, :);
medians_sla = medianvois(sla_idx, :);
medians_park = medianvois(park_idx, :);
[Y,rank_sla] = sort(abs(medians_sla),2,'descend');
[Y,rank_park] = sort(abs(medians_park),2,'descend');
[Y,rank_ctrl] = sort(abs(medians_ctrl),2,'descend');
save md_anatvoi_test_ranks rank_sla rank_park rank_ctrl