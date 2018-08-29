load md_filters filters
load multi_disease_DMN
load md_split_train
train_idx = train_idx+1;
test_idx=test_idx+1;
load md_kmeans_sol_consensus
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
%% abs median
[~,ranks_all,~] = perSubjectRankAbs(DMNtest(:,idxs),cl_filt(idxs),'median');
[~,ranks_controls,~] = perSubjectRankAbs(DMNtest(control_idx,idxs),cl_filt(idxs),'median');
[~,ranks_sla,~] = perSubjectRankAbs(DMNtest(sla_idx,idxs),cl_filt(idxs),'median');
[~,ranks_park,~] = perSubjectRankAbs(DMNtest(park_idx,idxs),cl_filt(idxs),'median');
save md_kmeans_consensus_filter_ranks_abs_median ranks_all ranks_controls ranks_sla ranks_park
