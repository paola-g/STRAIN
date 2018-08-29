library(R.matlab)
library(TopKLists)
setwd(".")

# output of input_topklists.m
conn <- readMat("md_kmeans_consensus_filter_ranks_abs_median.mat")
ranks_all <- conn$ranks.all
ranks_controls <- conn$ranks.controls
ranks_sla <- conn$ranks.sla
ranks_park <- conn$ranks.park
data_all = as.data.frame(t(ranks_all))
data_controls = as.data.frame(t(ranks_controls))
data_sla = as.data.frame(t(ranks_sla))
data_park = as.data.frame(t(ranks_park))
# delta plots
a=deltaplot(data_all, deltas = 1:30, subset.lists=50, directory='topklists/median_abs/all')
a=deltaplot(data_controls, deltas = 1:40, subset.lists=50, directory='topklists/median_abs/controls')
a=deltaplot(data_sla, deltas = 1:40, subset.lists=50, directory='topklists/median_abs/sla')
a=deltaplot(data_park, deltas = 1:40, subset.lists=50, directory='topklists/median_abs/park')
calculate.maxK(data_all, 34, 7, 10, 50) -> results_all
calculate.maxK(data_controls, 11, 8, 10, 50) -> results_controls
calculate.maxK(data_sla, 12, 10, 10, 50) -> results_sla
calculate.maxK(data_park, 11, 8, 10, 50) -> results_park
writeMat("md_cluster_abs_median_topranks_d8-10-8v10.mat", topkspace=results_all$topkspace, topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)


# ANATOMIC REGIONS
# output of anat_parcellation.m
conn <- readMat("md_anatvoi_abs_median_ranks.mat")
ranks_controls <- conn$A
ranks_sla <- conn$rank.sla
ranks_park <- conn$rank.park
data_controls = as.data.frame(t(ranks_controls))
data_sla = as.data.frame(t(ranks_sla))
data_park = as.data.frame(t(ranks_park))
# delta plots
a=deltaplot(data_controls, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_abs/controls')
a=deltaplot(data_sla, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_abs/sla')
a=deltaplot(data_park, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_abs/park')
calculate.maxK(data_controls, 36, 8, 10, 50) -> results_controls
calculate.maxK(data_sla, 41, 8, 10, 50) -> results_sla
calculate.maxK(data_park, 37, 8, 10, 50) -> results_park
writeMat("md_anatvoi_abs_median_topranks_d8-8-8v10.mat", topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)