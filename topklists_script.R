library(R.matlab)
library(TopKLists)
setwd("C:/Users/Paola/Dropbox/neuronelab/")



# absolute value
conn <- readMat("md_ext76_kmeans_f05s3_doc_ranks_abs_median.mat")
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
writeMat("md_ext76_kmeans_f05s3_toprank_d8-10-8v10_abs_median.mat", topkspace=results_all$topkspace, topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)


# descending order
conn <- readMat("md_ext76_kmeans_f05s3_doc_ranks_desc_median.mat")
ranks_all <- conn$ranks.all
ranks_controls <- conn$ranks.controls
ranks_sla <- conn$ranks.sla
ranks_park <- conn$ranks.park
data_all = as.data.frame(t(ranks_all))
data_controls = as.data.frame(t(ranks_controls))
data_sla = as.data.frame(t(ranks_sla))
data_park = as.data.frame(t(ranks_park))
# delta plots
a=deltaplot(data_all, deltas = 1:40, subset.lists=50, directory='topklists/median_desc/all')
a=deltaplot(data_controls, deltas = 1:40, subset.lists=50, directory='topklists/median_desc/controls')
a=deltaplot(data_sla, deltas = 1:40, subset.lists=50, directory='topklists/median_desc/sla')
a=deltaplot(data_park, deltas = 1:40, subset.lists=50, directory='topklists/median_desc/park')
calculate.maxK(data_all, 34, 7, 10, 50) -> results_all
calculate.maxK(data_controls, 11, 7, 10, 50) -> results_controls
calculate.maxK(data_sla, 12, 7, 10, 50) -> results_sla
calculate.maxK(data_park, 11, 9, 10, 50) -> results_park
writeMat("md_ext76_kmeans_f05s3_toprank_d7-7-9v10_desc_median.mat", topkspace=results_all$topkspace, topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)



# ascending order
conn <- readMat("md_ext76_kmeans_f05s3_doc_ranks_asc_median.mat")
ranks_all <- conn$ranks.all
ranks_controls <- conn$ranks.controls
ranks_sla <- conn$ranks.sla
ranks_park <- conn$ranks.park
data_all = as.data.frame(t(ranks_all))
data_controls = as.data.frame(t(ranks_controls))
data_sla = as.data.frame(t(ranks_sla))
data_park = as.data.frame(t(ranks_park))
# delta plots
a=deltaplot(data_all, deltas = 1:40, subset.lists=50, directory='topklists/median_asc/all')
a=deltaplot(data_controls, deltas = 1:40, subset.lists=50, directory='topklists/median_asc/controls')
a=deltaplot(data_sla, deltas = 1:40, subset.lists=50, directory='topklists/median_asc/sla')
a=deltaplot(data_park, deltas = 1:40, subset.lists=50, directory='topklists/median_asc/park')
calculate.maxK(data_all, 34, 7, 10, 50) -> results_all
calculate.maxK(data_controls, 11, 10, 10, 50) -> results_controls
calculate.maxK(data_sla, 12, 10, 10, 50) -> results_sla
calculate.maxK(data_park, 11, 11, 10, 50) -> results_park
writeMat("md_ext76_kmeans_f05s3_toprank_d10-10-11v10_asc_median.mat", topkspace=results_all$topkspace, topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)



# ANATOMIC REGIONS (test set)
conn <- readMat("ranks_median_sla_test.mat")
myranks <- conn$I
mydata = as.data.frame(t(myranks))
myresults <- calculate.maxK(mydata, 12, 15, 10, 50)
writeMat("medianvois_sla_topklists.mat",topksla=myresults$topkspace)

conn <- readMat("ranks_median_park_test.mat")
myranks <- conn$I
mydata = as.data.frame(t(myranks))
myresults <- calculate.maxK(mydata, 11, 5, 10, 50)
writeMat("medianvois_park_topklists.mat",topkpark=myresults$topkspace)

conn <- readMat("ranks_median_control_test.mat")
myranks <- conn$I
mydata = as.data.frame(t(myranks))
myresults <- calculate.maxK(mydata, 11, 5, 10, 50)
writeMat("medianvois_control_topklists.mat",topkcontrol=myresults$topkspace)


# ANATOMIC REGIONS (all data)
# absolute value
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

# descending value
conn <- readMat("md_anatvoi_desc_median_ranks.mat")
ranks_controls <- conn$A
ranks_sla <- conn$rank.sla
ranks_park <- conn$rank.park
data_controls = as.data.frame(t(ranks_controls))
data_sla = as.data.frame(t(ranks_sla))
data_park = as.data.frame(t(ranks_park))
# delta plots
a=deltaplot(data_controls, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_desc/controls')
a=deltaplot(data_sla, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_desc/sla')
a=deltaplot(data_park, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_desc/park')
calculate.maxK(data_controls, 36, 7, 10, 50) -> results_controls
calculate.maxK(data_sla, 41, 8, 10, 50) -> results_sla
calculate.maxK(data_park, 37, 8, 10, 50) -> results_park
writeMat("md_anatvoi_desc_median_topranks_d7-8-8v10.mat", topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)

# ascending value
conn <- readMat("md_anatvoi_asc_median_ranks.mat")
ranks_controls <- conn$A
ranks_sla <- conn$rank.sla
ranks_park <- conn$rank.park
data_controls = as.data.frame(t(ranks_controls))
data_sla = as.data.frame(t(ranks_sla))
data_park = as.data.frame(t(ranks_park))
# delta plotscompu
a=deltaplot(data_controls, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_asc/controls')
a=deltaplot(data_sla, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_asc/sla')
a=deltaplot(data_park, deltas = 1:40, subset.lists=90, directory='topklists/anat/median_asc/park')
calculate.maxK(data_controls, 36, 6, 10, 50) -> results_controls
calculate.maxK(data_sla, 41, 7, 10, 50) -> results_sla
calculate.maxK(data_park, 37, 8, 10, 50) -> results_park
writeMat("md_anatvoi_asc_median_topranks_d6-7-8v10.mat", topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)