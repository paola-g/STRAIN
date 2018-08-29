library(R.matlab)
library(TopKLists)
library(parallel)
setwd(".")

# output of input_topklists_perm.m
conn <- readMat("md_ext76_kmeans_f05s3_doc_permranks.mat") 

ranks_controls <- conn$ranks.controls
ranks_sla <- conn$ranks.sla
ranks_park <- conn$ranks.park

inputs <- 1:1000
processInput <- function(j) 
{
data_controls = as.data.frame(t(ranks_controls[,,j]))
data_sla = as.data.frame(t(ranks_sla[,,j]))
data_park = as.data.frame(t(ranks_park[,,j]))
calculate.maxK(data_controls, 11, 8, 10, 50) -> results_controls
calculate.maxK(data_sla, 12, 10, 10, 50) -> results_sla
calculate.maxK(data_park, 11, 8, 10, 50) -> results_park
writeMat(paste("consensus_toprank_d8-10-8v10_perm",j,".mat",sep=''), topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)
}

numCores <- detectCores()

results = mclapply(inputs, processInput, mc.cores = 20)



# ANATOMIC REGIONS 
# output of input_topklists_perm.m
conn <- readMat("md_anatvoi_permranks.mat")
ranks_controls <- conn$ranks.ctrl
ranks_sla <- conn$ranks.sla
ranks_park <- conn$ranks.park
inputs <- 1:1000
processInput <- function(j) 
{
  try(
  {
  data_controls = as.data.frame(t(ranks_controls[,,j]))
  data_sla = as.data.frame(t(ranks_sla[,,j]))
  data_park = as.data.frame(t(ranks_park[,,j]))
  calculate.maxK(data_controls, 37, 8, 10, 50) -> results_controls
  calculate.maxK(data_sla, 41, 8, 10, 50) -> results_sla
  calculate.maxK(data_park, 37, 8, 10, 50) -> results_park
  writeMat(paste("anatvoi_toprank_d8-8-8v10_perm",j,".mat",sep=''), topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)
  }
)}
results = mclapply(inputs, processInput, mc.cores = 16)
