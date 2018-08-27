library(R.matlab)
library(TopKLists)
library(parallel)
setwd("/mnt/data/pgaldi/topklists/")
source('TopKInference.R')

files = list.files(pattern='subsamples_f.*.mat')

inputs <- 1:length(files)
processInput <- function(j) {
try({
  conn <- readMat(files[j])
  
  ranks_controls <- conn$ranks.controls
  ranks_sla <- conn$ranks.sla
  ranks_park <- conn$ranks.park
  data_controls = as.data.frame(t(ranks_controls))
  data_sla = as.data.frame(t(ranks_sla))
  data_park = as.data.frame(t(ranks_park))
  calculate.maxK(data_controls, 17, 8, 10, 50) -> results_controls
  calculate.maxK(data_sla, 19, 10, 10, 50) -> results_sla
  calculate.maxK(data_park, 17, 8, 10, 50) -> results_park
  writeMat(paste("consensus_toprank_d8-10-8v10_smp",j,".mat",sep=''), topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)
})
}


results = mclapply(inputs, processInput, mc.cores = 8)



# ANATOMIC REGIONS (all data)
# absolute value
conn <- readMat("md_anatvoi_subranks.mat")
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
  calculate.maxK(data_controls, 25, 8, 10, 50) -> results_controls
  calculate.maxK(data_sla, 27, 8, 10, 50) -> results_sla
  calculate.maxK(data_park, 25, 8, 10, 50) -> results_park
  writeMat(paste("anatvoi_toprank_d8-8-8v10_smp",j,".mat",sep=''), topkcontrols=results_controls$topkspace, topksla=results_sla$topkspace, topkpark=results_park$topkspace)
  }
)}
results = mclapply(inputs, processInput, mc.cores = 16)
