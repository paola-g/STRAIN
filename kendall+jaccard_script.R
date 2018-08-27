library(R.matlab)
library(TopKLists)
setwd('.')

# PERMUTATIONS
# anatvoi
files = list.files(pattern='anatvoi_toprank_d8-8-8v10_perm.*.mat')
nspace = 90
kdists_an = matrix(0,length(files),3)
for (i in 1:length(files)){
  conn <- readMat(files[i])
  controls <- conn$topkcontrols
  sla <- conn$topksla
  park <- conn$topkpark
  n_c = length(controls)
  n_s = length(sla)
  n_p = length(park)
  rank_c = rep(n_c + 1, nspace)
  rank_s = rep(n_s + 1, nspace)
  rank_p = rep(n_p + 1, nspace)
  rank_c[match(controls, 1:nspace)] = 1:n_c
  rank_s[match(sla, 1:nspace)] = 1:n_s
  rank_p[match(park, 1:nspace)] = 1:n_p
  # 1:CtrlSla, 2: CtrlPark, 3: SlaPark
  kdists_an[i,1] = Kendall2Lists.c(rank_c, rank_s, n_c, n_s,nspace)
  kdists_an[i,2] = Kendall2Lists.c(rank_c, rank_p, n_c, n_p,nspace)
  kdists_an[i,3] = Kendall2Lists.c(rank_s, rank_p, n_s, n_p,nspace)
}
conn1 <- readMat('md_anatvoi_abs_median_topranks_d8-8-8v10.mat')
controls <- conn1$topkcontrols
sla <- conn1$topksla
park <- conn1$topkpark
n_c = length(controls)
n_s = length(sla)
n_p = length(park)
rank_c = rep(n_c + 1, nspace)
rank_s = rep(n_s + 1, nspace)
rank_p = rep(n_p + 1, nspace)
rank_c[match(controls, 1:nspace)] = 1:n_c
rank_s[match(sla, 1:nspace)] = 1:n_s
rank_p[match(park, 1:nspace)] = 1:n_p
# 1:CtrlSla, 2: CtrlPark, 3: SlaPark
CtrlSla_an = Kendall2Lists.c(rank_c, rank_s, n_c, n_s,nspace)
CtrlPark_an = Kendall2Lists.c(rank_c, rank_p, n_c, n_p,nspace)
SlaPark_an = Kendall2Lists.c(rank_s, rank_p, n_s, n_p,nspace)
pCtrlSla_an = length(which(kdists_an[,1]>CtrlSla_an))/length(files)
pCtrlPark_an = length(which(kdists_an[,2]>CtrlPark_an))/length(files)
pSlaPark_an = length(which(kdists_an[,3]>SlaPark_an))/length(files)

# clustering
files = list.files(pattern='consensus_toprank_d8-10-8v10_perm.*.mat')
nspace = 405
kdists_cl = matrix(0,length(files),3)
for (i in 1:length(files)){
  conn <- readMat(files[i])
  controls <- conn$topkcontrols
  sla <- conn$topksla
  park <- conn$topkpark
  n_c = length(controls)
  n_s = length(sla)
  n_p = length(park)
  rank_c = rep(n_c + 1, nspace)
  rank_s = rep(n_s + 1, nspace)
  rank_p = rep(n_p + 1, nspace)
  rank_c[match(controls, 1:nspace)] = 1:n_c
  rank_s[match(sla, 1:nspace)] = 1:n_s
  rank_p[match(park, 1:nspace)] = 1:n_p
  # 1:CtrlSla, 2: CtrlPark, 3: SlaPark
  kdists_cl[i,1] = Kendall2Lists.c(rank_c, rank_s, n_c, n_s,nspace)
  kdists_cl[i,2] = Kendall2Lists.c(rank_c, rank_p, n_c, n_p,nspace)
  kdists_cl[i,3] = Kendall2Lists.c(rank_s, rank_p, n_s, n_p,nspace)

}
conn0 <- readMat('md_ext76_kmeans_f05s3_toprank_d8-10-8v10_abs_median.mat')
controls <- conn1$topkcontrols
sla <- conn1$topksla
park <- conn1$topkpark
n_c = length(controls)
n_s = length(sla)
n_p = length(park)
rank_c = rep(n_c + 1, nspace)
rank_s = rep(n_s + 1, nspace)
rank_p = rep(n_p + 1, nspace)
rank_c[match(controls, 1:nspace)] = 1:n_c
rank_s[match(sla, 1:nspace)] = 1:n_s
rank_p[match(park, 1:nspace)] = 1:n_p
# 1:CtrlSla, 2: CtrlPark, 3: SlaPark
CtrlSla_cl = Kendall2Lists.c(rank_c, rank_s, n_c, n_s,nspace)
CtrlPark_cl = Kendall2Lists.c(rank_c, rank_p, n_c, n_p,nspace)
SlaPark_cl = Kendall2Lists.c(rank_s, rank_p, n_s, n_p,nspace)
pCtrlSla_cl = length(which(kdists_cl[,1]>CtrlSla_cl))/length(files)
pCtrlPark_cl = length(which(kdists_cl[,2]>CtrlPark_cl))/length(files)
pSlaPark_cl = length(which(kdists_cl[,3]>SlaPark_cl))/length(files)


writeMat("kendall_perm.mat", CtrlSla_an=CtrlSla_an, CtrlPark_an=CtrlPark_an, 
         SlaPark_an=SlaPark_an, kdists_an=kdists_an, CtrlSla_cl=CtrlSla_cl, 
         CtrlPark_cl=CtrlPark_cl, SlaPark_cl=SlaPark_cl, kdists_cl=kdists_cl, 
         pCtrlSla_an=pCtrlSla_an, pCtrlPark_an=pCtrlPark_an, pSlaPark_an=pSlaPark_an,
         pCtrlSla_cl=pCtrlSla_cl, pCtrlPark_cl=pCtrlPark_cl, pSlaPark_cl=pSlaPark_cl)


# SUBSAMPLES JACCARD
# anatvoi
files = list.files(pattern='anatvoi_toprank_d8-8-8v10_smp_norep.*.mat')
nspace = 90
nfiles = length(files)
jdists_an_smp = matrix(0,nfiles*(nfiles-1)/2,3)

index = 0
for (i in 1:nfiles){
  for (j in 1:nfiles){
    if (i<j){
      print(index)
      conn_i <- readMat(files[i])
      conn_j <- readMat(files[j])
      controls_i <- conn_i$topkcontrols
      sla_i <- conn_i$topksla
      park_i <- conn_i$topkpark
      controls_j <- conn_j$topkcontrols
      sla_j <- conn_j$topksla
      park_j <- conn_j$topkpark

      # 1: Ctrl, 2: Sla, 3: Park
      I <- length(intersect(controls_i,controls_j))
      S <- I/(length(controls_i)+length(controls_j)-I)
      jdists_an_smp[index,1] = S
      I <- length(intersect(sla_i,sla_j))
      S <- I/(length(sla_i)+length(sla_j)-I)
      jdists_an_smp[index,2] = S
      I <- length(intersect(park_i,park_j))
      S <- I/(length(park_i)+length(park_j)-I)
      jdists_an_smp[index,3] = S
      index = index+1
    }
  }
}
writeMat('jdists_an_smp', jdists_an_smp=jdists_an_smp)

# PERM JACCARD
# anatvoi
files = list.files(pattern='anatvoi_toprank_d8-8-8v10_perm.*.mat')
nspace = 90
nfiles = length(files)
jdists_an_perm = matrix(0,nfiles*(nfiles-1)/2,3)

index = 0
for (i in 1:nfiles){
  for (j in 1:nfiles){
    if (i<j){
      print(index)
      conn_i <- readMat(files[i])
      conn_j <- readMat(files[j])
      controls_i <- conn_i$topkcontrols
      sla_i <- conn_i$topksla
      park_i <- conn_i$topkpark
      controls_j <- conn_j$topkcontrols
      sla_j <- conn_j$topksla
      park_j <- conn_j$topkpark
      
      # 1: Ctrl, 2: Sla, 3: Park
      I <- length(intersect(controls_i,controls_j))
      S <- I/(length(controls_i)+length(controls_j)-I)
      jdists_an_perm[index,1] = S
      I <- length(intersect(sla_i,sla_j))
      S <- I/(length(sla_i)+length(sla_j)-I)
      jdists_an_perm[index,2] = S
      I <- length(intersect(park_i,park_j))
      S <- I/(length(park_i)+length(park_j)-I)
      jdists_an_perm[index,3] = S
      index = index+1
    }
  }
}
writeMat('jdists_an_perm', jdists_an_perm=jdists_an_perm)