source("Utils3_.R")
library(reticulate)
library(zoo)
import("scipy")
source_python('SJ.py')

# 164207 ------------------------------------------------------------------

df164207=read.table("propagation_164207_new_v2.txt",header=T)
true_states=order_states(df164207$type)

df164207=df164207[,c("a","e","theta","omega")]


df164207=compute_feat(df164207)
N=dim(df164207)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df164207)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df164207)
Lnsat=sat_mod$Lnsat

start_est164207=Sys.time()
est164207 <- parallel::mclapply(1:nrow(hp),
                                      function(x)
                                        SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                        kappa=hp[x,]$kappa,
                                                        df=df164207,
                                                        Lnsat=Lnsat,
                                                        true_states=true_states),
                                      mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(reticulate)
#   library(pdfCluster)
#   library(ggplot2)
#   library(xtable)
#   library(dplyr))
# est164207 <- clusterApply(cl,
#                     1:nrow(hp),
#                     function(x)
#                       SJM_lambdakappa(lambda=hp[x,]$lambda,
#                                       kappa=hp[x,]$kappa,
#                                       df=df164207,
#                                       Lnsat=Lnsat,
#                                       true_states=true_states)
# )
# stopCluster(cl)

end_est164207=Sys.time()
elapsed_est164207=end_est164207-start_est164207
save(est164207,elapsed_est164207,file="est164207.RData")


# 2001GO2 -------------------------------------------------------------------------

df2001GO2=read.table("propagation_2001GO2_new_v2.txt",header=T)
true_states=order_states(df2001GO2$type)

df2001GO2=df2001GO2[,c("a","e","theta","omega")]

df2001GO2=compute_feat(df2001GO2)
N=dim(df2001GO2)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2001GO2)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2001GO2)
Lnsat=sat_mod$Lnsat

start_est2001GO2=Sys.time()
est2001GO2 <- parallel::mclapply(1:nrow(hp),
                                 function(x)
                                   SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   df=df2001GO2,
                                                   Lnsat=Lnsat,
                                                   true_states=true_states),
                                 mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(reticulate)
#   library(pdfCluster)
#   library(ggplot2)
#   library(xtable)
#   library(dplyr))
# est2001GO2 <- clusterApply(cl,
#                     1:nrow(hp),
#                     function(x)
#                       SJM_lambdakappa(lambda=hp[x,]$lambda,
#                                       kappa=hp[x,]$kappa,
#                                       df=df2001GO2,
#                                       Lnsat=Lnsat,
#                                       true_states=true_states)
# )
# stopCluster(cl)

end_est2001GO2=Sys.time()
elapsed_est2001GO2=end_est2001GO2-start_est2001GO2
save(est2001GO2,elapsed_est2001GO2,file="est2001GO2.RData")

# 2002AA29 --------------------------------------------------------------------

df2002AA29=read.table("propagation_2002AA29_new_v2.txt",header=T)
true_states=order_states(df2002AA29$type)

df2002AA29=df2002AA29[,c("a","e","theta","omega")]

df2002AA29=compute_feat(df2002AA29)
N=dim(df2002AA29)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2002AA29)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2002AA29)
Lnsat=sat_mod$Lnsat

start_est2002AA29=Sys.time()
est2002AA29 <- parallel::mclapply(1:nrow(hp),
                                  function(x)
                                    SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    df=df2002AA29,
                                                    Lnsat=Lnsat,
                                                    true_states=true_states),
                                  mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(reticulate)
#   library(pdfCluster)
#   library(ggplot2)
#   library(xtable)
#   library(dplyr))
# est2002AA29 <- clusterApply(cl,
#                     1:nrow(hp),
#                     function(x)
#                       SJM_lambdakappa(lambda=hp[x,]$lambda,
#                                       kappa=hp[x,]$kappa,
#                                       df=df2002AA29,
#                                       Lnsat=Lnsat,
#                                       true_states=true_states)
# )
# stopCluster(cl)

end_est2002AA29=Sys.time()
elapsed_est2002AA29=end_est2002AA29-start_est2002AA29
save(est2002AA29,elapsed_est2002AA29,file="est2002AA29.RData")


# 2015SO2 --------------------------------------------------------------------

df2015SO2=read.table("propagation_2015SO2_new_v2.txt",header=T)
true_states=order_states(df2015SO2$type)

df2015SO2=df2015SO2[,c("a","e","theta","omega")]

df2015SO2=compute_feat(df2015SO2)
N=dim(df2015SO2)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2015SO2)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2015SO2)
Lnsat=sat_mod$Lnsat

start_est2015SO2=Sys.time()
est2015SO2 <- parallel::mclapply(1:nrow(hp),
                                 function(x)
                                   SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   df=df2015SO2,
                                                   Lnsat=Lnsat,
                                                   true_states=true_states),
                                 mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(reticulate)
#   library(pdfCluster)
#   library(ggplot2)
#   library(xtable)
#   library(dplyr))
# est2015SO2 <- clusterApply(cl,
#                     1:nrow(hp),
#                     function(x)
#                       SJM_lambdakappa(lambda=hp[x,]$lambda,
#                                       kappa=hp[x,]$kappa,
#                                       df=df2015SO2,
#                                       Lnsat=Lnsat,
#                                       true_states=true_states)
# )
# stopCluster(cl)

end_est2015SO2=Sys.time()
elapsed_est2015SO2=end_est2015SO2-start_est2015SO2
save(est2015SO2,elapsed_est2015SO2,file="est2015SO2.RData")


# 2016HO3 -----------------------------------------------------------------

df2016HO3=read.table("propagation_2016HO3_new_v2.txt",header=T)
true_states=order_states(df2016HO3$type)

df2016HO3=df2016HO3[,c("a","e","theta","omega")]

df2016HO3=compute_feat(df2016HO3)
N=dim(df2016HO3)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2016HO3)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2016HO3)
Lnsat=sat_mod$Lnsat

start_est2016HO3=Sys.time()
est2016HO3 <- parallel::mclapply(1:nrow(hp),
                                 function(x)
                                   SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   df=df2016HO3,
                                                   Lnsat=Lnsat,
                                                   true_states=true_states),
                                 mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(reticulate)
#   library(pdfCluster)
#   library(ggplot2)
#   library(xtable)
#   library(dplyr))
# est2016HO3 <- clusterApply(cl,
#                     1:nrow(hp),
#                     function(x)
#                       SJM_lambdakappa(lambda=hp[x,]$lambda,
#                                       kappa=hp[x,]$kappa,
#                                       df=df2016HO3,
#                                       Lnsat=Lnsat,
#                                       true_states=true_states)
# )
# stopCluster(cl)

end_est2016HO3=Sys.time()
elapsed_est2016HO3=end_est2016HO3-start_est2016HO3
save(est2016HO3,elapsed_est2016HO3,file="est2016HO3.RData")


# 2019GM1 -----------------------------------------------------------------

df2019GM1=read.table("propagation_2019GM1_new_v2.txt",header=T)
true_states=order_states(df2019GM1$type)

df2019GM1=df2019GM1[,c("a","e","theta","omega")]

df2019GM1=compute_feat(df2019GM1)
N=dim(df2019GM1)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2019GM1)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2019GM1)
Lnsat=sat_mod$Lnsat

start_est2019GM1=Sys.time()
est2019GM1 <- parallel::mclapply(1:nrow(hp),
                                 function(x)
                                   SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   df=df2019GM1,
                                                   Lnsat=Lnsat,
                                                   true_states=true_states),
                                 mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(reticulate)
#   library(pdfCluster)
#   library(ggplot2)
#   library(xtable)
#   library(dplyr))
# est2019GM1 <- clusterApply(cl,
#                     1:nrow(hp),
#                     function(x)
#                       SJM_lambdakappa(lambda=hp[x,]$lambda,
#                                       kappa=hp[x,]$kappa,
#                                       df=df2019GM1,
#                                       Lnsat=Lnsat,
#                                       true_states=true_states)
# )
# stopCluster(cl)

end_est2019GM1=Sys.time()
elapsed_est2019GM1=end_est2019GM1-start_est2019GM1
save(est2019GM1,elapsed_est2019GM1,file="est2019GM1.RData")


# 2020PP1 -----------------------------------------------------------------

df2020PP1=read.table("propagation_2020PP1_new_v2.txt",header=T)
true_states=order_states(df2020PP1$type)

df2020PP1=df2020PP1[,c("a","e","theta","omega")]

df2020PP1=compute_feat(df2020PP1)
N=dim(df2020PP1)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2020PP1)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2020PP1)
Lnsat=sat_mod$Lnsat

start_est2020PP1=Sys.time()
est2020PP1 <- parallel::mclapply(1:nrow(hp),
                                 function(x)
                                   SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   df=df2020PP1,
                                                   Lnsat=Lnsat,
                                                   true_states=true_states),
                                 mc.cores = parallel::detectCores()-1)

# cl<-makeCluster(parallel::detectCores(),type="SOCK")
# parallel::clusterExport(cl,ls())
# parallel::clusterEvalQ(cl,{
#   library(reticulate)
#   library(pdfCluster)
#   library(ggplot2)
#   library(xtable)
#   library(dplyr))
# est2020PP1 <- clusterApply(cl,
#                     1:nrow(hp),
#                     function(x)
#                       SJM_lambdakappa(lambda=hp[x,]$lambda,
#                                       kappa=hp[x,]$kappa,
#                                       df=df2020PP1,
#                                       Lnsat=Lnsat,
#                                       true_states=true_states)
# )
# stopCluster(cl)

end_est2020PP1=Sys.time()
elapsed_est2020PP1=end_est2020PP1-start_est2020PP1
save(est2020PP1,elapsed_est2020PP1,file="est2020PP1.RData")
