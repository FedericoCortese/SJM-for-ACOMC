source("Utils3_.R")
library(reticulate)
library(zoo)
import("scipy")
source_python('SJ.py')

# 2015XX169 ---------------------------------------------------------------
df2015XX169=read.table("propagation_2015XX169_new_v2.txt",header=T)
unique(df2015XX169$type)

true_states=order_states(df2015XX169$type)

df2015XX169=df2015XX169[,c("a","e","theta","omega")]

df2015XX169=compute_feat(df2015XX169,wdn=10,am1=T)
N=dim(df2015XX169)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2015XX169)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2015XX169)
Lnsat=sat_mod$Lnsat

start_est2015XX169=Sys.time()
est2015XX169 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=df2015XX169,
                                                     Lnsat=Lnsat,
                                                     K=3,
                                                     true_states=true_states),
                                   mc.cores = parallel::detectCores()-1)
end_est2015XX169=Sys.time()
elapsed_est2015XX169=end_est2015XX169-start_est2015XX169
save(est2015XX169,elapsed_est2015XX169,file="est2015XX169.RData")


# 2016CA138 ---------------------------------------------------------------
df2016CA138=read.table("propagation_2016CA138_new_v2.txt",header=T)
unique(df2016CA138$type)

true_states=order_states(df2016CA138$type)

df2016CA138=df2016CA138[,c("a","e","theta","omega")]

df2016CA138=compute_feat(df2016CA138,wdn=10,am1=T)
N=dim(df2016CA138)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2016CA138)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2016CA138)
Lnsat=sat_mod$Lnsat

start_est2016CA138=Sys.time()
est2016CA138 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=df2016CA138,
                                                     Lnsat=Lnsat,
                                                     K=4,
                                                     true_states=true_states),
                                   mc.cores = parallel::detectCores()-1)
end_est2016CA138=Sys.time()
elapsed_est2016CA138=end_est2016CA138-start_est2016CA138
save(est2016CA138,elapsed_est2016CA138,file="est2016CA138.RData")

# 2016CO246 ---------------------------------------------------------------
df2016CO246=read.table("propagation_2016CO246_new_v2.txt",header=T)
unique(df2016CO246$type)

windows()
par(mfrow=c(2,1))
plot(df2016CO246$t,df2016CO246$a,type='l')
plot(df2016CO246$t,df2016CO246$type,type='l',col='red')

true_states=order_states(df2016CO246$type)

df2016CO246=df2016CO246[,c("a","e","theta","omega")]

df2016CO246=compute_feat(df2016CO246,wdn=75,am1=T)
N=dim(df2016CO246)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2016CO246)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2016CO246)
Lnsat=sat_mod$Lnsat

start_est2016CO246=Sys.time()
est2016CO246 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=df2016CO246,
                                                     Lnsat=Lnsat,
                                                     K=4,
                                                     true_states=true_states),
                                   mc.cores = parallel::detectCores()-1)
end_est2016CO246=Sys.time()
elapsed_est2016CO246=end_est2016CO246-start_est2016CO246
save(est2016CO246,elapsed_est2016CO246,file="est2016CO246.RData")

# 2014OL339 ---------------------------------------------------------------
df2014OL339=read.table("propagation_2014OL339_new_v2.txt",header=T)
unique(df2014OL339$type)

true_states=order_states(df2014OL339$type)

df2014OL339=df2014OL339[,c("a","e","theta","omega")]

df2014OL339=compute_feat(df2014OL339,wdn=10,am1=T)
N=dim(df2014OL339)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2014OL339)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2014OL339)
Lnsat=sat_mod$Lnsat

start_est2014OL339=Sys.time()
est2014OL339 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=df2014OL339,
                                                     Lnsat=Lnsat,
                                                     K=4,
                                                     true_states=true_states),
                                   mc.cores = parallel::detectCores()-1)
end_est2014OL339=Sys.time()
elapsed_est2014OL339=end_est2014OL339-start_est2014OL339
save(est2014OL339,elapsed_est2014OL339,file="est2014OL339.RData")

# 2020PN1 -----------------------------------------------------------------
df2020PN1=read.table("propagation_2020PN1_new_v2.txt",header=T)
unique(df2020PN1$type)

true_states=order_states(df2020PN1$type)

df2020PN1=df2020PN1[,c("a","e","theta","omega")]

df2020PN1=compute_feat(df2020PN1,wdn=10,am1=T)
N=dim(df2020PN1)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2020PN1)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

true_states=tail(true_states,N)

sat_mod=SJM_sat(df2020PN1)
Lnsat=sat_mod$Lnsat

start_est2020PN1=Sys.time()
est2020PN1 <- parallel::mclapply(1:nrow(hp),
                                 function(x)
                                   SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   df=df2020PN1,
                                                   Lnsat=Lnsat,
                                                   K=4,
                                                   true_states=true_states),
                                 mc.cores = parallel::detectCores()-1)
end_est2020PN1=Sys.time()
elapsed_est2020PN1=end_est2020PN1-start_est2020PN1
save(est2020PN1,elapsed_est2020PN1,file="est2020PN1.RData")
