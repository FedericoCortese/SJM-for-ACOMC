source("Utils3_.R")
library(reticulate)
library(zoo)
library(ggplot2)
library(ggpubr)
library(dplyr)
import("scipy")
source_python('SJ.py')

# 2015XX169 ---------------------------------------------------------------
df2015XX169=read.table("propagation_2015XX169_new_v2.txt",header=T)
p2015XX169=plot_real4(df2015XX169)

Ptrue2015XX169=ggarrange(p2015XX169$Pa,
                         p2015XX169$Pe,
                         p2015XX169$Ptheta,
                         p2015XX169$Pomega,
                         nrow=4
                         #,main="164207"
)
Ptrue2015XX169

true_states=order_states(df2015XX169$type)
length(unique(true_states))

df2015XX169=df2015XX169[,c("a","e","theta","omega")]

df2015XX169=compute_feat(df2015XX169)
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

modsel2015XX169=data.frame(hp,FTIC=unlist(lapply(est2015XX169,function(x)x$FTIC)),
                           BAC=unlist(lapply(est2015XX169,function(x)x$BAC))
)
sel=1
df2015XX169=read.table("propagation_2015XX169_new_v2.txt",header=T)
df2015XX169=df2015XX169[,c("a","e","theta","omega")]

df2015XX169=tail(df2015XX169,N)
df2015XX169_b=data.frame(df2015XX169,type=est2015XX169[[sel]]$est_states)
p2015XX169_b=plot_real4(df2015XX169_b)

Pest2015XX169_b=ggarrange(p2015XX169_b$Pa,
                          p2015XX169_b$Pe,
                          p2015XX169_b$Ptheta,
                          p2015XX169_b$Pomega,
                          nrow=4
                          #,main="2015XX169"
)
Pest2015XX169_b

gridExtra::grid.arrange(Ptrue2015XX169,Pest2015XX169_b)

# 2016CA138 ---------------------------------------------------------------
df2016CA138=read.table("propagation_2016CA138_new_v2.txt",header=T)
unique(df2016CA138$type)

p2016CA138=plot_real4(df2016CA138)

Ptrue2016CA138=ggarrange(p2016CA138$Pa,
                         p2016CA138$Pe,
                         p2016CA138$Ptheta,
                         p2016CA138$Pomega,
                         nrow=4
                         #,main="164207"
)
Ptrue2016CA138

true_states=order_states(df2016CA138$type)

df2016CA138=df2016CA138[,c("a","e","theta","omega")]

df2016CA138=compute_feat(df2016CA138,wdn=10,a=F,e=F,omega=F)
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

modsel2016CA138=data.frame(hp,FTIC=unlist(lapply(est2016CA138,function(x)x$FTIC)),
                           BAC=unlist(lapply(est2016CA138,function(x)x$BAC))
)
sel=8
df2016CA138=read.table("propagation_2016CA138_new_v2.txt",header=T)
df2016CA138=df2016CA138[,c("a","e","theta","omega")]

df2016CA138=tail(df2016CA138,N)
df2016CA138_b=data.frame(df2016CA138,type=est2016CA138[[sel]]$est_states)
p2016CA138_b=plot_real4(df2016CA138_b)

Pest2016CA138_b=ggarrange(p2016CA138_b$Pa,
                          p2016CA138_b$Pe,
                          p2016CA138_b$Ptheta,
                          p2016CA138_b$Pomega,
                          nrow=4
                          #,main="2015XX169"
)

gridExtra::grid.arrange(Ptrue2016CA138,
                        Pest2016CA138_b)

# 2016CO246 ---------------------------------------------------------------
df2016CO246=read.table("propagation_2016CO246_new_v2.txt",header=T)

true_states=order_states(df2016CO246$type)

p2016CO246=plot_real4(df2016CO246)

Ptrue2016CO246=ggarrange(p2016CO246$Pa,
                         p2016CO246$Pe,
                         p2016CO246$Ptheta,
                         p2016CO246$Pomega,
                         nrow=4
                         #,main="164207"
)
Ptrue2016CO246

df2016CO246=df2016CO246[,c("a","e","theta","omega")]

df2016CO246=compute_feat(df2016CO246,wdn=10,e=F,omega=F)
N=dim(df2016CO246)[1]

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

modsel2016CO246=data.frame(hp,FTIC=unlist(lapply(est2016CO246,function(x)x$FTIC)),
                           BAC=unlist(lapply(est2016CO246,function(x)x$BAC))
)
sel=8
df2016CO246=read.table("propagation_2016CO246_new_v2.txt",header=T)
df2016CO246=df2016CO246[,c("a","e","theta","omega")]

df2016CO246=tail(df2016CO246,N)
df2016CO246_b=data.frame(df2016CO246,type=est2016CO246[[sel]]$est_states)
p2016CO246_b=plot_real4(df2016CO246_b)

Pest2016CO246_b=ggarrange(p2016CO246_b$Pa,
                          p2016CO246_b$Pe,
                          p2016CO246_b$Ptheta,
                          p2016CO246_b$Pomega,
                          nrow=4
                          #,main="2015XX169"
)

gridExtra::grid.arrange(Ptrue2016CO246,
                        Pest2016CO246_b)

# 2014OL339 ---------------------------------------------------------------
df2014OL339=read.table("propagation_2014OL339_new_v2.txt",header=T)
unique(df2014OL339$type)

true_states=order_states(df2014OL339$type)

p2014OL339=plot_real4(df2014OL339)

Ptrue2014OL339=ggarrange(p2014OL339$Pa,
                         p2014OL339$Pe,
                         p2014OL339$Ptheta,
                         p2014OL339$Pomega,
                         nrow=4
                         #,main="164207"
)
Ptrue2014OL339

df2014OL339=df2014OL339[,c("a","e","theta","omega")]

df2014OL339=compute_feat(df2014OL339,am1=T)
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

p2020PN1=plot_real4(df2020PN1)

Ptrue2020PN1=ggarrange(p2020PN1$Pa,
                       p2020PN1$Pe,
                       p2020PN1$Ptheta,
                       p2020PN1$Pomega,
                       nrow=4
                       #,main="164207"
)
Ptrue2020PN1

df2020PN1=df2020PN1[,c("a","e","theta","omega")]

df2020PN1=compute_feat(df2020PN1,e=F,theta=F,wdn=50)
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

modsel2020PN1=data.frame(hp,FTIC=unlist(lapply(est2020PN1,function(x)x$FTIC)),
                           BAC=unlist(lapply(est2020PN1,function(x)x$BAC))
)
sel=7
df2020PN1=read.table("propagation_2020PN1_new_v2.txt",header=T)
df2020PN1=df2020PN1[,c("a","e","theta","omega")]

df2020PN1=tail(df2020PN1,N)
df2020PN1_b=data.frame(df2020PN1,type=est2020PN1[[sel]]$est_states)
p2020PN1_b=plot_real4(df2020PN1_b)

Pest2020PN1_b=ggarrange(p2020PN1_b$Pa,
                          p2020PN1_b$Pe,
                          p2020PN1_b$Ptheta,
                          p2020PN1_b$Pomega,
                          nrow=4
                          #,main="2020PN1"
)
Pest2020PN1_b

gridExtra::grid.arrange(Ptrue2020PN1,Pest2020PN1_b)
