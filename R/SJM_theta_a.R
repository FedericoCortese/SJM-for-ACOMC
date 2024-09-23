source("Utils3_.R")
library(reticulate)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(zoo)
import("scipy")
source_python('SJ.py')

# 164207 ------------------------------------------------------------------

df164207_true=read.table("propagation_164207_new_v2.txt",header=T)
wdn=c(10,50)
df164207=comp_feat_theta_a(df164207_true,wdn)
sat_mod=SJM_sat(df164207)
Lnsat=sat_mod$Lnsat
lambda=.1
kappa=1
prova=SJM_lambdakappa(lambda,kappa,K=2,df164207,Lnsat,Ksat=6,
                      alpha0=NULL,K0=NULL,pers0=.95,
                      true_states=tail(df164207_true$type,dim(df164207)[1]))

plot(df164207$t,df164207$theta,type="p",col=prova$est_states)

# 2016CA138 ------------------------------------------------------------------
df2016CA138_true=read.table("propagation_2016CA138_new_v2.txt",header=T)
wdn=c(10,50)
df2016CA138=comp_feat_theta_a(df2016CA138_true,wdn)
sat_mod=SJM_sat(df2016CA138)
Lnsat=sat_mod$Lnsat
lambda=.2
kappa=2
prova=SJM_lambdakappa(lambda,kappa,K=4,df2016CA138,Lnsat,Ksat=6,
                      alpha0=NULL,K0=NULL,pers0=.95,
                      true_states=tail(df2016CA138_true$type,dim(df2016CA138)[1]))
plot(df2016CA138$t,df2016CA138$a,type="p",col=prova$est_states)


# 2015XX169 ---------------------------------------------------------------

df2015XX169_true=read.table("propagation_2015XX169_new_v2.txt",header=T)
wdn=c(10,50)
df2015XX169=comp_feat_theta_a(df2015XX169_true,wdn)
sat_mod=SJM_sat(df2015XX169)
Lnsat=sat_mod$Lnsat
lambda=0
kappa=2
prova=SJM_lambdakappa(lambda,kappa,K=3,df2015XX169,Lnsat,Ksat=6,
                      alpha0=NULL,K0=NULL,pers0=.95,
                      true_states=tail(df2015XX169_true$type,dim(df2015XX169)[1]))
plot(df2015XX169$t,df2015XX169$theta,type="p",col=prova$est_states)
