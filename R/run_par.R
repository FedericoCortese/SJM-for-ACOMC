source("Utils3_.R")
library(reticulate)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(zoo)
import("scipy")
source_python('SJ.py')

# 164207 ------------------------------------------------------------------

df164207=read.table("propagation_164207_new_v2.txt",header=T)
p164207=plot_real(df164207,"blue")
Ptrue164207=ggarrange(p164207$Pa,
                      p164207$Pe,
                      p164207$Ptheta,
                      p164207$Pomega,
                      nrow=4
                      #,main="164207"
)
Ptrue164207

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
#   library(dplyr)})
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
save(df164207,est164207,elapsed_est164207,file="est164207.RData")

modsel164207=data.frame(hp,FTIC=unlist(lapply(est164207,function(x)x$FTIC)),
                        BAC=unlist(lapply(est164207,function(x)x$BAC))
                        #ARI=unlist(lapply(est164207,function(x)x$ARI))
)
sel=25
modsel164207

estw164207=data.frame(weight=est164207[[sel]]$est_weights,feat=colnames(df164207))

# Sort by weight
estw164207=estw164207[order(estw164207$weight,decreasing = T),]
estw164207

df164207_b=data.frame(df164207,type=est164207[[sel]]$est_states)
p164207_b=plot_real(df164207_b,"red")

Pest164207_b=ggarrange(p164207_b$Pa,
                       p164207_b$Pe,
                       p164207_b$Ptheta,
                       p164207_b$Pomega,
                       nrow=4
                       #,main="164207"
)

png("164207EST.png",width =1600,height=900)
annotate_figure(Pest164207_b,
                top = text_grob("164207 - est - 96% acc.", color = "black", face = "bold", size = 14))
dev.off()

# 2001GO2 -------------------------------------------------------------------------

df2001GO2=read.table("propagation_2001GO2_new_v2.txt",header=T)
true_states=order_states(df2001GO2$type)

p2001GO2=plot_real(df2001GO2,"blue")

Ptrue2001GO2=ggarrange(p2001GO2$Pa,
                        p2001GO2$Pe,
                        p2001GO2$Ptheta,
                        p2001GO2$Pomega,
                        nrow=4
                        #,main="2002AA29"
)
Ptrue2001GO2
#true_states=df2001GO2$type

df2001GO2=df2001GO2[,c("a","e","theta","omega")]

df2001GO2=compute_feat(df2001GO2,wdn=10,a=F)
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
save(df2001GO2,est2001GO2,elapsed_est2001GO2,file="est2001GO2.RData")

modsel2001GO2=data.frame(hp,FTIC=unlist(lapply(est2001GO2,function(x)x$FTIC)),
                         BAC=unlist(lapply(est2001GO2,function(x)x$BAC))
                         #ARI=unlist(lapply(est2001GO2,function(x)x$ARI))
)

sel=14
df2001GO2=read.table("propagation_2001GO2_new_v2.txt",header=T)
df2001GO2=df2001GO2[,c("a","e","theta","omega")]
df2001GO2=tail(df2001GO2,N)
df2001GO2_b=data.frame(df2001GO2,type=est2001GO2[[sel]]$est_states)
p2001GO2_b=plot_real(df2001GO2_b,"red")
Pest2001GO2_b=ggarrange(p2001GO2_b$Pa,
                        p2001GO2_b$Pe,
                        p2001GO2_b$Ptheta,
                        p2001GO2_b$Pomega,
                        nrow=4
                        #,main="2001GO2"
)

png("2001GO2EST.png",width =1600,height=900)
annotate_figure(Pest2001GO2_b,
                top = text_grob("2001GO2 - est - 90% acc.", color = "black", face = "bold", size = 14))
dev.off()


# 2002AA29 --------------------------------------------------------------------

df2002AA29=read.table("propagation_2002AA29_new_v2.txt",header=T)
true_states=order_states(df2002AA29$type)

p2002AA29=plot_real(df2002AA29,"blue")

Ptrue2002AA29=ggarrange(p2002AA29$Pa,
                        p2002AA29$Pe,
                        p2002AA29$Ptheta,
                        p2002AA29$Pomega,
                        nrow=4
                        #,main="2002AA29"
)
Ptrue2002AA29
true_states=df2002AA29$type

df2002AA29=df2002AA29[,c("a","e","theta","omega")]

df2002AA29=compute_feat(df2002AA29,wdn=10)
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
save(df2002AA29,est2002AA29,elapsed_est2002AA29,file="est2002AA29.RData")


modsel2002AA29=data.frame(hp,FTIC=unlist(lapply(est2002AA29,function(x)x$FTIC)),
                         BAC=unlist(lapply(est2002AA29,function(x)x$BAC))
                         #ARI=unlist(lapply(est2001GO2,function(x)x$ARI))
)

sel=13
df2002AA29_b=data.frame(df2002AA29,type=est2002AA29[[sel]]$est_states)
p2002AA29_b=plot_real(df2002AA29_b,"red")
Pest2002AA29_b=ggarrange(p2002AA29_b$Pa,
                        p2002AA29_b$Pe,
                        p2002AA29_b$Ptheta,
                        p2002AA29_b$Pomega,
                        nrow=4
                        #,main="2001GO2"
)

png("2002AA29EST.png",width =1600,height=900)
annotate_figure(Pest2002AA29_b,
                top = text_grob("2002AA29 - est - 98% acc.", color = "black", face = "bold", size = 14))
dev.off()

# 2015SO2 --------------------------------------------------------------------

df2015SO2=read.table("propagation_2015SO2_new_v2.txt",header=T)
true_states=order_states(df2015SO2$type)

p2015SO2=plot_real(df2015SO2,"blue")

Ptrue2015SO2=ggarrange(p2015SO2$Pa,
                       p2015SO2$Pe,
                        p2015SO2$Ptheta,
                        p2015SO2$Pomega,
                        nrow=4
                        #,main="2002AA29"
)
Ptrue2015SO2
true_states=df2015SO2$type

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
save(df2015SO2,est2015SO2,elapsed_est2015SO2,file="est2015SO2.RData")

modsel2015SO2=data.frame(hp,FTIC=unlist(lapply(est2015SO2,function(x)x$FTIC)),
                         BAC=unlist(lapply(est2015SO2,function(x)x$BAC))
                         #ARI=unlist(lapply(est2001GO2,function(x)x$ARI))
)

sel=25
df2015SO2_b=data.frame(df2015SO2,type=est2015SO2[[sel]]$est_states)
p2015SO2_b=plot_real(df2015SO2_b,"red")
Pest2015SO2_b=ggarrange(p2015SO2_b$Pa,
                        p2015SO2_b$Pe,
                        p2015SO2_b$Ptheta,
                        p2015SO2_b$Pomega,
                        nrow=4
                        #,main="2001GO2"
)

png("2015SO2EST.png",width =1600,height=900)
annotate_figure(Pest2015SO2_b,
                top = text_grob("2015SO2 - est - 85% acc.", color = "black", face = "bold", size = 14))
dev.off()

# 2016HO3 -----------------------------------------------------------------

df2016HO3=read.table("propagation_2016HO3_new_v2.txt",header=T)
true_states=order_states(df2016HO3$type)

p2016HO3=plot_real(df2016HO3,"blue")

Ptrue2016HO3=ggarrange(p2016HO3$Pa,
                       p2016HO3$Pe,
                       p2016HO3$Ptheta,
                       p2016HO3$Pomega,
                       nrow=4
                       #,main="2016HO3"
)
Ptrue2016HO3

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
save(df2016HO3,est2016HO3,elapsed_est2016HO3,file="est2016HO3.RData")

modseldf2016HO3=data.frame(hp,FTIC=unlist(lapply(est2016HO3,function(x)x$FTIC)),
                          BAC=unlist(lapply(est2016HO3,function(x)x$BAC))
                          #ARI=unlist(lapply(est2001GO2,function(x)x$ARI))
)

sel=25
df2016HO3_b=data.frame(df2016HO3,type=est2016HO3[[sel]]$est_states)
p2016HO3_b=plot_real(df2016HO3_b,"red")
Pest2016HO3_b=ggarrange(p2016HO3_b$Pa,
                         p2016HO3_b$Pe,
                         p2016HO3_b$Ptheta,
                         p2016HO3_b$Pomega,
                         nrow=4
                         #,main="2001GO2"
)

png("2016HO3EST.png",width =1600,height=900)
annotate_figure(Pest2016HO3_b,
                top = text_grob("2016HO3 - est - 94% acc.", color = "black", face = "bold", size = 14))
dev.off()

# 2019GM1 -----------------------------------------------------------------

df2019GM1=read.table("propagation_2019GM1_new_v2.txt",header=T)
true_states=order_states(df2019GM1$type)

p2019GM1=plot_real(df2019GM1,"blue")

Ptrue2019GM1=ggarrange(p2019GM1$Pa,
                       p2019GM1$Pe,
                       p2019GM1$Ptheta,
                       p2019GM1$Pomega,
                       nrow=4
                       #,main="2019GM1"
)
Ptrue2019GM1

df2019GM1=df2019GM1[,c("a","e","theta","omega")]

df2019GM1=compute_feat(df2019GM1,wdn=10,a=F,e=F)
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
save(df2019GM1,est2019GM1,elapsed_est2019GM1,file="est2019GM1.RData")

modsel2019GM1=data.frame(hp,FTIC=unlist(lapply(est2019GM1,function(x)x$FTIC)),
                           BAC=unlist(lapply(est2019GM1,function(x)x$BAC))
                           #ARI=unlist(lapply(est2001GO2,function(x)x$ARI))
)

sel=7

df2019GM1=read.table("propagation_2019GM1_new_v2.txt",header=T)
df2019GM1=df2019GM1[,c("a","e","theta","omega")]

df2019GM1=tail(df2019GM1,N)
df2019GM1_b=data.frame(df2019GM1,type=est2019GM1[[sel]]$est_states)
p2019GM1_b=plot_real(df2019GM1_b,"red")
Pest2019GM1_b=ggarrange(p2019GM1_b$Pa,
                        p2019GM1_b$Pe,
                        p2019GM1_b$Ptheta,
                        p2019GM1_b$Pomega,
                        nrow=4
                        #,main="2001GO2"
)

png("2019GM1EST.png",width =1600,height=900)
annotate_figure(Pest2019GM1_b,
                top = text_grob("2019GM1 - est - 95%", color = "black", face = "bold", size = 14))
dev.off()

# 2020PP1 -----------------------------------------------------------------

df2020PP1=read.table("propagation_2020PP1_new_v2.txt",header=T)
true_states=order_states(df2020PP1$type)

p2020PP1=plot_real(df2020PP1,"blue")

Ptrue2020PP1=ggarrange(p2020PP1$Pa,
                       p2020PP1$Pe,
                       p2020PP1$Ptheta,
                       p2020PP1$Pomega,
                       nrow=4
                       #,main="2020PP1"
)
Ptrue2020PP1

df2020PP1=df2020PP1[,c("a","e","theta","omega")]

df2020PP1=compute_feat(df2020PP1,wdn=10,a=F)
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
#   library(dplyr)})
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
save(df2020PP1,est2020PP1,elapsed_est2020PP1,file="est2020PP1.RData")

modsel2020PP1=data.frame(hp,FTIC=unlist(lapply(est2020PP1,function(x)x$FTIC)),
                           BAC=unlist(lapply(est2020PP1,function(x)x$BAC))
                           #ARI=unlist(lapply(est2001GO2,function(x)x$ARI))
)

sel=7
df2020PP1=read.table("propagation_2020PP1_new_v2.txt",header=T)
df2020PP1=df2020PP1[,c("a","e","theta","omega")]

df2020PP1=tail(df2020PP1,N)
df2020PP1_b=data.frame(df2020PP1,type=est2020PP1[[sel]]$est_states)
p2020PP1_b=plot_real(df2020PP1_b,"red")
Pest2020PP1_b=ggarrange(p2020PP1_b$Pa,
                        p2020PP1_b$Pe,
                        p2020PP1_b$Ptheta,
                        p2020PP1_b$Pomega,
                        nrow=4
                        #,main="2001GO2"
)

png("2020PP1EST.png",width =1600,height=900)
annotate_figure(Pest2020PP1_b,
                top = text_grob("2020PP1 - est - 95% acc.", color = "black", face = "bold", size = 14))
dev.off()

Pest2020PP1_b
