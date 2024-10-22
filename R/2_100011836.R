source("Utils3_.R")
library(plotly)

df100011836=read.table("100011836.txt",header = T)

start=which.min(abs(df100011836$t-7.9e8))
end=which.min(abs(df100011836$t-8.4e8))

df100011836=df100011836[start:end,]

summary(diff(df100011836$t))
plot(diff(df100011836$t),type='p')

#df100011836$type=-1
data=df100011836
data_name="100011836"

names(df100011836)

df100011836=df100011836[,c("t","a","theta")]

WW=c(10,30)
tt_100011836=feat_comput_theta_a(df100011836, "100011836",
                                 l=5,
                                 wdn=WW,
                              l2=WW,
                              l3=WW,
                              tt_thres_maxmin=2.7,
                              tt_thres_diffmaxmin=0.25)

Y100011836_final=tt_100011836[complete.cases(tt_100011836),]
time_100011836=Y100011836_final$t
Y100011836_final=subset(Y100011836_final,select=-t)

lambda=c(0,5,10,15,20,30)
K=3:5
kappa=seq(1,ceiling(sqrt(dim(Y100011836_final)[2])),by=1.5)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y100011836_final)
Lnsat=sat_mod$Lnsat
start_100011836=Sys.time()
est100011836 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y100011836_final,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100011836=Sys.time()
elapsed_100011836=end_100011836-start_100011836
save(time_100011836,Y100011836_final,est100011836,file="100011836.RData")

modsel100011836=data.frame(hp,
                           FTIC=unlist(lapply(est100011836,function(x)x$FTIC))
)
best_mod=modsel100011836[which.min(modsel100011836$FTIC),]
best_mod