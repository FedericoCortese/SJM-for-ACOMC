source("Utils4_.R")
# 2002AA29 ----------------------------------------------------------------

df2002AA29=read.table("propagation_2002AA29_new_v2.txt",header=T)

summary(diff(df2002AA29$t))
# Step 1000

# 2016HO3 -----------------------------------------------------------------

df2016HO3=read.table("propagation_2016HO3_new_v2.txt",header=T)

summary(diff(df2016HO3$t))

# Step 1000
min_lag=10
max_lag=75

g_truth=df2016HO3[,c("t","type")]

df2016HO3=subset(df2016HO3,select=-c(type,i,u))

data_old_feat=compute_feat(df2016HO3,
                           wdn=c(min_lag,max_lag),
                           Omega=F)

data_a=a_feat(df2016HO3,l3=c(min_lag,max_lag))

data_maxmin=max_min_feat(df2016HO3,
                         l=5,
                         tt_thres_maxmin = 2.4,
                         l2=c(min_lag,max_lag),
                         tt_thres_diffmaxmin=.7)

names(data_old_feat)

# Y=merge(data_old_feat,data_a,by='t')
# Y=merge(Y,data_maxmin,by='t')
# Y=Y[complete.cases(Y),]
# head(Y)
Y=data_old_feat
Y=Y[complete.cases(Y),]

#Y=subset(Y,select=-c(I_diffmaxmin,diffmaxmin,maxs,mins))

timesY=Y[,1]
Y=Y[,-1]

lambda=c(0,5,10,15,20,30)
K=2:6
kappa=seq(1,ceiling(sqrt(dim(Y)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y)
Lnsat=sat_mod$Lnsat
start_2016HO3=Sys.time()
est2016HO3 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_2016HO3=Sys.time()
elapsed_2016HO3=end_2016HO3-start_2016HO3
save(timesY,Y,est2016HO3,file="3_2016HO3.RData")

modsel2016HO3=data.frame(hp,
                           FTIC=unlist(lapply(est2016HO3,function(x)x$FTIC))
)

best_mod=modsel2016HO3[which.min(modsel2016HO3$FTIC),]
best_mod

sel=190
#sel=68
estw2016HO3=data.frame(var=colnames(Y),
                         weight=est2016HO3[[sel]]$est_weights)

estw2016HO3=estw2016HO3[order(estw2016HO3$weight,decreasing = T),]
head(estw2016HO3,15)

# 1000006174 --------------------------------------------------------------

df100006174=read.table("100006174.txt",header = T)

start=which.min(abs(df100006174$t-1.6e8))
end=which.min(abs(df100006174$t-1.97e8))

df100006174=df100006174[start:end,]
summary(diff(df100006174$t))
# Step 7500

# 100011836 ---------------------------------------------------------------





