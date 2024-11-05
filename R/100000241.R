df100000241=read.table("100000241.txt",header = T)

start=which.min(abs(df100000241$t-83183570))
end=which.min(abs(df100000241$t-244373585))

df100000241=df100000241[start:end,]

summary(diff(df100000241$t))
which(diff(df100000241$t)!=1000)
#plot(diff(df100000241$t),type='p')


df100000241=df100000241[-(1:4850),]
summary(diff(df100000241$t))


#df100000241$type=-1
data=df100000241
data_name="100000241"

library(ggplot2)
library(dplyr)
library(plotly)

source("Utils4_.R")

max_lag=10
min_lag=5
# %FC importante che questa soglia sia abbastanza grande da distinguere QS dal resto
tt_thres_diffmaxmin=0.07

data_thetavol=thetavol_feat(data,
                            wdn=c(min_lag,max_lag))
data_a=a_feat(data,l3=c(min_lag,max_lag))

data_maxmin=max_min_feat(data,
                         tt_thres_maxmin = 2.4,
                         l2=c(min_lag,max_lag),
                         tt_thres_diffmaxmin=tt_thres_diffmaxmin)
names(data_maxmin)

#data_maxmin=data_maxmin[,c("t","I_TD","I_HS","I_QS","mean_osc")]

# data_maxmin=data_maxmin[,c("t","mean_osc")]
# %FC provo a reintrodurre value max min, e diff max min
# data_maxmin=data_maxmin[,c("t","mean_osc","value_min" ,"value_max","diffmaxmin","I_diffmaxmin")]

data_fin=merge(data_thetavol,data_a,by="t")
data_fin=merge(data_fin,data_maxmin,by="t")

names(data_fin)

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]

data2=merge(data_fin,df100000241,by="t")
data2 <- data2 %>%
  mutate(Segment_max = cumsum(maxs != lag(maxs, default = first(maxs))),
         Segment_min = cumsum(mins != lag(mins, default = first(mins))),
         maxs_mins=mins+maxs*2)

data2 <- data2 %>%
  group_by(Segment_max) %>%
  mutate(next_t_max = dplyr::lead(t), next_theta_max = dplyr::lead(theta))

data2 <- data2 %>%
  group_by(Segment_min) %>%
  mutate(next_t_min = dplyr::lead(t), next_theta_min = dplyr::lead(theta))


p_I <- ggplot(data2) + 
  labs(title = "100000241", x = "t", y = "Theta") +
  geom_line(aes(x=t,y=theta),col='grey')+
  geom_point(aes(x=t,y=theta),col='grey40',size=1)+
  geom_point(aes(x = t, y = theta, color = factor(maxs_mins)), size = 1) +
  scale_color_manual(values = c("0" = "grey50", "1"='blue' ,"2" = "red")) +
  geom_line(aes(x=t,y=I_TD),col='green4')+
  geom_line(aes(x=t,y=I_HS),col='cyan3')+
  geom_line(aes(x=t,y=I_QS),col="violet")+
  #geom_line(aes(x=t,y=I_diffmaxmin*2),col='black')+
  theme_minimal() +
  theme(legend.position = "none") 
#p_I
# Convert the ggplot object to a plotly interactive plot
ggplotly(p_I)

# fit ---------------------------------------------------------------------

lambda=c(0,5,10,15,20,30)
K=2:7
# %FC Con K da 2 a 6 seleziona K=3
# K=5
# %FC Provato a fissare K=5 ma non vede compound
# K=4

kappa=seq(1,ceiling(sqrt(dim(Y)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y)
Lnsat=sat_mod$Lnsat
start_100000241=Sys.time()
est100000241 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100000241=Sys.time()
elapsed_100000241=end_100000241-start_100000241
save(data_fin,Y,est100000241,file="100000241.RData")

modsel100000241=data.frame(hp,
                           FTIC=unlist(lapply(est100000241,function(x)x$FTIC))
)

best_mod=modsel100000241[which.min(modsel100000241$FTIC),]
best_mod
