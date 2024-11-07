#source("Utils3_.R")
library(plotly)

df100006174=read.table("100006174.txt",header = T)

start=which.min(abs(df100006174$t-1.6e8))
end=which.min(abs(df100006174$t-1.97e8))

df100006174=df100006174[start:end,]

summary(diff(df100006174$t))
plot(diff(df100006174$t),type='p')

#df100006174$type=-1
data=df100006174
data_name="100006174"

names(df100006174)

source("Utils4_.R")

#max_lag=30
max_lag=10
min_lag=5
# %FC importante che questa soglia sia abbastanza grande da distinguere QS dal resto
tt_thres_diffmaxmin=pi/4

data_thetavol=thetavol_feat(data,
                            wdn=c(min_lag,max_lag))
data_a=a_feat(data,l3=c(min_lag,max_lag))

data_maxmin=max_min_feat(data,
                         l=3,
                         tt_thres_maxmin = 2.4,
                         l2=c(min_lag,max_lag),
                         tt_thres_diffmaxmin=tt_thres_diffmaxmin)
names(data_maxmin)

data_maxmin_full=data_maxmin
data_maxmin=data_maxmin[,c("t","I_TD","I_HS","I_QS","mean_osc")]

# data_maxmin=data_maxmin[,c("t","mean_osc")]
# %FC provo a reintrodurre value max min, e diff max min
# data_maxmin=data_maxmin[,c("t","mean_osc","value_min" ,"value_max","diffmaxmin","I_diffmaxmin")]

data_fin=merge(data_thetavol,data_a,by="t")
data_fin=merge(data_fin,data_maxmin,by="t")
data_fin=merge(data_fin,data[,c("t","a","theta")],by="t")

names(data_fin)

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]

lambda=c(0,5,10,15,20,30)
K=2:6
# %FC Con K da 2 a 6 seleziona K=3
# K=5
# %FC Provato a fissare K=5 ma non vede compound
# K=4

kappa=seq(1,ceiling(sqrt(dim(Y)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y)
Lnsat=sat_mod$Lnsat
start_100006174=Sys.time()
est100006174 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100006174=Sys.time()
elapsed_100006174=end_100006174-start_100006174
save(data_fin,Y,est100006174,file="3_100006174.RData")

modsel100006174=data.frame(hp,
                           FTIC=unlist(lapply(est100006174,function(x)x$FTIC))
)

best_mod=modsel100006174[which.min(modsel100006174$FTIC),]
best_mod

sel=47
#sel=68
estw100006174=data.frame(var=colnames(Y),
                         weight=est100006174[[sel]]$est_weights)

estw100006174=estw100006174[order(estw100006174$weight,decreasing = T),]
head(estw100006174,15)

# Merge with a theta and t for plotting
data_a_theta=tail(data[,c("t","a","theta")],dim(Y)[1])

df_res_100006174=data.frame(data_a_theta,
                            Y,
                            State=est100006174[[sel]]$est_states
)

df_res_100006174 <- df_res_100006174 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100006174 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a))

p_a_res <- ggplot(data = df_segments_a) + 
  geom_segment(aes(x = t, y = a, 
                   xend = next_t, yend = next_a), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=a,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(df_res_100006174$State)) +
  labs(title = "100006174", 
       x = "Time (t)", 
       y = "Values (a)") +
  theme_minimal()
ggplotly(p_a_res)

df_segments_theta <- df_res_100006174 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))

p_theta_res <- ggplot(data = df_segments_theta) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=theta,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(df_res_100006174$State)) +
  labs(title = "100006174", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)
