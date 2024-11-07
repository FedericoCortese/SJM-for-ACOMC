#source("Utils3_.R")
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

names(Y)

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
start_100011836=Sys.time()
est100011836 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100011836=Sys.time()
elapsed_100011836=end_100011836-start_100011836
save(data_fin,Y,est100011836,file="3_100011836.RData")

modsel100011836=data.frame(hp,
                           FTIC=unlist(lapply(est100011836,function(x)x$FTIC))
)

best_mod=modsel100011836[which.min(modsel100011836$FTIC),]
best_mod

sel=69
#sel=68
estw100011836=data.frame(var=colnames(Y),
                         weight=est100011836[[sel]]$est_weights)

estw100011836=estw100011836[order(estw100011836$weight,decreasing = T),]
head(estw100011836,15)

# Merge with a theta and t for plotting
data_a_theta=tail(data[,c("t","a","theta")],dim(Y)[1])

df_res_100011836=data.frame(data_a_theta,
                            Y,
                            State=est100011836[[sel]]$est_states
)

df_res_100011836 <- df_res_100011836 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100011836 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a))

p_a_res <- ggplot(data = df_segments_a) + 
  geom_segment(aes(x = t, y = a, 
                   xend = next_t, yend = next_a), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=a,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = "Values (a)") +
  theme_minimal()
ggplotly(p_a_res)

df_segments_theta <- df_res_100011836 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))

p_theta_res <- ggplot(data = df_segments_theta) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=theta,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)
