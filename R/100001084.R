df100001084=read.table("100001084.txt",header = T)

# TP, QS and HS

start=which.min(abs(df100001084$t-69770790))
end=which.min(abs(df100001084$t-86694783.10))

df100001084=df100001084[start:end,]

summary(diff(df100001084$t))
#plot(diff(df100001084$t),type='p')

#df100001084$type=-1
data=df100001084
data_name="100001084"

library(ggplot2)
library(dplyr)
library(plotly)

source("Utils4_.R")

# Meglio lag contenuti
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
# 
# Y=data_fin[complete.cases(data_fin),]
# Y=Y[,-1]

data2=merge(data_fin,df100001084,by="t")

data2$tt=1:dim(data2)[1]
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



########
# Adjust TP, say set TP=0 whenever it comes after a QS

# temp=data2$I_TD+2*data2$I_HS+3*data2$I_QS
# for(i in (tail(which(is.na(temp)),1)+2):length(temp)){
#   if(temp[i]==1&temp[i-1]==3){
#     temp[i]=0
#   }
#   else if(temp[i]==1&temp[i-1]==0){
#     temp[i]=0
#   }
# }
# 
# data2$I_TD=as.numeric(temp==1)
########

p_I <- ggplot(data2) + 
  labs(title = "100001084", x = "tt", y = "Theta") +
  geom_line(aes(x=tt,y=theta),col='grey')+
  geom_point(aes(x=tt,y=theta),col='grey40',size=1)+
  geom_point(aes(x = tt, y = theta, color = factor(maxs_mins)), size = 1) +
  scale_color_manual(values = c("0" = "grey50", "1"='blue' ,"2" = "red")) +
  geom_line(aes(x=tt,y=I_TD),col='green4')+
  geom_line(aes(x=tt,y=I_HS),col='cyan3')+
  geom_line(aes(x=tt,y=I_QS),col="violet")+
  geom_line(aes(x=tt,y=mean_osc),col='black')+
  theme_minimal() +
  theme(legend.position = "none") 
#p_I
# Convert the ggplot object to a plotly interactive plot
ggplotly(p_I)


data2=data.frame(data2)
Y=data2[complete.cases(data2),]
Y=Y[,colnames(data_fin)[-1]]


# fit ---------------------------------------------------------------------

# Tolgo I_TD perche a volte succede che ho un massimo in HS e un minimo in QS che hanno stesso segno e vengono intepretati come TP
# Y=subset(Y,select=-I_TD) 
# Purtroppo cosi seleziona solo due regimi

Y=subset(Y,select=-c(I_diffmaxmin,diffmaxmin,maxs,mins,value_max,value_min)) 

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
start_100001084=Sys.time()
est100001084 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100001084=Sys.time()
elapsed_100001084=end_100001084-start_100001084
save(data_fin,Y,est100001084,file="100001084.RData")

modsel100001084=data.frame(hp,
                           FTIC=unlist(lapply(est100001084,function(x)x$FTIC))
)

best_mod=modsel100001084[which.min(modsel100001084$FTIC),]
best_mod

sel=57
#sel=68
estw100001084=data.frame(var=colnames(Y),
                         weight=est100001084[[sel]]$est_weights)

estw100001084=estw100001084[order(estw100001084$weight,decreasing = T),]
head(estw100001084,15)

# Merge with a theta and t for plotting
data_a_theta=tail(data[,c("t","a","theta")],dim(Y)[1])

df_res_100001084=data.frame(data_a_theta,
                            Y,
                            State=est100001084[[sel]]$est_states
)

df_res_100001084 <- df_res_100001084 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100001084 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a))

p_a_res <- ggplot(data = df_segments_a) + 
  geom_segment(aes(x = t, y = a, 
                   xend = next_t, yend = next_a), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=a,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(df_res_100001084$State)) +
  labs(title = "100001084", 
       x = "Time (t)", 
       y = "Values (a)") +
  theme_minimal()
ggplotly(p_a_res)

df_segments_theta <- df_res_100001084 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))

p_theta_res <- ggplot(data = df_segments_theta) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=theta,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(df_res_100001084$State)) +
  labs(title = "100001084", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)

# Check feat by feat
pp <- ggplot(data=df_res_100001084,aes(x=t)) + 
  geom_point(aes(y = sd_dtheta_long, 
                 color = as.factor(State)), 
             size = 1) +
  scale_color_manual(values = 1:max(df_res_100001084$State)) +
  labs(title = "100001084", 
       x = "Time (t)", 
       y = " ") +
  theme_minimal()

pp
