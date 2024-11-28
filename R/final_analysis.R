source("Utils5_.R")
# 2002AA29 ----------------------------------------------------------------

df2002AA29=read.table("propagation_2002AA29_new_v2.txt",header=T)
#df2002AA29=trans_theta(df2002AA29)

ground_truth=df2002AA29$type
df2002AA29=subset(df2002AA29,select=-c(type,u,i))
names(df2002AA29)

summary(diff(df2002AA29$t))
# Step 1000

summary(df2002AA29$theta)
plot(df2002AA29$theta,type='l')

Y=compute_feat(df2002AA29,wdn=10,wdn_decomp=10,
               a=T,e=T,theta=T,omega=T)
names(Y)
head(Y)

timesY=Y[,1]
Y=Y[,-1]
Y=Y[complete.cases(Y),]
N=dim(Y)[1]
ground_truth=tail(ground_truth,N)

lambda=c(0,5,10,15,20,30)
K=2
kappa=seq(1,ceiling(sqrt(dim(Y)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y)
Lnsat=sat_mod$Lnsat
start_2002AA29=Sys.time()
est2002AA29 <- parallel::mclapply(1:nrow(hp),
                                  function(x)
                                    SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    df=Y,
                                                    Lnsat=Lnsat),
                                  mc.cores = parallel::detectCores())

end_2002AA29=Sys.time()
elapsed_2002AA29=end_2002AA29-start_2002AA29
save(timesY,Y,est2002AA29,file="3_2002AA29.RData")

modsel2002AA29=data.frame(hp,
                          FTIC=unlist(lapply(est2002AA29,function(x)x$FTIC))
)
best_mod=modsel2002AA29[which.min(modsel2002AA29$FTIC),]
best_mod

sel=13
estw2002AA29=data.frame(var=colnames(Y),
                         weight=est2002AA29[[sel]]$est_weights)

estw2002AA29=estw2002AA29[order(estw2002AA29$weight,decreasing = T),]
head(estw2002AA29,15)

est_states2002AA29=est2002AA29[[sel]]$est_states
ground_truth=order_states(ground_truth)
est_states2002AA29=order_states(est_states2002AA29)

BAC(est_states2002AA29,ground_truth)

# mclust::adjustedRandIndex(ground_truth,est_states2002AA29)

dfresINV_2002AA29=data.frame(t=tail(timesY,N),
                          Y,
                          State=est_states2002AA29)

dfres_2002AA29=data.frame(t=tail(timesY,N),
                          trans_theta(Y),
                          State=est_states2002AA29)


dfres_2002AA29 <-  dfres_2002AA29 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

df_segments_theta_2002AA29 <- dfres_2002AA29 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))


label_size=18
p_theta_res_2002AA29 <- ggplot(data = df_segments_theta_2002AA29) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=theta,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(dfres_2002AA29$State),
                     labels = c("HS", "QS"),
                     name="Orbital regime") +
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(theta ~ "(rad)")) +
  scale_y_continuous(
    breaks = c(-pi, 0, pi),  # Specify where to place the labels
    labels = c(expression(-pi), expression(0), expression(pi))  # Use LaTeX-style labels
  )+
  theme_minimal()+
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    #plot.title = element_text(size = 18),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "right"                  
  )+
  guides(color = guide_legend(override.aes = list(size = 5)))

p_theta_res_2002AA29

# State conditional summ stat
tapply(dfres_2002AA29$theta,dfres_2002AA29$State,mean)
tapply(dfres_2002AA29$theta,dfres_2002AA29$State,sd)

tapply(c(NA,diff(dfres_2002AA29$theta)),dfres_2002AA29$State,mean,na.rm=T)
tapply(c(NA,diff(dfres_2002AA29$theta)),dfres_2002AA29$State,sd,na.rm=T)

tapply(dfres_2002AA29$omega,dfres_2002AA29$State,mean)
tapply(dfres_2002AA29$omega,dfres_2002AA29$State,sd)

tapply(dfres_2002AA29$e,dfres_2002AA29$State,mean)
tapply(dfres_2002AA29$e,dfres_2002AA29$State,sd)

table(est_states2002AA29)/length(est_states2002AA29)*100
# Identify transitions
sequence=est_states2002AA29
transitions <- table(sequence[-length(sequence)], sequence[-1])
transition_probs <- prop.table(transitions, 1)
print(transition_probs)
run_lengths <- rle(sequence)
avg_duration <- tapply(run_lengths$lengths, run_lengths$values, mean)
print(avg_duration)


# 2016HO3 -----------------------------------------------------------------

df2016HO3=read.table("propagation_2016HO3_new_v2.txt",header=T)
#df2016HO3=trans_theta(df2016HO3)

ground_truth=df2016HO3$type
df2016HO3=subset(df2016HO3,select=-c(type,u,i))
names(df2016HO3)

summary(diff(df2016HO3$t))
# Step 1000

summary(df2016HO3$theta)
plot(df2016HO3$theta,type='l')

Y=compute_feat(df2016HO3,wdn=10,wdn_decomp=10,
               a=T,e=T,theta=T,omega=T)
names(Y)
head(Y)

timesY=Y[,1]
Y=Y[,-1]
Y=Y[complete.cases(Y),]
N=dim(Y)[1]
ground_truth=tail(ground_truth,N)

lambda=c(0,5,10,15,20,30)
K=2
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
                                 mc.cores = parallel::detectCores())

end_2016HO3=Sys.time()
elapsed_2016HO3=end_2016HO3-start_2016HO3
save(timesY,Y,est2016HO3,file="3_2016HO3.RData")

modsel2016HO3=data.frame(hp,
                         FTIC=unlist(lapply(est2016HO3,function(x)x$FTIC))
)
best_mod=modsel2016HO3[which.min(modsel2016HO3$FTIC),]
best_mod

sel=15
estw2016HO3=data.frame(var=colnames(Y),
                       weight=est2016HO3[[sel]]$est_weights)

estw2016HO3=estw2016HO3[order(estw2016HO3$weight,decreasing = T),]
head(estw2016HO3,15)

est_states2016HO3=est2016HO3[[sel]]$est_states
ground_truth=order_states(ground_truth)
est_states2016HO3=order_states(est_states2016HO3)

BAC(est_states2016HO3,ground_truth)

dfresINV_2016HO3=data.frame(t=tail(timesY,N),
                            Y,
                            State=est_states2016HO3)

dfres_2016HO3=data.frame(t=tail(timesY,N),
                         trans_theta(Y),
                         State=est_states2016HO3)


dfres_2016HO3 <-  dfres_2016HO3 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

df_segments_theta_2016HO3 <- dfres_2016HO3 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))


label_size=18
p_theta_res_2016HO3 <- ggplot(data = df_segments_theta_2016HO3) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=theta,
                 color=as.factor(State)))+
  # scale_color_manual(values = 1:max(dfres_2016HO3$State),
  #                    labels = c("HS", "QS"),
  #                    name="Orbital regime") +
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(theta ~ "(rad)")) +
  scale_y_continuous(
    breaks = c(-pi, 0, pi),  # Specify where to place the labels
    labels = c(expression(-pi), expression(0), expression(pi))  # Use LaTeX-style labels
  )+
  theme_minimal()+
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "right"                  
  )+
  guides(color = guide_legend(override.aes = list(size = 5)))

p_theta_res_2016HO3

# State conditional summ stat
tapply(dfres_2016HO3$theta,dfres_2016HO3$State,mean)
tapply(dfres_2016HO3$theta,dfres_2016HO3$State,sd)

tapply(c(NA,diff(dfres_2016HO3$theta)),dfres_2016HO3$State,mean,na.rm=T)
tapply(c(NA,diff(dfres_2016HO3$theta)),dfres_2016HO3$State,sd,na.rm=T)

tapply(dfres_2016HO3$omega,dfres_2016HO3$State,mean)
tapply(dfres_2016HO3$omega,dfres_2016HO3$State,sd)

tapply(dfres_2016HO3$e,dfres_2016HO3$State,mean)
tapply(dfres_2016HO3$e,dfres_2016HO3$State,sd)

table(est_states2016HO3)/length(est_states2016HO3)*100
# Identify transitions
sequence=est_states2016HO3
transitions <- table(sequence[-length(sequence)], sequence[-1])
transition_probs <- prop.table(transitions, 1)
print(transition_probs)
run_lengths <- rle(sequence)
avg_duration <- tapply(run_lengths$lengths, run_lengths$values, mean)
print(avg_duration)


table(est_states2016HO3)/length(est_states2016HO3)*100
# Identify transitions
sequence=est_states2016HO3
transitions <- table(sequence[-length(sequence)], sequence[-1])
transition_probs <- prop.table(transitions, 1)
print(transition_probs)
run_lengths <- rle(sequence)
avg_duration <- tapply(run_lengths$lengths, run_lengths$values, mean)
print(avg_duration)


# Combined plot 2002AA29 - 2016HO3 ----------------------------------------

library(ggpubr)
ggarrange(p_theta_res_2002AA29,p_theta_res_2016HO3,
          ncol=1,nrow=2,common.legend = T)


# 2015XX169 ---------------------------------------------------------------

df2015XX169=read.table("propagation_2015XX169_new_v2.txt",header=T)
summary(diff(df2015XX169$t))
unique(df2015XX169$type)

# Step 1000

#df2015XX169=trans_theta(df2015XX169)
summary(df2015XX169)

ground_truth_2015XX169=df2015XX169$type
times_2015XX169=df2015XX169$t
df2015XX169=subset(df2015XX169,select=-c(i,u,type))
names(df2015XX169)

data=df2015XX169
max_lag=10
min_lag=5
tt_thres_diffmaxmin=pi/4

data_thetavol=thetavol_feat(data,
                            wdn=c(min_lag,max_lag))
data_a=a_feat(data,l3=c(min_lag,max_lag))

data_maxmin=max_min_feat(data,
                         l=5,
                         tt_thres_maxmin = 3,
                         l2=c(min_lag,max_lag),
                         tt_thres_diffmaxmin=tt_thres_diffmaxmin)
names(data_maxmin)

data_maxmin_full=data_maxmin
data_maxmin=data_maxmin[,c("t","mean_osc")]
data_fin=merge(data_thetavol,data_a,by="t")
data_fin=merge(data_fin,data_maxmin,by="t")
data_fin=merge(data_fin,data[,c("t","a","theta")],by="t")

names(data_fin)

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]
N=dim(Y)[1]

plot(Y$mean_osc,type='l')

lambda=c(0,5,10,15,20,30)
K=3
kappa=seq(1,ceiling(sqrt(dim(Y)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y)
Lnsat=sat_mod$Lnsat
start_2015XX169=Sys.time()
est2015XX169 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_2015XX169=Sys.time()
elapsed_2015XX169=end_2015XX169-start_2015XX169
save(data_fin,Y,est2015XX169,file="3_2015XX169.RData")

modsel2015XX169=data.frame(hp,
                           FTIC=unlist(lapply(est2015XX169,function(x)x$FTIC))
)

best_mod=modsel2015XX169[which.min(modsel2015XX169$FTIC),]
best_mod

sel=7
#sel=68
estw2015XX169=data.frame(var=colnames(Y),
                         weight=est2015XX169[[sel]]$est_weights)

estw2015XX169=estw2015XX169[order(estw2015XX169$weight,decreasing = T),]
head(estw2015XX169,15)

ground_truth_2015XX169=tail(ground_truth_2015XX169,N)
ground_truth_2015XX169=order_states(ground_truth_2015XX169)
est_states_2015XX169=order_states(est2015XX169[[sel]]$est_states)

BAC(ground_truth_2015XX169,est_states_2015XX169)

# 2020PN1 -------------------------------------------------------------------------

df2020PN1=read.table("propagation_2020PN1_new_v2.txt",header=T)
summary(diff(df2020PN1$t))
names(df2020PN1)
unique(df2020PN1$type)

df2020PN1=trans_theta(df2020PN1)

ground_truth_2020PN1=df2020PN1$type
times_2020PN1=df2020PN1$t
df2020PN1=subset(df2020PN1,select=-c(i,u,type))
names(df2020PN1)

data=df2020PN1
summary(data)
max_lag=10
min_lag=5
tt_thres_diffmaxmin=pi/4

data_thetavol=thetavol_feat(data,
                            wdn=c(min_lag,max_lag))
data_a=a_feat(data,l3=c(min_lag,max_lag))

data_maxmin=max_min_feat(data,
                         l=5,
                         tt_thres_maxmin = 3,
                         l2=c(min_lag,max_lag),
                         tt_thres_diffmaxmin=tt_thres_diffmaxmin)
names(data_maxmin)

data_maxmin_full=data_maxmin
data_maxmin=data_maxmin[,c("t","mean_osc")]
data_fin=merge(data_thetavol,data_a,by="t")
data_fin=merge(data_fin,data_maxmin,by="t")
data_fin=merge(data_fin,data[,c("t","a","theta")],by="t")

names(data_fin)

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]
N=dim(Y)[1]

plot(Y$ind_a_short,type='l')

lambda=c(0,5,10,15,20,30)
K=4
kappa=seq(1,ceiling(sqrt(dim(Y)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y)
Lnsat=sat_mod$Lnsat
start_2020PN1=Sys.time()
est2020PN1 <- parallel::mclapply(1:nrow(hp),
                                 function(x)
                                   SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   df=Y,
                                                   Lnsat=Lnsat),
                                 mc.cores = parallel::detectCores()-1)

end_2020PN1=Sys.time()
elapsed_2020PN1=end_2020PN1-start_2020PN1
save(data_fin,Y,est2020PN1,file="3_2020PN1.RData")

modsel2020PN1=data.frame(hp,
                         FTIC=unlist(lapply(est2020PN1,function(x)x$FTIC))
)

best_mod=modsel2020PN1[which.min(modsel2020PN1$FTIC),]
best_mod

sel=8
#sel=68
estw2020PN1=data.frame(var=colnames(Y),
                       weight=est2020PN1[[sel]]$est_weights)

estw2020PN1=estw2020PN1[order(estw2020PN1$weight,decreasing = T),]
head(estw2020PN1,15)

ground_truth_2020PN1=tail(ground_truth_2020PN1,N)
ground_truth_2020PN1=order_states(ground_truth_2020PN1)
est_states_2020PN1=order_states(est2020PN1[[sel]]$est_states)

BAC(ground_truth_2020PN1,est_states_2020PN1)

plot(Y$theta,type='l')
plot(ground_truth_2020PN1,type='l')
lines(est_states_2020PN1,col='red')


# 2015SO2 -----------------------------------------------------------------

df2015SO2=read.table("propagation_2015SO2_new_v2.txt",header=T)
unique(df2015SO2$type)

# 100006174 --------------------------------------------------------------

source("Utils5_.R")
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

#max_lag=30
max_lag=10
# min_lag=5
# %FC importante che questa soglia sia abbastanza grande da distinguere QS dal resto
tt_thres_diffmaxmin=pi/4

# data_thetavol=thetavol_feat(inv_trans_theta(data),
#                             wdn=c(min_lag,max_lag))
# 
# data_a=a_feat(data,l3=c(min_lag,max_lag))

data_thetavol=thetavol_feat(inv_trans_theta(data),
                            wdn=max_lag)



data_a=a_feat(data,l3=max_lag)

# data_maxmin=max_min_feat(data,
#                          l=3,
#                          tt_thres_maxmin = 2.4,
#                          l2=c(min_lag,max_lag),
#                          tt_thres_diffmaxmin=tt_thres_diffmaxmin)
data_maxmin=max_min_feat(data,
                         l=3,
                         tt_thres_maxmin = 2.4,
                         l2=max_lag,
                         tt_thres_diffmaxmin=tt_thres_diffmaxmin)
names(data_maxmin)

data_maxmin_full=data_maxmin
data_maxmin=data_maxmin[,c("t","I_TD","I_HS","I_QS","mean_osc")]

data_fin=merge(data_thetavol,data_a,by="t")
data_fin=merge(data_fin,data_maxmin,by="t")
data_fin=merge(data_fin,inv_trans_theta(data[,c("t","a","theta")]),by="t")

names(data_fin)
summary(data_fin)

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]

lambda=c(0,5,10,15,20,30)
# K=2:6
# %FC Con K da 2 a 6 seleziona K=3
# K=5
# %FC Provato a fissare K=5 ma non vede compound
# K=4
K=4

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

sel=10
#sel=68
estw100006174=data.frame(var=colnames(Y),
                         weight=est100006174[[sel]]$est_weights)

estw100006174=estw100006174[order(estw100006174$weight,decreasing = T),]
head(estw100006174,15)

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

#df_segments_a$zoom_group <- ifelse(df_segments_a$t < df_segments_a$t[dim(Y)[1] / 2], "First Half", "Second Half")

zoom=300:1600
label_size=18
p_a_100006174=ggplot(data = df_segments_a[zoom,]) + 
  geom_segment(aes(x = t, y = a, 
                   xend = next_t, yend = next_a), 
               size = 1, color = 'grey80') +
  geom_point(aes(x = t, y = a, color = as.factor(State))) +
  scale_color_manual(values = c(4, 1, 2, 3),
                     labels = c("NR", "HS", "QS", "CP"),
                     name = "Orbital regime") +
  scale_x_continuous(labels = label_scientific())+
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(a ~ "(AU)")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "top"                  
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
# +
#   facet_zoom(x = t < df_segments_a$t[dim(Y)[1] / 2], zoom.size = 1) +
#   facet_zoom(x = t >= df_segments_a$t[dim(Y)[1] / 2], zoom.size = 1)

#ggplotly(p_a_res)



df_segments_theta <- df_res_100006174 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))

#df_segments_theta$zoom_group <- ifelse(df_segments_theta$t < df_segments_theta$t[dim(Y)[1] / 2], "First Half", "Second Half")

# Plot con facet_zoom
p_theta_100006174=ggplot(data = df_segments_theta[zoom,]) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1, color = 'grey80') +
  geom_point(aes(x = t, y = theta, color = as.factor(State))) +
  # scale_color_manual(values = c(4, 1, 2, 3),
  #                    labels = c("NR", "HS", "QS", "CP"),
  #                    name = "Orbital regime") +
  scale_x_continuous(labels = label_scientific())+
  scale_y_continuous(
    breaks = c(-pi, 0, pi),  # Specify where to place the labels
    labels = c(expression(-pi), expression(0), expression(pi))  # Use LaTeX-style labels
  )+
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(theta ~ "(rad)")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "top"                  
  ) + guides(color = guide_legend(override.aes = list(size = 5)))
# +
#   facet_zoom(x = t < df_segments_theta$t[dim(Y)[1] / 2], zoom.size = 1) +
#   facet_zoom(x = t >= df_segments_theta$t[dim(Y)[1] / 2], zoom.size = 1)

X11()
ggarrange(p_a_100006174,p_theta_100006174,nrow=2,common.legend = T)

#1:NR 2:HS 3:QS 4:CP

tapply(df_res_100006174$theta,df_res_100006174$State,mean)
tapply(df_res_100006174$theta,df_res_100006174$State,sd)

tapply(c(NA,diff(df_res_100006174$theta)),df_res_100006174$State,mean,na.rm=T)
tapply(c(NA,diff(df_res_100006174$theta)),df_res_100006174$State,sd,na.rm=T)

tapply(df_res_100006174$a,df_res_100006174$State,mean)
tapply(df_res_100006174$a,df_res_100006174$State,sd)

# tapply(df_res_100006174$omega,df_res_100006174$State,mean)
# tapply(df_res_100006174$omega,df_res_100006174$State,sd)
# 
# tapply(df_res_100006174$e,df_res_100006174$State,mean)
# tapply(df_res_100006174$e,df_res_100006174$State,sd)

est_states100006174=est100006174[[sel]]$est_states

round(table(est_states100006174)/length(est_states100006174)*100,2)
# Identify transitions
sequence=est_states100006174
transitions <- table(sequence[-length(sequence)], sequence[-1])
transition_probs <- prop.table(transitions, 1)
print(transition_probs)
run_lengths <- rle(sequence)
avg_duration <- tapply(run_lengths$lengths, run_lengths$values, mean)
print(avg_duration)

# 100011836 ---------------------------------------------------------------

source("Utils5_.R")
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

#max_lag=30
max_lag=10
#min_lag=5
# %FC importante che questa soglia sia abbastanza grande da distinguere QS dal resto
tt_thres_diffmaxmin=pi/4

data_thetavol=thetavol_feat(inv_trans_theta(data),
                            wdn=max_lag)

data_a=a_feat(data,l3=max_lag)

data_maxmin=max_min_feat(data,
                         l=3,
                         tt_thres_maxmin = 2.4,
                         l2=max_lag,
                         tt_thres_diffmaxmin=tt_thres_diffmaxmin)
names(data_maxmin)

data_maxmin_full=data_maxmin
data_maxmin=data_maxmin[,c("t","I_TD","I_HS","I_QS","mean_osc")]

data_fin=merge(data_thetavol,data_a,by="t")
data_fin=merge(data_fin,data_maxmin,by="t")
data_fin=merge(data_fin,inv_trans_theta(data[,c("t","a","theta")]),by="t")

names(data_fin)
summary(data_fin)

ggplotly(ggplot(data_fin)+
           geom_line(aes(x=t,y=theta)))

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]

lambda=c(0,5,10,15,20,30)
#K=2:6
# %FC Con K da 2 a 6 seleziona K=3
# K=5
# %FC Provato a fissare K=5 ma non vede compound
# K=4
K=5

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

sel=15
#sel=68
estw100011836=data.frame(var=colnames(Y),
                         weight=est100011836[[sel]]$est_weights)

estw100011836=estw100011836[order(estw100011836$weight,decreasing = T),]
estw100011836

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

#df_segments_a$zoom_group <- ifelse(df_segments_a$t < df_segments_a$t[dim(Y)[1] / 2], "First Half", "Second Half")

N=dim(Y)[1]
zoom=400:(N/3)

p_a_100011836=ggplot(data = df_segments_a[zoom,]) + 
  geom_segment(aes(x = t, y = a, 
                   xend = next_t, yend = next_a), 
               size = 1, color = 'grey80') +
  geom_point(aes(x = t, y = a, color = as.factor(State))) +
  scale_color_manual(values = c(4, 2, 3, 1,6),
                     labels = c("NR", "QS", "CP", "HS","TP"),
                     name = "Orbital regime") +
  scale_x_continuous(breaks=seq(range(df_segments_a$t[zoom])[1],
                                range(df_segments_a$t[zoom])[2],
                                length.out=4),
                     labels = label_scientific())+
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(a ~ "(AU)")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "top"                  
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))
# +
#   facet_zoom(x = t < df_segments_a$t[dim(Y)[1] / 2], zoom.size = 1) +
#   facet_zoom(x = t >= df_segments_a$t[dim(Y)[1] / 2], zoom.size = 1)

#ggplotly(p_a_res)

p_a_100011836

df_segments_theta <- df_res_100011836 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))

#df_segments_theta$zoom_group <- ifelse(df_segments_theta$t < df_segments_theta$t[dim(Y)[1] / 2], "First Half", "Second Half")

label_size=18
# Plot con facet_zoom
p_theta_100011836=ggplot(data = df_segments_theta[zoom,]) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1, color = 'grey80') +
  geom_point(aes(x = t, y = theta, color = as.factor(State))) +
  scale_color_manual(values = c(4, 2, 3, 1,6),
                     labels = c("NR", "QS", "CP", "HS","TP"),
                     name = "Orbital regime") +
  labs(title = " ", 
       x = "t (year)", 
       y = bquote(theta ~ "(rad)")) +
  scale_y_continuous(
    breaks = c(-pi, 0, pi),  # Specify where to place the labels
    labels = c(expression(-pi), expression(0), expression(pi))  # Use LaTeX-style labels
  )+
  scale_x_continuous(breaks=seq(range(df_segments_theta$t[zoom])[1],
                                range(df_segments_theta$t[zoom])[2],
                                length.out=4),
                     labels = label_scientific())+
  theme_minimal() +
  theme(
    axis.text = element_text(size = label_size),        
    axis.title = element_text(size = label_size),       
    plot.title = element_text(size = label_size),      
    legend.text = element_text(size = label_size),      
    legend.title = element_text(size = label_size),    
    legend.key.size = unit(1.5, "cm"),          
    legend.position = "top"                  
  )+
  guides(color = guide_legend(override.aes = list(size = 5)))
# +
#   facet_zoom(x = t < df_segments_theta$t[dim(Y)[1] / 2], zoom.size = 1) +
#   facet_zoom(x = t >= df_segments_theta$t[dim(Y)[1] / 2], zoom.size = 1)

ggarrange(p_a_100011836,p_theta_100011836,nrow=2,common.legend = T)



#4:NR, 2:QS, 3:CP, 1:HS, 5:TP


tapply(df_res_100011836$theta,df_res_100011836$State,mean)
tapply(df_res_100011836$theta,df_res_100011836$State,sd)

tapply(c(NA,diff(df_res_100011836$theta)),df_res_100011836$State,mean,na.rm=T)
tapply(c(NA,diff(df_res_100011836$theta)),df_res_100011836$State,sd,na.rm=T)

tapply(df_res_100011836$a,df_res_100011836$State,mean)
tapply(df_res_100011836$a,df_res_100011836$State,sd)

# tapply(df_res_100011836$omega,df_res_100011836$State,mean)
# tapply(df_res_100011836$omega,df_res_100011836$State,sd)
# 
# tapply(df_res_100011836$e,df_res_100011836$State,mean)
# tapply(df_res_100011836$e,df_res_100011836$State,sd)

est_states100011836=est100011836[[sel]]$est_states

round(table(est_states100011836)/length(est_states100011836)*100,2)
# Identify transitions
sequence=est_states100011836
transitions <- table(sequence[-length(sequence)], sequence[-1])
transition_probs <- prop.table(transitions, 1)
print(transition_probs)
run_lengths <- rle(sequence)
avg_duration <- tapply(run_lengths$lengths, run_lengths$values, mean)
print(avg_duration*7500)





