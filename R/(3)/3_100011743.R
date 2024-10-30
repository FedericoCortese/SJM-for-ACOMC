df100011743=read.table("100011743.txt",header = T)

start=which.min(abs(df100011743$t-1.5e8))
end=which.min(abs(df100011743$t-2.e8))

df100011743=df100011743[start:end,]

summary(diff(df100011743$t))
#plot(diff(df100011743$t),type='p')

#df100011743$type=-1
data=df100011743
data_name="100011743"

library(ggplot2)
library(dplyr)
library(plotly)

source("Utils4_.R")

# %FC 30 e' troppo lungo arriva tardi

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

names(data_fin)

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]

# fit ---------------------------------------------------------------------

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
start_100011743=Sys.time()
est100011743 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100011743=Sys.time()
elapsed_100011743=end_100011743-start_100011743
save(data_fin,Y,est100011743,file="3_100011743.RData")

modsel100011743=data.frame(hp,
                           FTIC=unlist(lapply(est100011743,function(x)x$FTIC))
)

best_mod=modsel100011743[which.min(modsel100011743$FTIC),]
best_mod

sel=11
#sel=68
estw100011743=data.frame(var=colnames(Y),
                         weight=est100011743[[sel]]$est_weights)

estw100011743=estw100011743[order(estw100011743$weight,decreasing = T),]
head(estw100011743,15)

# Merge with a theta and t for plotting
data_a_theta=tail(data[,c("t","a","theta")],dim(Y)[1])

df_res_100011743=data.frame(data_a_theta,
                            Y,
                            State=est100011743[[sel]]$est_states
)

df_res_100011743 <- df_res_100011743 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100011743 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a))

p_a_res <- ggplot(data = df_segments_a) + 
  geom_segment(aes(x = t, y = a, 
                   xend = next_t, yend = next_a), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=a,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(df_res_100011743$State)) +
  labs(title = "100011743", 
       x = "Time (t)", 
       y = "Values (a)") +
  theme_minimal()
ggplotly(p_a_res)

df_segments_theta <- df_res_100011743 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta))

p_theta_res <- ggplot(data = df_segments_theta) + 
  geom_segment(aes(x = t, y = theta, 
                   xend = next_t, yend = next_theta), 
               size = 1,color='grey80') +
  geom_point(aes(x=t,y=theta,
                 color=as.factor(State)))+
  scale_color_manual(values = 1:max(df_res_100011743$State)) +
  labs(title = "100011743", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)

# Check feat by feat
pp <- ggplot(data=df_res_100011743,aes(x=t)) + 
  geom_point(aes(y = sd_dtheta_long, 
                 color = as.factor(State)), 
             size = 1) +
  scale_color_manual(values = 1:max(df_res_100011743$State)) +
  labs(title = "100011743", 
       x = "Time (t)", 
       y = " ") +
  theme_minimal()

pp

# results -----------------------------------------------------------------



data2 <- data %>%
  mutate(Segment_max = cumsum(maxs != lag(maxs, default = first(maxs))),
         Segment_min = cumsum(mins != lag(mins, default = first(mins))),
         maxs_mins=mins+maxs*2)

data2 <- data2 %>%
  group_by(Segment_max) %>%
  mutate(next_t_max = dplyr::lead(t), next_theta_max = dplyr::lead(theta))

data2 <- data2 %>%
  group_by(Segment_min) %>%
  mutate(next_t_min = dplyr::lead(t), next_theta_min = dplyr::lead(theta))

# Create the plot with grey segments and points colored by the 'maxs' value
p_maxs_mins <- ggplot(data2) + 
  geom_segment(aes(x = t, y = theta, xend = next_t_max, yend = next_theta_max), 
               color = 'grey80', size = 1) +
  geom_segment(aes(x = t, y = theta, xend = next_t_min, yend = next_theta_min),
               color = 'grey80', size = 1) +
  geom_point(aes(x = t, y = theta, color = factor(maxs_mins)), size = 1) +
  scale_color_manual(values = c("0" = "grey50", "1"='blue' ,"2" = "red")) +
  #geom_point(aes(x = t, y = theta, color = factor(mins)), size = 1) +
  #scale_color_manual(values = c("0" = "grey50", "1" = "blue")) +
  labs(title = "Plot of Segment Maxs and Mins", x = "t", y = "Theta") +
  geom_line(aes(x=t,y=mean_osc))+
  theme_minimal() +
  theme(legend.position = "none") 

# Convert the ggplot object to a plotly interactive plot
ggplotly(p_maxs_mins)



# Check max min diff etc --------------------------------------------------

data2=merge(data_fin,df100011743,by="t")
data2=merge(data2,data_maxmin_full[,c("t","maxs","mins","I_diffmaxmin")],by="t")
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
  labs(title = "Plot of Segment Maxs and Mins", x = "t", y = "Theta") +
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

