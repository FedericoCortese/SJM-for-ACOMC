source("Utils3_.R")
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

df100006174=df100006174[,c("t","a","theta")]

WW=c(10,30)
tt_100006174=feat_comput_theta_a(df100006174, "100006174",
                                 l=5,
                                 wdn=WW,
                                 l2=WW,
                                 l3=WW,
                                 tt_thres_maxmin=2.7,
                                 tt_thres_diffmaxmin=0.25)

Y100006174_final=tt_100006174[complete.cases(tt_100006174),]
time_100006174=Y100006174_final$t
Y100006174_final=subset(Y100006174_final,select=-t)

lambda=c(0,5,10,15,20,30)
K=2:5
kappa=seq(1,ceiling(sqrt(dim(Y100006174_final)[2])),by=1.5)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y100006174_final)
Lnsat=sat_mod$Lnsat
start_100006174=Sys.time()
est100006174 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y100006174_final,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100006174=Sys.time()
elapsed_100006174=end_100006174-start_100006174
save(time_100006174,Y100006174_final,est100006174,file="2_100006174.RData")

modsel100006174=data.frame(hp,
                           FTIC=unlist(lapply(est100006174,function(x)x$FTIC))
)
best_mod=modsel100006174[which.min(modsel100006174$FTIC),]
best_mod

sel=38
estw100006174=data.frame(var=colnames(Y100006174_final),
                         weight=est100006174[[sel]]$est_weights)

estw100006174=estw100006174[order(estw100006174$weight,decreasing = T),]
head(estw100006174,15)

df_res_100006174=data.frame(t=time_100006174,
                            Y100006174_final,
                            State=est100006174[[sel]]$est_states)
df_res_100006174 <- df_res_100006174 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100006174 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_a_res <- ggplot() + 
  geom_segment(data = df_segments_a, aes(x = t, y = a, xend = next_t, yend = next_a, color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100006174$State)) +
  labs(title = "100006174", 
       x = "Time (t)", 
       y = "Values (a)") +
  theme_minimal()
ggplotly(p_a_res)

# Theta
df_res_100006174 <- df_res_100006174 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_theta <- df_res_100006174 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_theta_res <- ggplot() + 
  geom_segment(data = df_segments_theta, aes(x = t, y = theta,
                                             xend = next_t, yend = next_theta, color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100006174$State)) +
  labs(title = "100006174", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)

head(estw100006174,15)

pp <- ggplot(data=df_res_100006174,aes(x=t)) + 
  geom_point(aes(y = sd_w1_dtheta, 
                 color = as.factor(State)), 
             size = 1) +
  scale_color_manual(values = 1:max(df_res_100006174$State)) +
  labs(title = "100006174", 
       x = "Time (t)", 
       y = " ") +
  theme_minimal()

library(ggpubr)
ggarrange(pp,p_theta_res,nrow=2)
