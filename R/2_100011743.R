source("Utils3_.R")
library(plotly)

df100011743=read.table("100011743.txt",header = T)

start=which.min(abs(df100011743$t-1.5e8))
end=which.min(abs(df100011743$t-2.e8))

df100011743=df100011743[start:end,]

summary(diff(df100011743$t))
plot(diff(df100011743$t),type='p')

#df100011743$type=-1
data=df100011743
data_name="100011743"

names(df100011743)

df100011743=df100011743[,c("t","a","theta")]

WW=c(10,30)
#WW=10
tt_100011743=feat_comput_theta_a(df100011743, "100011743",
                                 l=5,
                                 wdn=WW,
                                 l2=WW,
                                 l3=WW,
                                 tt_thres_maxmin=2.5,
                                 tt_thres_diffmaxmin=0.25)

Y100011743_final=tt_100011743[complete.cases(tt_100011743),]
time_100011743=Y100011743_final$t
Y100011743_final=subset(Y100011743_final,select=-t)

lambda=c(0,5,10,15,20,30)
K=3:5
kappa=seq(1,ceiling(sqrt(dim(Y100011743_final)[2])),by=1.5)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(Y100011743_final)
Lnsat=sat_mod$Lnsat
start_100011743=Sys.time()
est100011743 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y100011743_final,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100011743=Sys.time()
elapsed_100011743=end_100011743-start_100011743
save(time_100011743,Y100011743_final,est100011743,file="2_100011743.RData")

modsel100011743=data.frame(hp,
                           FTIC=unlist(lapply(est100011743,function(x)x$FTIC))
)
best_mod=modsel100011743[which.min(modsel100011743$FTIC),]
best_mod

sel=25
estw100011743=data.frame(var=colnames(Y100011743_final),
                         weight=est100011743[[sel]]$est_weights)

estw100011743=estw100011743[order(estw100011743$weight,decreasing = T),]
head(estw100011743,15)

df_res_100011743=data.frame(t=time_100011743,
                            Y100011743_final,
                            State=est100011743[[sel]]$est_states)
df_res_100011743 <- df_res_100011743 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100011743 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_a_res <- ggplot() + 
  geom_segment(data = df_segments_a, aes(x = t, y = a, xend = next_t, yend = next_a, color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100011743$State)) +
  labs(title = "100011743", 
       x = "Time (t)", 
       y = "Values (a)") +
  theme_minimal()
ggplotly(p_a_res)

# Theta
df_res_100011743 <- df_res_100011743 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_theta <- df_res_100011743 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_theta_res <- ggplot() + 
  geom_segment(data = df_segments_theta, aes(x = t, y = theta,
                                             xend = next_t, yend = next_theta, color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100011743$State)) +
  labs(title = "100011743", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)

head(estw100011743,15)

pp <- ggplot(data=df_res_100011743,aes(x=t)) + 
  geom_point(aes(y = sd_dtheta_long, 
                 color = as.factor(State)), 
             size = 1) +
  scale_color_manual(values = 1:max(df_res_100011743$State)) +
  labs(title = "100011743", 
       x = "Time (t)", 
       y = " ") +
  theme_minimal()

library(ggpubr)
ggarrange(pp,p_theta_res,nrow=2)
