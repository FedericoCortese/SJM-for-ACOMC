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
K=2:5
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
save(time_100011836,Y100011836_final,est100011836,file="2_100011836.RData")

modsel100011836=data.frame(hp,
                           FTIC=unlist(lapply(est100011836,function(x)x$FTIC))
)
best_mod=modsel100011836[which.min(modsel100011836$FTIC),]
best_mod

sel=54
#sel=68
estw100011836=data.frame(var=colnames(Y100011836_final),
                         weight=est100011836[[sel]]$est_weights)

estw100011836=estw100011836[order(estw100011836$weight,decreasing = T),]
head(estw100011836,15)

df_res_100011836=data.frame(t=time_100011836,
                            Y100011836_final,
                            State=est100011836[[sel]]$est_states)
df_res_100011836 <- df_res_100011836 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_a <- df_res_100011836 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_a = dplyr::lead(a)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_a_res <- ggplot() + 
  geom_segment(data = df_segments_a, aes(x = t, y = a, xend = next_t, yend = next_a, color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = "Values (a)") +
  theme_minimal()
ggplotly(p_a_res)

# Theta
df_res_100011836 <- df_res_100011836 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
df_segments_theta <- df_res_100011836 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_theta_res <- ggplot() + 
  geom_segment(data = df_segments_theta, aes(x = t, y = theta,
                                             xend = next_t, yend = next_theta, color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)

head(estw100011836,15)

pp <- ggplot(data=df_res_100011836,aes(x=t)) + 
  geom_point(aes(y = sd_w1_dtheta, 
                 color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = " ") +
  theme_minimal()

library(ggpubr)
ggarrange(pp,p_theta_res,nrow=2)


# max min -----------------------------------------------------------------

data=df100011836
l=5
tt_thres_maxmin=2.7

maxs <- as.numeric(peaks(data$theta, span = l))
mins <- as.numeric(peaks(-data$theta, span = l))

# Values above tt_thres_maxmin and below -tt_thres_maxmin are not considered as max and min
maxs[which(data$theta>(tt_thres_maxmin))] <- 0
mins[which(data$theta<(-tt_thres_maxmin))] <- 0

data2=data.frame(data,maxs,mins)

df_res_100011836 <- df_res_100011836 %>%
  mutate(Segment = cumsum(State != lag(State, default = first(State))))

# Step 2: Create a new dataframe with start and end points for each line segment
data2 <- data2  %>%
  mutate(Segment = cumsum(maxs != lag(maxs, default = first(maxs))))

data2 <- data2 %>%
  group_by(Segment) %>%
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_maxs <- ggplot() + 
  geom_segment(data = data2, aes(x = t, y = theta,
                                 xend = next_t, yend = next_theta, color = as.factor(maxs)), 
               size = 1) +
  #scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_maxs)


p_maxs <- ggplot() + 
  geom_point(data = data2, aes(x = t, y = theta, color = as.factor(maxs)), 
             size = 1) +
  #scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_maxs)

p_mins <- ggplot() + 
  geom_point(data = data2, aes(x = t, y = theta, color = as.factor(mins)), 
             size = 1) +
  #scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_mins)


# signal processing -------------------------------------------------------

library(signal)

example=c(.9*pi,pi,-.9*pi)
plot(example,type='l',ylim=c(-2,5))
example_unwrap=unwrap(example)
lines(example_unwrap,col='red')

data$temp=unwrap(data$theta,tol=.5)
par(mfrow=c(1,1))
plot(scale(data$theta),type='l')
lines(scale(temp),col='red')

data$dtemp=c(NA,diff(temp))

p1=ggplot(data=data,aes(x=t))+
  geom_line(aes(y=dtemp),col='red')
p2=ggplot(data=data,aes(x=t))+
  geom_line(aes(y=theta),col='grey40')
library(ggpubr)
ggarrange(p1,p2,nrow=2)
  
