source("Utils3_.R")
library(plotly)

df100006174=read.table("100006174.txt",header = T)

start=which.min(abs(df100006174$t-1.6e8))
end=which.min(abs(df100006174$t-1.97e8))

df100006174=df100006174[start:end,]

summary(diff(df100006174$t))
plot(diff(df100006174$t),type='p')

df100006174$type=-1

tt_100006174=theta_trans_plot(df100006174, "100006174",l=5,
                              l2=10,
                              l3=30,tt_thres_maxmin=2.5,
                              tt_thres_diffmaxmin=0.25)
d_temp=tt_100006174$data
# p <- ggplot(d_temp, aes(x = seq_along(theta_trans))) +
#   geom_line(aes(y = theta_trans), color = 'blue', size = 1) +
#   #geom_line(aes(y = I_diffmaxmin), color = 'red', size = 1) +
#   geom_line(aes(y = as.numeric(ind_a)), color = 'violet', size = 1)+
#   labs(title = "Line Plot of theta_trans and I_diffmaxmin",
#        x = "Index",
#        y = "Values") +
#   theme_minimal()
# ggplotly(p)

df100006174_2.2=subset(d_temp,select=c(theta_trans,count_min,count_max,
                                       value_min,
                                       value_max,
                                       diffmaxmin,
                                       I_diffmaxmin,
                                       ind_a))

df100006174_2.2$ind_a=as.numeric(df100006174_2.2$ind_a)

df100006174_1=subset(df100006174,select=-c(t,type))
df100006174_1=compute_feat(df100006174_1,omega=F,e=F,a=F)
df100006174_1.2=subset(df100006174_1,
                       #select=-c(a, e, i, Omega, omega, theta)
                       select=-c(i,Omega, theta)
)

# merge by rownames with by=0

#df100006174_final = merge(df100006174_1.2,df100006174[,c("t","a","e","i","Omega","omega")],by=0, all=TRUE)

df100006174_final = merge(df100006174_1.2,df100006174[,c("t","a")],by=0, all=TRUE)
rownames(df100006174_final) = df100006174_final$Row.names
df100006174_final = df100006174_final[,-1]

df100006174_final = merge(df100006174_final, 
                          df100006174_2.2, by=0, all=TRUE)
rownames(df100006174_final) = df100006174_final$Row.names
df100006174_final = df100006174_final[,-1]

Y100006174_final=df100006174_final[complete.cases(df100006174_final),]
time_100006174=Y100006174_final$t
Y100006174_final=subset(Y100006174_final,select=-t)

lambda=c(0,5,10,15,20,30)
K=3:5
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
save(time_100006174,Y100006174_final,est100006174,file="100006174.RData")

modsel100006174=data.frame(hp,
                           FTIC=unlist(lapply(est100006174,function(x)x$FTIC))
)
best_mod=modsel100006174[which.min(modsel100006174$FTIC),]
best_mod
sel=26
estw100006174=data.frame(var=colnames(Y100006174_final),
                         weight=est100006174[[sel]]$est_weights)

estw100006174=estw100006174[order(estw100006174$weight,decreasing = T),]
head(estw100006174,15)

df_res_100006174=data.frame(t=time_100006174,
                            Y100006174_final,
                            State=est100006174[[sel]]$est_states)

# df_res_100006174 <- df_res_100006174 %>%
#   mutate(xmax = dplyr::lead(t, default = max(t) + 1))  # Define xmax as the next t value or the max t + 1
#
# Create the plot with geom_tile
# p_a_res <- ggplot(df_res_100006174, aes(x = t)) + 
#   geom_line(aes(y = a), color = 'blue', size = 1) + 
#   geom_tile(aes(x = (t + xmax) / 2, width = xmax - t, 
#                 y = mean(a), 
#                 height = max(a), fill = as.factor(State)), 
#             alpha = 0.2) +
#   scale_fill_manual(values = 1:max(df_res_100006174$State)) +
#   labs(title = "Line Plot of a and states", 
#        x = "Index", 
#        y = "Values") +
#   theme_minimal()
# p_a_res

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
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta_trans)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_theta_res <- ggplot() + 
  geom_segment(data = df_segments_theta, aes(x = t, y = theta_trans,
                                             xend = next_t, yend = next_theta, color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100006174$State)) +
  labs(title = "100006174", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)


# Omega
# df_res_100006174 <- df_res_100006174 %>%
#   mutate(Segment = cumsum(State != lag(State, default = first(State))))
# 
# # Step 2: Create a new dataframe with start and end points for each line segment
# df_segments_omega<- df_res_100006174 %>%
#   group_by(Segment) %>%
#   mutate(next_t = dplyr::lead(t), next_omega = dplyr::lead(omega)) %>%
#   filter(!is.na(next_t))  # Remove rows where the next point is missing
# 
# p_omega_res <- ggplot() + 
#   geom_segment(data = df_segments_omega, aes(x = t, y = omega,
#                                              xend = next_t, yend = next_omega, color = as.factor(State)), 
#                size = 1) +
#   scale_color_manual(values = 1:max(df_res_100006174$State)) +
#   labs(title = "100006174", 
#        x = "Time (t)", 
#        y = "Values (omega)") +
#   theme_minimal()
# ggplotly(p_omega_res)
