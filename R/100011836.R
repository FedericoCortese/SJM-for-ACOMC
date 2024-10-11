source("Utils3_.R")
library(plotly)

df100011836=read.table("100011836.txt",header = T)

start=which.min(abs(df100011836$t-7.9e8))
end=which.min(abs(df100011836$t-8.4e8))

df100011836=df100011836[start:end,]

summary(diff(df100011836$t))
plot(diff(df100011836$t),type='p')

df100011836$type=-1

tt_100011836=theta_trans_plot(df100011836, "100011836",l=5,
                              l2=10,
                              l3=30,tt_thres_maxmin=2.5,
                              tt_thres_diffmaxmin=0.25)
d_temp=tt_100011836$data
p <- ggplot(d_temp, aes(x = seq_along(a))) +
  geom_line(aes(y = a), color = 'blue', size = 1) +
  #geom_line(aes(y = I_diffmaxmin), color = 'red', size = 1) +
  geom_line(aes(y = as.numeric(ind_a)), color = 'violet', size = 1)+
  labs(title = "Line Plot of theta_trans and I_diffmaxmin",
       x = "Index",
       y = "Values") +
  theme_minimal()
ggplotly(p)



df100011836_2.2=subset(d_temp,select=c(theta_trans,count_min,count_max,
                                       value_min,
                                       value_max,
                                       diffmaxmin,
                                       I_diffmaxmin
                                       ,
                                       ind_a
                                       ))

df100011836_2.2$ind_a=as.numeric(df100011836_2.2$ind_a)

# Drop belowabove 1
df100011836_1=subset(df100011836,select=-c(t,type))
df100011836_1=compute_feat(df100011836_1,omega=F,e=F,a=F)
df100011836_1.2=subset(df100011836_1,
                       #select=-c(a, e, i, Omega, omega, theta)
                       select=-c(i,Omega, theta)
)

# merge by rownames with by=0

#df100011836_final = merge(df100011836_1.2,df100011836[,c("t","a","e","i","Omega","omega")],by=0, all=TRUE)
df100011836_final = merge(df100011836_1.2,df100011836[,c("t","a")],by=0, all=TRUE)
rownames(df100011836_final) = df100011836_final$Row.names
df100011836_final = df100011836_final[,-1]

df100011836_final = merge(df100011836_final, 
                          df100011836_2.2, by=0, all=TRUE)
rownames(df100011836_final) = df100011836_final$Row.names
df100011836_final = df100011836_final[,-1]

Y100011836_final=df100011836_final[complete.cases(df100011836_final),]
time_100011836=Y100011836_final$t
Y100011836_final=subset(Y100011836_final,select=-t)

lambda=c(0,5,10,15,20,30)
K=3:5
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
save(time_100011836,Y100011836_final,est100011836,file="100011836.RData")

modsel100011836=data.frame(hp,
                           FTIC=unlist(lapply(est100011836,function(x)x$FTIC))
)
best_mod=modsel100011836[which.min(modsel100011836$FTIC),]
best_mod
sel=41
estw100011836=data.frame(var=colnames(Y100011836_final),
                         weight=est100011836[[sel]]$est_weights)

estw100011836=estw100011836[order(estw100011836$weight,decreasing = T),]
head(estw100011836,15)

df_res_100011836=data.frame(t=time_100011836,
                            Y100011836_final,
                            State=est100011836[[sel]]$est_states)

# df_res_100011836 <- df_res_100011836 %>%
#   mutate(xmax = dplyr::lead(t, default = max(t) + 1))  # Define xmax as the next t value or the max t + 1
#
# Create the plot with geom_tile
# p_a_res <- ggplot(df_res_100011836, aes(x = t)) + 
#   geom_line(aes(y = a), color = 'blue', size = 1) + 
#   geom_tile(aes(x = (t + xmax) / 2, width = xmax - t, 
#                 y = mean(a), 
#                 height = max(a), fill = as.factor(State)), 
#             alpha = 0.2) +
#   scale_fill_manual(values = 1:max(df_res_100011836$State)) +
#   labs(title = "Line Plot of a and states", 
#        x = "Index", 
#        y = "Values") +
#   theme_minimal()
# p_a_res

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
  mutate(next_t = dplyr::lead(t), next_theta = dplyr::lead(theta_trans)) %>%
  filter(!is.na(next_t))  # Remove rows where the next point is missing

p_theta_res <- ggplot() + 
  geom_segment(data = df_segments_theta, aes(x = t, y = theta_trans,
                                             xend = next_t, yend = next_theta, color = as.factor(State)), 
               size = 1) +
  scale_color_manual(values = 1:max(df_res_100011836$State)) +
  labs(title = "100011836", 
       x = "Time (t)", 
       y = "Values (theta)") +
  theme_minimal()
ggplotly(p_theta_res)


# # Omega
# df_res_100011836 <- df_res_100011836 %>%
#   mutate(Segment = cumsum(State != lag(State, default = first(State))))
# 
# # Step 2: Create a new dataframe with start and end points for each line segment
# df_segments_omega<- df_res_100011836 %>%
#   group_by(Segment) %>%
#   mutate(next_t = dplyr::lead(t), next_omega = dplyr::lead(omega)) %>%
#   filter(!is.na(next_t))  # Remove rows where the next point is missing
# 
# p_omega_res <- ggplot() + 
#   geom_segment(data = df_segments_omega, aes(x = t, y = omega,
#                                              xend = next_t, yend = next_omega, color = as.factor(State)), 
#                size = 1) +
#   scale_color_manual(values = 1:max(df_res_100011836$State)) +
#   labs(title = "100011836", 
#        x = "Time (t)", 
#        y = "Values (omega)") +
#   theme_minimal()
# ggplotly(p_omega_res)
