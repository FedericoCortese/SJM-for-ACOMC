df100011836=read.table("100011836.txt",header = T)

start=which.min(abs(df100011836$t-7.9e8))
end=which.min(abs(df100011836$t-8.4e8))

df100011836=df100011836[start:end,]

summary(diff(df100011836$t))
#plot(diff(df100011836$t),type='p')

#df100011836$type=-1
data=df100011836
data_name="100011836"

library(ggplot2)
library(dplyr)
library(plotly)

source("Utils4_.R")

data_thetavol=thetavol_feat(data)
data_a=a_feat(data,l3=c(5,30))
data_maxmin=max_min_feat(data,tt_thres_maxmin = 2.4)
names(data_maxmin)
#data_maxmin=data_maxmin[,c("t","I_TD","I_HS","I_QS","mean_osc")]
data_maxmin=data_maxmin[,c("t","mean_osc")]

data_fin=merge(data_thetavol,data_a,by="t")
data_fin=merge(data_fin,data_maxmin,by="t")

names(data_fin)

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]


# fit ---------------------------------------------------------------------

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
                                                     df=Y,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_100011836=Sys.time()
elapsed_100011836=end_100011836-start_100011836
save(time_100011836,Y100011836_final,est100011836,file="2_100011836.RData")


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