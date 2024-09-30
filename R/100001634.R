source("Utils3_.R")

df100001634=read.table("100001634.txt",header = T)

start=which(df100001634$t==77540640.20)
end=which(df100001634$t==198815647.70)

df100001634=df100001634[start:end,]



# EDA ---------------------------------------------------------------------
summary(diff(df100001634$t))

plot(diff(df100001634$t))

ind_max=which.max(diff(df100001634$t))
ind_min=which.min(diff(df100001634$t))

# drop outliers
df100001634=df100001634[-(1:ind_min),]          

df100001634$type=-1

library(plotly)
tt_100001634=theta_trans_plot(df100001634,"100001634",tt_thres_maxmin=1.5,tt_thres_diffmaxmin = 0.05)

tt_100001634$P_theta
ggplotly(tt_100001634$P_thetatrans)
ggplotly(tt_100001634$P_count_min)
ggplotly(tt_100001634$P_count_max)

ggplotly(tt_100001634$P_I_diff_maxmin)

d_temp=tt_100001634$data
d_temp=d_temp[complete.cases(d_temp),]

#d_temp$I_diffmaxmin=as.factor(d_temp$I_diffmaxmin)
d_temp$dtt=c(NA,diff(d_temp$theta_trans))

ggplotly(
  ggplot(d_temp, aes(x = t, y = theta_trans, color = I_diffmaxmin)) +
    geom_point() +
    theme_minimal()
)

library(plotly)
p <- ggplot(d_temp, aes(x = seq_along(theta_trans))) + 
  geom_line(aes(y = theta_trans), color = 'blue', size = 1) + 
    geom_line(aes(y = I_diffmaxmin), color = 'red', size = 1) + 
    labs(title = "Line Plot of theta_trans and I_diffmaxmin", 
       x = "Index", 
       y = "Values") +
  theme_minimal()
ggplotly(p)

ggplotly(ggplot(d_temp, aes(x = t, y = theta_trans))+
  geom_line()+
  theme_bw()
)




# First approach (Berlin) -------------------------------------------------

df100001634_1=compute_feat(df100001634)
N=dim(df100001634)[1]

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df100001634_1)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(df100001634_1)
Lnsat=sat_mod$Lnsat

start_est164207=Sys.time()
est100001634 <- parallel::mclapply(1:nrow(hp),
                                function(x)
                                  SJM_lambdakappa(lambda=hp[x,]$lambda,
                                                  kappa=hp[x,]$kappa,
                                                  df=df100001634_1,
                                                  Lnsat=Lnsat,
                                                  true_states=true_states),
                                mc.cores = parallel::detectCores()-1)

end_est100001634=Sys.time()
elapsed_est100001634=end_est100001634-start_est100001634
save(df100001634,est100001634,elapsed_est100001634,file="est100001634.RData")