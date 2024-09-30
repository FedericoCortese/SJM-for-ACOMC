source("Utils3_.R")

df100011712=read.table("100011712.txt",header = T)

start=which(df100011712$t==99619763.30)
end=which(df100011712$t==1000114901.1)

df100011712=df100011712[start:end,]

summary(diff(df100011712$t))
plot(diff(df100011712$t),type='p')

ind_max=tail(which(diff(df100011712$t)==max(diff(df100011712$t))),1)
ind_min=tail(which(diff(df100011712$t)==min(diff(df100011712$t))),1)
df100011712=df100011712[-(1:ind_min),]    

summary(diff(df100011712$t))

df100011712$type=-1

library(plotly)
tt_100011712=theta_trans_plot(df100011712,"100011712",tt_thres_maxmin=1.5,tt_thres_diffmaxmin = 0.05)

plot(tt_100011712$data$theta_trans,type='l')

d_temp=tt_100011712$data

library(plotly)
p <- ggplot(d_temp, aes(x = seq_along(theta_trans))) + 
  geom_line(aes(y = theta_trans), color = 'blue', size = 1) + 
  geom_line(aes(y = I_diffmaxmin), color = 'red', size = 1) + 
  labs(title = "Line Plot of theta_trans and I_diffmaxmin", 
       x = "Index", 
       y = "Values") +
  theme_minimal()
ggplotly(p)

# First approach (Berlin) -------------------------------------------------
df100011712_1=compute_feat(df100011712)
N=dim(df100011712)[1]

lambda=c(0,5,10,15,20,30)
K=c(2,3,4)
kappa=seq(1,ceiling(sqrt(dim(df100011712_1)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(df100011712_1[,-1])
Lnsat=sat_mod$Lnsat

start_est100011712=Sys.time()
est100011712 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=df100011712_1[,-1],
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_est100011712=Sys.time()
elapsed_est100011712=end_est100011712-start_est100011712
save(df100011712_1,est100011712,file="est100011712.RData")
