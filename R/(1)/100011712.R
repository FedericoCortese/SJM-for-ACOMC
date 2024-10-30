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
df100011712_1=subset(df100011712,select=-c(t,type))
df100011712_1=compute_feat(df100011712_1)
N=dim(df100011712)[1]

lambda=c(0,5,10,15,20,30)
K=3
kappa=seq(1,ceiling(sqrt(dim(df100011712_1)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

sat_mod=SJM_sat(df100011712_1)
Lnsat=sat_mod$Lnsat

start_est100011712=Sys.time()
est100011712 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=df100011712_1,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_est100011712=Sys.time()
elapsed_est100011712=end_est100011712-start_est100011712
save(df100011712_1,est100011712,file="est100011712.RData")

modsel100011712=data.frame(hp,FTIC=unlist(lapply(est100011712,function(x)x$FTIC))
)
sel=24

estw100011712=data.frame(var=colnames(df100011712_1),
                         weight=est100011712[[sel]]$est_weights)

estw100011712=estw100011712[order(estw100011712$weight,decreasing = T),]
head(estw100011712)

plot(df100011712_1$omega,type='l')
lines(est100011712[[sel]]$est_states,col='red')


# Second approach ---------------------------------------------------------


df100011712_2=theta_trans_plot(df100011712, "100011712",l=5,l2=100,
                               l3=70,tt_thres_maxmin=2.5,
                             tt_thres_diffmaxmin=0.25)
library(plotly)
ggplotly(df100011712_2[[2]])

pp=ggplot(df100011712_2$data[1:25000,],aes(x=t))+
  geom_line(aes(y=theta_trans),col='blue')+
  geom_line(aes(y=count_max),col="red")+
  theme_classic()
ggplotly(pp)

Y_100011712_2=select(df100011712_2$data,subset=-c(t,t_orig,type,ind_a))
Y_100011712_2=Y_100011712_2[complete.cases(Y_100011712_2),]

sat_mod=SJM_sat(Y_100011712_2)
Lnsat=sat_mod$Lnsat


lambda=c(0,5,10,15,20,30)
K=3
kappa=seq(1,ceiling(sqrt(dim(Y_100011712_2)[2])),by=1)
hp=expand.grid(K=K,lambda=lambda,kappa=kappa)

start_est100011712=Sys.time()
est100011712_2 <- parallel::mclapply(1:nrow(hp),
                                   function(x)
                                     SJM_lambdakappa(K=hp[x,]$K,lambda=hp[x,]$lambda,
                                                     kappa=hp[x,]$kappa,
                                                     df=Y_100011712_2,
                                                     Lnsat=Lnsat),
                                   mc.cores = parallel::detectCores()-1)

end_est100011712=Sys.time()
elapsed_est100011712=end_est100011712-start_est100011712
save(df100011712_2,est100011712_2,file="est100011712_2.RData")

modsel100011712_2=data.frame(hp,FTIC=unlist(lapply(est100011712_2,function(x)x$FTIC))
)
sel=12

estw100011712=data.frame(var=colnames(Y_100011712_2),
                         weight=est100011712_2[[sel]]$est_weights)

estw100011712=estw100011712[order(estw100011712$weight,decreasing = T),]
head(estw100011712)

df100011712_3=select(df100011712_2$data,subset=-c(t,t_orig,type,ind_a))
df100011712_3=df100011712_3[complete.cases(df100011712_3),]
df100011712_3$t=1:dim(df100011712_3)[1]
df100011712_3$state=est100011712_2[[sel]]$est_states

pp=ggplot(data=df100011712_3,aes(x=t))+
  geom_line(aes(y=theta_trans),col='grey40')+
  geom_line(aes(y=state),col='red')
ggplotly(pp)


