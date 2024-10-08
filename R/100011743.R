source("Utils3_.R")
library(plotly)

df100011743=read.table("100011743.txt",header = T)

start=which.min(abs(df100011743$t-1.5e8))
end=which.min(abs(df100011743$t-2.e8))

df100011743=df100011743[start:end,]

summary(diff(df100011743$t))
plot(diff(df100011743$t),type='p')

df100011743$type=-1

tt_100011743=theta_trans_plot(df100011743, "100011743",l=5,
                              l2=10,
                              l3=50,tt_thres_maxmin=2.5,
                              tt_thres_diffmaxmin=0.25)
d_temp=tt_100011743$data
p <- ggplot(d_temp, aes(x = seq_along(theta_trans))) + 
  geom_line(aes(y = theta_trans), color = 'blue', size = 1) + 
  #geom_line(aes(y = I_diffmaxmin), color = 'red', size = 1) + 
  geom_line(aes(y = as.numeric(ind_a)), color = 'violet', size = 1)+
  labs(title = "Line Plot of theta_trans and I_diffmaxmin", 
       x = "Index", 
       y = "Values") +
  theme_minimal()
ggplotly(p)

df100011743_2.2=subset(d_temp,select=c(theta_trans,count_min,count_max,
                                              value_min,
                                              value_max,
                                              diffmaxmin,
                                              I_diffmaxmin,
                                              ind_a))

df100011743_2.2$ind_a=as.numeric(df100011743_2.2$ind_a)

df100011743_1=subset(df100011743,select=-c(t,type))
df100011743_1=compute_feat(df100011743_1)
df100011743_1.2=subset(df100011743_1,
                       select=c(dtheta,sd_w1_dtheta,sd_w2_dtheta,domega,sd_w1_domega,sd_w2_domega))

# merge by rownames with by=0

df100011743_final = merge(df100011743_1.2,df100011743[,c("t","a","e","i","Omega","omega")],by=0, all=TRUE)
rownames(df100011743_final) = df100011743_final$Row.names
df100011743_final = df100011743_final[,-1]

df100011743_final = merge(df100011743_final, 
                          df100011743_2.2, by=0, all=TRUE)
rownames(df100011743_final) = df100011743_final$Row.names
df100011743_final = df100011743_final[,-1]

Y100011743_final=df100011743_final[complete.cases(df100011743_final),]
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
save(time_100011743,Y100011743_final,est100011743,file="100011743.RData")

modsel100011743=data.frame(hp,FTIC=unlist(lapply(est100011743,function(x)x$FTIC))
)
sel=32
estw100011743=data.frame(var=colnames(Y100011743_final),
                         weight=est100011743[[sel]]$est_weights)

estw100011743=estw100011743[order(estw100011743$weight,decreasing = T),]
head(estw100011743)

df_res_100011743=data.frame(t=time_100011743,
                            Y100011743_final,
                            State=est100011743[[sel]]$est_states)

p_a_res <- ggplot(df_res_100011743, aes(x = t)) + 
  geom_line(aes(y = a), color = 'blue', size = 1) + 
  geom_line(aes(y = State), color = 'violet', size = 1)+
  labs(title = "Line Plot of a and states", 
       x = "Index", 
       y = "Values") +
  theme_minimal()
ggplotly(p_a_res)
