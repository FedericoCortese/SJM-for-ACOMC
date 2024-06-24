source("Utils3_.R")
df164207=read.table("propagation_164207_new_v2.txt",header=T)


df=df164207[,c("a","e","theta","omega")]
p164207=plot_real(df164207,"blue")
p164207$Pa
p164207$Pe
p164207$Ptheta
p164207$Pomega

# Compute features
df_feat=compute_feat(df,wdn=100)

# Estimate saturated model
est_sat=SJM_sat(df_feat)

# Estimate SJM model
N=dim(df_feat)[1]
true_states=tail(df164207$type,N)

est_lk=SJM_lambdakappa(lambda=15,kappa=4,df=df_feat,Lnsat=est_sat$Lnsat,true_states=true_states)
est_lk2=SJM_lambdakappa(lambda=5,kappa=4,df=df_feat,Lnsat=est_sat$Lnsat,true_states=true_states)
est_lk3=SJM_lambdakappa(lambda=15,kappa=6,df=df_feat,Lnsat=est_sat$Lnsat,true_states=true_states)

est_lk$ARI;est_lk$overlap;est_lk$FTIC
est_lk2$ARI;est_lk2$overlap;est_lk2$FTIC
est_lk3$ARI;est_lk3$overlap;est_lk3$FTIC

p1=plot_real(data.frame(df_feat,type=est_lk$est_states))
p1$Pa
p1$Pe
p1$Ptheta
p1$Pomega
