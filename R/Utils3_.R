library(ggplot2)
library(zoo)

order_states=function(states){
  
  # This function organizes states by assigning 1 to the first observed state and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  # states is a vector of observed states
  
  N=length(states)
  states_temp=rep(0,N)
  new=1
  states_temp[1]=new
  for(i in 2:N){
    if(sum(states[i]==states[1:(i-1)])==0){
      # we enter this is if-stat. whenever a new state appeares
      states_temp[i]=new+1
      new=new+1
    }
    else{
      states_temp[i]=states_temp[which(states[1:(i-1)]==states[i])[1]]
    }
  }
  return(states_temp)
}

order_states_freq=function(states){
  
  # This function organizes states by assigning 1 to the mostly observed state and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  # states is a vector of observed states
  
  states_temp=match(states,names(sort(table(states))))
  
  return(states_temp)
}


SJM_sat=function(df,Ksat=6){
  
  #  df is a data.frame WITHOUT time column
  
  # library(reticulate)
  # import("scipy")
  # source_python('SJ.py')
  # Remove time
  #Y=df[,-1]
  Y=apply(df,2,scale)
  
  res_sat=sparse_jump(Y=as.matrix(Y), n_states=as.integer(Ksat), 
                      max_features=sqrt(dim(Y)[2]), 
                      jump_penalty=0,
                      max_iter=as.integer(10), 
                      tol=1e-4, 
                      n_init=as.integer(10), 
                      verbose=F)
 
  
  est_states_sat=order_states(res_sat[[1]])
  
  Lnsat=sum(get_BCSS(as.matrix(Y),est_states_sat))
  
  return(list(Lnsat=Lnsat,Ksat=Ksat))
  
}


library(reticulate)
import("scipy")
source_python('SJ.py')

SJM_lambdakappa=function(lambda,kappa,K=2,df,Lnsat,Ksat=6,alpha0=NULL,K0=NULL,pers0=.95,
                         true_states=NULL){
  
  #  df is a data.frame WITHOUT time column
  
  if(is.null(K0)){
    K0=K
  }
  
  #Y=df[,-1]
  
  YY=apply(df,2,scale)
  PP=dim(YY)[2]
  
  if(is.null(alpha0)){
    alpha0=PP
  }
  
  N=dim(YY)[1]
  
  res=sparse_jump(Y=as.matrix(YY), 
                  n_states=as.integer(K), 
                  max_features=kappa, 
                  jump_penalty=lambda,
                  max_iter=as.integer(10), 
                  tol=1e-4, 
                  n_init=as.integer(10), 
                  verbose=F)
  
  
  est_weights=res[[2]]
  est_weights=est_weights/sum(est_weights)
  
  est_states=order_states(res[[1]])
  
  indx=which( est_weights!=0)
  XX=YY[,indx]
  ### First version: compute BCSS as L1 norm
  Ln=sum(get_BCSS(as.matrix(XX),est_states))
  
  pen=sum(est_states[1:(N-1)]!=est_states[2:N])
  alphak=length(which(est_weights!=0))
  
  CKp=-log(K)-log(2*pi)*kappa/2
  CKp_sat=-log(Ksat)-log(2*pi)*sqrt(PP)/2
  
  anFTIC=log(log(N))*log(PP)
  anAIC=rep(2,length(N))
  anBIC=log(N)
  
  pen0=(1-pers0)*N*(K0-1)
  
  TotalPenalty=(alpha0+pen0)*K+K0*(alphak-alpha0+pen-pen0) 
  Ln_diff=Lnsat-Ln
  CKp_diff=CKp_sat-CKp
  
  FTIC=2*CKp_diff+(Ln_diff+anFTIC*TotalPenalty)/N
  BIC=2*CKp_diff+(Ln_diff+anBIC*TotalPenalty)/N
  AIC=2*CKp_diff+(Ln_diff+anAIC*TotalPenalty)/N
  
  if(!is.null(true_states)){
    true_states=order_states_freq(true_states)
    true_states=factor(true_states,levels=1:K)
    est_states=order_states_freq(est_states)
    est_states=factor(est_states,levels=1:K)
    BAC=caret::confusionMatrix(true_states,est_states)$overall[1]
    #overlap=sum(true_states==est_states)/N
    ARI=pdfCluster::adj.rand.index(true_states,est_states)
    return(list(FTIC=FTIC,
                BIC=BIC,
                AIC=AIC,
                ARI=ARI,
                #overlap=overlap,
                BAC=BAC,
                est_states=est_states,
                est_weights=est_weights))
  }
  else{
    return(list(FTIC=FTIC,
                BIC=BIC,
                AIC=AIC,
                est_states=est_states,
                est_weights=est_weights))
  }
  
}

compute_feat_old=function(dat,wdn=c(10,75),am1=F,wdn_decomp=10){
  
  # Add first differences for each variable to dat2
  #dat$da=c(NA,diff(dat$a))
  dat$de=c(NA,diff(dat$e))
  dat$dtheta=c(NA,diff(dat$theta))
  dat$domega=c(NA,diff(dat$omega))
  
  ###
  # Add moving average
  
  if(length(wdn)>1){
    #dat$ma_w2_a=rollapply(dat$a, wdn[2], mean, fill=NA)
    dat$ma_w2_e=rollapply(dat$e, wdn[2], mean, fill=NA)
    dat$ma_w2_theta=rollapply(dat$theta, wdn[2], mean, fill=NA)
    dat$ma_w2_omega=rollapply(dat$omega, wdn[2], mean, fill=NA)
    
    #dat$sd_w2_a=rollapply(dat$a, wdn[2], sd, fill=NA)
    dat$sd_w2_e=rollapply(dat$e, wdn[2], sd, fill=NA)
    dat$sd_w2_theta=rollapply(dat$theta, wdn[2], sd, fill=NA)
    dat$sd_w2_omega=rollapply(dat$omega, wdn[2], sd, fill=NA)
    
    #dat$sd_w2_da=rollapply(dat$da, wdn[2], sd, fill=NA)
    dat$sd_w2_de=rollapply(dat$de, wdn[2], sd, fill=NA)
    dat$sd_w2_dtheta=rollapply(dat$dtheta, wdn[2], sd, fill=NA)
    dat$sd_w2_domega=rollapply(dat$domega, wdn[2], sd, fill=NA)
    
    # dat$corr_w2_da_de=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","de")],
    #                                                  width=wdn[2], function(x) cor(x[,1],x[,2]),
    #                                                  by.column=FALSE))
    # dat$corr_w2_da_dtheta=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","dtheta")],
    #                                                      width=wdn[2], function(x) cor(x[,1],x[,2]),
    #                                                      by.column=FALSE))
    # dat$corr_w2_da_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","domega")],
    #                                                      width=wdn[2], function(x) cor(x[,1],x[,2]),
    #                                                      by.column=FALSE))
    # dat$corr_w2_de_dtheta=c(rep(NA,wdn[2]-1),rollapply(dat[,c("de","dtheta")],
    #                                                      width=wdn[2], function(x) cor(x[,1],x[,2]), 
    #                                                      by.column=FALSE))
    # dat$corr_w2_de_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("de","domega")],
    #                                                      width=wdn[2], function(x) cor(x[,1],x[,2]), 
    #                                                      by.column=FALSE))
    # dat$corr_w2_dtheta_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("dtheta","domega")],
    #                                                          width=wdn[2], function(x) cor(x[,1],x[,2]), 
    #                                                          by.column=FALSE))
  }
  
  #dat$ma_w1_a=rollapply(dat$a, wdn[1], mean, fill=NA)
  dat$ma_w1_e=rollapply(dat$e, wdn[1], mean, fill=NA)
  dat$ma_w1_theta=rollapply(dat$theta, wdn[1], mean, fill=NA)
  dat$ma_w1_omega=rollapply(dat$omega, wdn[1], mean, fill=NA)
  
  # Add moving standard deviation for each variable to dat
  #dat$sd_w1_a=rollapply(dat$a, wdn[1], sd, fill=NA)
  dat$sd_w1_e=rollapply(dat$e, wdn[1], sd, fill=NA)
  dat$sd_w1_theta=rollapply(dat$theta, wdn[1], sd, fill=NA)
  dat$sd_w1_omega=rollapply(dat$omega, wdn[1], sd, fill=NA)
  
  # Add moving st. dev. of first differences
  
  # Add moving correlation 
  library(zoo)
  
  # Add moving correlations between first differences
  # dat$corr_w1_da_de=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","de")],
  #                                                   width=wdn[1], function(x) cor(x[,1],x[,2]),
  #                                                   by.column=FALSE))
  # dat$corr_w1_da_dtheta=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","dtheta")],
  #                                                       width=wdn[1], function(x) cor(x[,1],x[,2]),
  #                                                       by.column=FALSE))
  # dat$corr_w1_da_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","domega")],
  #                                                       width=wdn[1], function(x) cor(x[,1],x[,2]),
  #                                                       by.column=FALSE))
  # dat$corr_w1_de_dtheta=c(rep(NA,wdn[1]-1),rollapply(dat[,c("de","dtheta")],
  #                                                       width=wdn[1], function(x) cor(x[,1],x[,2]), 
  #                                                       by.column=FALSE))
  # dat$corr_w1_de_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("de","domega")],
  #                                                       width=wdn[1], function(x) cor(x[,1],x[,2]), 
  #                                                       by.column=FALSE))
  # dat$corr_w1_dtheta_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("dtheta","domega")],
  #                                                           width=wdn[1], function(x) cor(x[,1],x[,2]), 
  #                                                           by.column=FALSE))
  
  if(am1){
    dat$am1=dat$a-1
    dat$belowabove1=as.numeric(
      rollapply(dat$a, wdn[2], function(x){all(x<1)|all(x>1)}, fill=NA)
    )
  }
  
  # tsa=ts(dat$a,frequency = wdn_decomp)
  # tsa.stl <- stl(tsa,
  #                s.window=wdn_decomp,
  #                t.window=NULL,
  #                na.action = na.approx,
  #                robust=T)
  # dat$seas_a=tsa.stl$time.series[,1]
  # dat$trend_a=tsa.stl$time.series[,2]
  # dat$remainder_a=tsa.stl$time.series[,3]
  
  tse=ts(dat$e,frequency = wdn_decomp)
  tse.stl <- stl(tse, 
                 s.window=wdn_decomp,
                 t.window=NULL,
                 na.action = na.approx,
                 robust=T)
  
  dat$seas_e=tse.stl$time.series[,1]
  dat$trend_e=tse.stl$time.series[,2]
  dat$remainder_e=tse.stl$time.series[,3]

  tstheta=ts(dat$theta,frequency = wdn_decomp)
  tstheta.stl <- stl(tstheta, 
                     s.window=wdn_decomp,
                     t.window=NULL,
                     na.action = na.approx,
                     robust=T)
  
  dat$seas_theta=tstheta.stl$time.series[,1]
  dat$trend_theta=tstheta.stl$time.series[,2]
  dat$remainder_theta=tstheta.stl$time.series[,3]
  
  tsomega=ts(dat$omega,frequency = wdn_decomp)
  tsomega.stl <- stl(tsomega, 
                     s.window=wdn_decomp,
                     t.window=NULL,
                     na.action = na.approx,
                     robust=T)
  
  dat$seas_omega=tsomega.stl$time.series[,1]
  dat$trend_omega=tsomega.stl$time.series[,2]
  dat$remainder_omega=tsomega.stl$time.series[,3]
  
  # Add moving standard deviation for trend, seas, and reminder of each variable
  # dat$sd_w1_seas_a=rollapply(dat$seas_a, wdn[1], sd, fill=NA)
  # dat$sd_w1_trend_a=rollapply(dat$trend_a, wdn[1], sd, fill=NA)
  # dat$sd_w1_remainder_a=rollapply(dat$remainder_a, wdn[1], sd, fill=NA)
  
  dat$sd_w1_seas_e=rollapply(dat$seas_e, wdn[1], sd, fill=NA)
  dat$sd_w1_trend_e=rollapply(dat$trend_e, wdn[1], sd, fill=NA)
  dat$sd_w1_remainder_e=rollapply(dat$remainder_e, wdn[1], sd, fill=NA)
  
  dat$sd_w1_seas_theta=rollapply(dat$seas_theta, wdn[1], sd, fill=NA)
  dat$sd_w1_trend_theta=rollapply(dat$trend_theta, wdn[1], sd, fill=NA)
  dat$sd_w1_remainder_theta=rollapply(dat$remainder_theta, wdn[1], sd, fill=NA)

  dat$sd_w1_seas_omega=rollapply(dat$seas_omega, wdn[1], sd, fill=NA)
  dat$sd_w1_trend_omega=rollapply(dat$trend_omega, wdn[1], sd, fill=NA)
  dat$sd_w1_remainder_omega=rollapply(dat$remainder_omega, wdn[1], sd, fill=NA)
  
  if(length(wdn)>1){
    # dat$sd_w2_seas_a=rollapply(dat$seas_a, wdn[2], sd, fill=NA)
    # dat$sd_w2_trend_a=rollapply(dat$trend_a, wdn[2], sd, fill=NA)
    # dat$sd_w2_remainder_a=rollapply(dat$remainder_a, wdn[2], sd, fill=NA)
    
    dat$sd_w2_seas_e=rollapply(dat$seas_e, wdn[2], sd, fill=NA)
    dat$sd_w2_trend_e=rollapply(dat$trend_e, wdn[2], sd, fill=NA)
    dat$sd_w2_remainder_e=rollapply(dat$remainder_e, wdn[2], sd, fill=NA)
    
    dat$sd_w2_seas_theta=rollapply(dat$seas_theta, wdn[2], sd, fill=NA)
    dat$sd_w2_trend_theta=rollapply(dat$trend_theta, wdn[2], sd, fill=NA)
    dat$sd_w2_remainder_theta=rollapply(dat$remainder_theta, wdn[2], sd, fill=NA)
    
    dat$sd_w2_seas_omega=rollapply(dat$seas_omega, wdn[2], sd, fill=NA)
    dat$sd_w2_trend_omega=rollapply(dat$trend_omega, wdn[2], sd, fill=NA)
    dat$sd_w2_remainder_omega=rollapply(dat$remainder_omega, wdn[2], sd, fill=NA)
  }
  
  datC=dat[complete.cases(dat),]
  
  return(datC)
  
}

compute_feat=function(dat,wdn=c(10,75),wdn_decomp=10,
                      a=T,e=T,theta=T,omega=T){
  
  library(zoo)
  if(a){
    # a features
    dat$da=c(NA,diff(dat$a))
    
    tsa=ts(dat$a,frequency = wdn_decomp)
    tsa.stl <- stl(tsa,
                   s.window=wdn_decomp,
                   t.window=NULL,
                   na.action = na.approx,
                   robust=T)
    dat$seas_a=tsa.stl$time.series[,1]
    dat$trend_a=tsa.stl$time.series[,2]
    dat$remainder_a=tsa.stl$time.series[,3]
    
    # Add moving standard deviation for trend, seas, and reminder of each variable
    dat$sd_w1_seas_a=rollapply(dat$seas_a, wdn[1], sd, fill=NA)
    dat$sd_w1_trend_a=rollapply(dat$trend_a, wdn[1], sd, fill=NA)
    dat$sd_w1_remainder_a=rollapply(dat$remainder_a, wdn[1], sd, fill=NA)
    
    dat$ma_w1_a=rollapply(dat$a, wdn[1], mean, fill=NA)
    dat$sd_w1_a=rollapply(dat$a, wdn[1], sd, fill=NA)
    dat$sd_w1_da=rollapply(dat$da, wdn[1], sd, fill=NA)
    
    if(length(wdn)>1){
      dat$ma_w2_a=rollapply(dat$a, wdn[2], mean, fill=NA)
      dat$sd_w2_a=rollapply(dat$a, wdn[2], sd, fill=NA)
      dat$sd_w2_da=rollapply(dat$da, wdn[2], sd, fill=NA)
      
      dat$sd_w2_seas_a=rollapply(dat$seas_a, wdn[2], sd, fill=NA)
      dat$sd_w2_trend_a=rollapply(dat$trend_a, wdn[2], sd, fill=NA)
      dat$sd_w2_remainder_a=rollapply(dat$remainder_a, wdn[2], sd, fill=NA)
    }
    
    dat$am1=dat$a-1
    dat$belowabove1=as.numeric(
      rollapply(dat$a, wdn[1], function(x){all(x<1)|all(x>1)}, fill=NA)
    )
    
  }
  else{
    dat=subset(dat,select=-a)
  }
  
  if(e){
    dat$de=c(NA,diff(dat$e))
    
    dat$ma_w1_e=rollapply(dat$e, wdn[1], mean, fill=NA)
    dat$sd_w1_e=rollapply(dat$e, wdn[1], sd, fill=NA)
    dat$sd_w1_de=rollapply(dat$de, wdn[1], sd, fill=NA)
    
    tse=ts(dat$e,frequency = wdn_decomp)
    tse.stl <- stl(tse, 
                   s.window=wdn_decomp,
                   t.window=NULL,
                   na.action = na.approx,
                   robust=T)
    
    dat$seas_e=tse.stl$time.series[,1]
    dat$trend_e=tse.stl$time.series[,2]
    dat$remainder_e=tse.stl$time.series[,3]
    
    dat$sd_w1_seas_e=rollapply(dat$seas_e, wdn[1], sd, fill=NA)
    dat$sd_w1_trend_e=rollapply(dat$trend_e, wdn[1], sd, fill=NA)
    dat$sd_w1_remainder_e=rollapply(dat$remainder_e, wdn[1], sd, fill=NA)
    
    if(length(wdn)>1){
      dat$ma_w2_e=rollapply(dat$e, wdn[2], mean, fill=NA)
      dat$sd_w2_e=rollapply(dat$e, wdn[2], sd, fill=NA)
      dat$sd_w2_de=rollapply(dat$de, wdn[2], sd, fill=NA)
      
      dat$sd_w2_seas_e=rollapply(dat$seas_e, wdn[2], sd, fill=NA)
      dat$sd_w2_trend_e=rollapply(dat$trend_e, wdn[2], sd, fill=NA)
      dat$sd_w2_remainder_e=rollapply(dat$remainder_e, wdn[2], sd, fill=NA)
    }
    
  }
  else{
    dat=subset(dat,select=-e)
  }
  
  if(theta){
    # theta features
    dat$dtheta=c(NA,diff(dat$theta))
    
    tstheta=ts(dat$theta,frequency = wdn_decomp)
    tstheta.stl <- stl(tstheta, 
                       s.window=wdn_decomp,
                       t.window=NULL,
                       na.action = na.approx,
                       robust=T)
    
    dat$seas_theta=tstheta.stl$time.series[,1]
    dat$trend_theta=tstheta.stl$time.series[,2]
    dat$remainder_theta=tstheta.stl$time.series[,3]
    
    dat$sd_w1_seas_theta=rollapply(dat$seas_theta, wdn[1], sd, fill=NA)
    dat$sd_w1_trend_theta=rollapply(dat$trend_theta, wdn[1], sd, fill=NA)
    dat$sd_w1_remainder_theta=rollapply(dat$remainder_theta, wdn[1], sd, fill=NA)
    
    dat$ma_w1_theta=rollapply(dat$theta, wdn[1], mean, fill=NA)
    dat$sd_w1_theta=rollapply(dat$theta, wdn[1], sd, fill=NA)
    dat$sd_w1_dtheta=rollapply(dat$dtheta, wdn[1], sd, fill=NA)
    
    if(length(wdn)>1){
      dat$ma_w2_theta=rollapply(dat$theta, wdn[2], mean, fill=NA)
      dat$sd_w2_theta=rollapply(dat$theta, wdn[2], sd, fill=NA)
      dat$sd_w2_dtheta=rollapply(dat$dtheta, wdn[2], sd, fill=NA)
      
      dat$sd_w2_seas_theta=rollapply(dat$seas_theta, wdn[2], sd, fill=NA)
      dat$sd_w2_trend_theta=rollapply(dat$trend_theta, wdn[2], sd, fill=NA)
      dat$sd_w2_remainder_theta=rollapply(dat$remainder_theta, wdn[2], sd, fill=NA)
    }
    
  }
  else{
    dat=subset(dat,select=-theta)
  }
  
  if(omega){
    # omega features
    dat$domega=c(NA,diff(dat$omega))
    
    tsomega=ts(dat$omega,frequency = wdn_decomp)
    tsomega.stl <- stl(tsomega, 
                       s.window=wdn_decomp,
                       t.window=NULL,
                       na.action = na.approx,
                       robust=T)
    
    dat$seas_omega=tsomega.stl$time.series[,1]
    dat$trend_omega=tsomega.stl$time.series[,2]
    dat$remainder_omega=tsomega.stl$time.series[,3]
    
    dat$sd_w1_seas_omega=rollapply(dat$seas_omega, wdn[1], sd, fill=NA)
    dat$sd_w1_trend_omega=rollapply(dat$trend_omega, wdn[1], sd, fill=NA)
    dat$sd_w1_remainder_omega=rollapply(dat$remainder_omega, wdn[1], sd, fill=NA)
    
    dat$ma_w1_omega=rollapply(dat$omega, wdn[1], mean, fill=NA)
    dat$sd_w1_omega=rollapply(dat$omega, wdn[1], sd, fill=NA)
    dat$sd_w1_domega=rollapply(dat$domega, wdn[1], sd, fill=NA)
    
    if(length(wdn)>1){
      dat$ma_w2_omega=rollapply(dat$omega, wdn[2], mean, fill=NA)
      dat$sd_w2_omega=rollapply(dat$omega, wdn[2], sd, fill=NA)
      dat$sd_w2_domega=rollapply(dat$domega, wdn[2], sd, fill=NA)
      
      dat$sd_w2_seas_omega=rollapply(dat$seas_omega, wdn[2], sd, fill=NA)
      dat$sd_w2_trend_omega=rollapply(dat$trend_omega, wdn[2], sd, fill=NA)
      dat$sd_w2_remainder_omega=rollapply(dat$remainder_omega, wdn[2], sd, fill=NA)
    }
    
  }
  else{
    dat=subset(dat,select=-omega)
  }
  
  # dat$corr_w2_da_de=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","de")],
  #                                                  width=wdn[2], function(x) cor(x[,1],x[,2]),
  #                                                  by.column=FALSE))
  # dat$corr_w2_da_dtheta=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","dtheta")],
  #                                                      width=wdn[2], function(x) cor(x[,1],x[,2]),
  #                                                      by.column=FALSE))
  # dat$corr_w2_da_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","domega")],
  #                                                      width=wdn[2], function(x) cor(x[,1],x[,2]),
  #                                                      by.column=FALSE))
  # dat$corr_w2_de_dtheta=c(rep(NA,wdn[2]-1),rollapply(dat[,c("de","dtheta")],
  #                                                      width=wdn[2], function(x) cor(x[,1],x[,2]), 
  #                                                      by.column=FALSE))
  # dat$corr_w2_de_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("de","domega")],
  #                                                      width=wdn[2], function(x) cor(x[,1],x[,2]), 
  #                                                      by.column=FALSE))
  # dat$corr_w2_dtheta_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("dtheta","domega")],
  #                                                          width=wdn[2], function(x) cor(x[,1],x[,2]), 
  #                                                          by.column=FALSE))
  
  # Add moving correlations between first differences
  # dat$corr_w1_da_de=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","de")],
  #                                                   width=wdn[1], function(x) cor(x[,1],x[,2]),
  #                                                   by.column=FALSE))
  # dat$corr_w1_da_dtheta=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","dtheta")],
  #                                                       width=wdn[1], function(x) cor(x[,1],x[,2]),
  #                                                       by.column=FALSE))
  # dat$corr_w1_da_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","domega")],
  #                                                       width=wdn[1], function(x) cor(x[,1],x[,2]),
  #                                                       by.column=FALSE))
  # dat$corr_w1_de_dtheta=c(rep(NA,wdn[1]-1),rollapply(dat[,c("de","dtheta")],
  #                                                       width=wdn[1], function(x) cor(x[,1],x[,2]), 
  #                                                       by.column=FALSE))
  # dat$corr_w1_de_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("de","domega")],
  #                                                       width=wdn[1], function(x) cor(x[,1],x[,2]), 
  #                                                       by.column=FALSE))
  # dat$corr_w1_dtheta_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("dtheta","domega")],
  #                                                           width=wdn[1], function(x) cor(x[,1],x[,2]), 
  #                                                           by.column=FALSE))
  
  datC=dat[complete.cases(dat),]
  
  return(datC)
  
}

plot_real=function(df,color="red"){
  
  # This function plots time-series of a, e, theta, and omega highlighting the regimes with different colors
  
  # df is a data.frame with columns t, a, e, theta, omega, and type, and without time
  
  df$t=1:nrow(df)
  
  df$type=factor(df$type)
  levels(df$type)=c(1,2)
  
  back_res=df[df$type == 1,]
  back_res$tend=back_res$t + 1
  
  Pa=ggplot(df, aes(t, a)) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res, fill = color, alpha = 0.15) +
    geom_line(colour = "black") +
    theme_bw()
  
  Pe=ggplot(df, aes(t, e)) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res, fill = color, alpha = 0.15) +
    geom_line(colour = "black") +
    theme_bw()
  
  Ptheta=ggplot(df, aes(t, theta)) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res, fill = color, alpha = 0.15) +
    geom_line(colour = "black") +
    theme_bw()
  
  Pomega=ggplot(df, aes(t, omega)) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res, fill = color, alpha = 0.15) +
    geom_line(colour = "black") +
    theme_bw()
  
  return(list(Pa=Pa,Pe=Pe,Ptheta=Ptheta,Pomega=Pomega))
}

plot_real4=function(df,color=c("red3","orange","blue")){
  
  # This function plots time-series of a, e, theta, and omega highlighting the regimes with different colors
  
  # df is a data.frame with columns t, a, e, theta, omega, and type, and without time
  
  df$t=1:nrow(df)
  
  df$type=factor(df$type)
  levels(df$type)=c(1,2,3,4)
  
  back_res1=df[df$type == 1,]
  back_res1$tend=back_res1$t + 1
  
  back_res2=df[df$type == 2,]
  back_res2$tend=back_res2$t + 1
  
  back_res3=df[df$type == 3,]
  back_res3$tend=back_res3$t + 1
  
  Pa=ggplot(df, aes(t, a)) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res1, fill = color[1], alpha = 0.15) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res2, fill = color[2], alpha = 0.15) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res3, fill = color[3], alpha = 0.15) +
    geom_line(colour = "black") +
    theme_bw()
  
  Pe=ggplot(df, aes(t, e)) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res1, fill = color[1], alpha = 0.15) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res2, fill = color[2], alpha = 0.15) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res3, fill = color[3], alpha = 0.15) +
    geom_line(colour = "black") +
    theme_bw()
  
  Ptheta=ggplot(df, aes(t, theta)) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res1, fill = color[1], alpha = 0.15) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res2, fill = color[2], alpha = 0.15) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res3, fill = color[3], alpha = 0.15) +
    geom_line(colour = "black") +
    theme_bw()
  
  Pomega=ggplot(df, aes(t, omega)) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res1, fill = color[1], alpha = 0.15) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res2, fill = color[2], alpha = 0.15) +
    geom_rect(ymin = -Inf, ymax = Inf, aes(xmin = t, xmax = tend),
              data = back_res3, fill = color[3], alpha = 0.15) +
    geom_line(colour = "black") +
    theme_bw()
  
  return(list(Pa=Pa,Pe=Pe,Ptheta=Ptheta,Pomega=Pomega))
}

DD = function(x,type,out.plot=F){
  
  # This function computes the drawdowns based on local maxima of type 1,2 or 3.
  
  # Arguments:
  # x is the series of prices
  # type is the desired type of local maximum 
  # out.plot is a logical variable to plot the results
  
  # Value:
  # List with the drawdowns, if out.plot is TRUE, it also returns a plot for the drawdowns.
  # tmax is a vector with times of the local maxima, tmin is a vector with times of the local minima.
  
  TT=length(x)
  
  switch(type, 
         '1'={
           tmin = tmax = NULL
           dd = NULL
           for(t in 1:(TT-1)){
             if(x[t]>x[t+1]){
               tmax = c(tmax,t)
               tmin = c(tmin,t+1)
               dd = c(dd,100*(x[t+1]/x[t]-1))
             }
             
           }
         },
         '2'={
           tmin = tmax = NULL
           dd = NULL
           for(t in 2:(TT-1)){
             if(x[t]<x[t-1] & x[t]<x[t+1]){
               tmin = c(tmin,t)
               if(!is.null(tmax)) dd = c(dd,100*(x[t]/x[tpmax]-1))
             }
             if(x[t]>x[t-1] & x[t]>x[t+1]){
               tmax = c(tmax,t)
               tpmax = t
             }
           }
         },
         '3'={
           tmax = NULL
           for(t in 2:(TT-1)){
             if(all(x[t]>x[1:(t-1)]) & x[t]>x[t+1]) tmax = c(tmax,t)
           }
           if(is.null(tmax)){
             tmax=1
           }
           tmin = NULL
           dd = NULL
           if(length(tmax)==1){
             tmin=which.min(x)
             dd=100*(x[tmin]/x[tmax]-1)
           }
           else{
             for(j in 2:length(tmax)){
               t1 = tmax[j-1]+which.min(x[(tmax[j-1]+1):(tmax[j]-1)]) #time of the jth loc min
               tmin = c(tmin,t1)
               dd = c(dd,100*(x[t1]/x[tmax[j]]-1))
             }
           }
         },
         {
           return(print('Not available'))
         }
  )
  if(out.plot){
    maxTime=rep(0,length(x))
    maxTime[tmax]=1
    minTime=rep(0,length(x))
    minTime[tmin]=2
    data=data.frame(Price=x,Time=1:length(x),category=as.factor(maxTime+minTime))
    if(type==1){
      data$category=recode(data$category,"1"="max","2"="min","3"="both")
      plot=ggplot(data=data,aes(x=Time))+
        geom_line(aes(y=Price))+
        geom_point(aes(y=Price,color=category,shape=category),size=3)+
        scale_color_manual(breaks=c("max","min","both"),values = c("0"="grey1","max" = "green3", "min" = "red2","both"="yellow2"))+
        scale_shape_manual(guide="none",values = c(1,16,16,16))+ 
        theme_bw()+ theme(axis.text=element_text(size=18),
                          axis.title = element_text(size=18),
                          legend.position = "top",
                          legend.title = element_blank(),
                          legend.text = element_text(size=16))+ 
        scale_fill_discrete(labels=c("none",'max', 'min', 'both'))
    }
    else{
      data$category=recode(data$category,"1"="max","2"="min")
      plot=ggplot(data=data,aes(x=Time))+
        geom_line(aes(y=Price))+
        geom_point(aes(y=Price,color=category,shape=category),size=3)+
        scale_color_manual(breaks=c("max","min"),values = c("0"="grey1","max" = "green3", "min" = "red2"))+
        scale_shape_manual(guide="none",values = c(1,16,16,16))+ 
        theme_bw()+ theme(axis.text=element_text(size=18),
                          axis.title = element_text(size=18),
                          legend.position = "top",
                          legend.title = element_blank(),
                          legend.text = element_text(size=16))+ 
        scale_fill_discrete(labels=c("none",'max', 'min'))
    }
    
    return(list(dd=dd,tmin=tmin,tmax=tmax,plot=plot))
  }
  else{
    return(list(dd=dd,tmin=tmin,tmax=tmax))
  }
}

last_min_max <- function(last_min, last_max, TT,x) {
  # Initialize vectors to store the last observed min and max for each time step
  observed_min <- numeric(TT)
  observed_max <- numeric(TT)
  last_observed_min <- numeric(TT)
  last_observed_max <- numeric(TT)
  
  # For each time step t from 1 to TT, find the last observed min and max
  for (t in 1:TT) {
    observed_min[t] <- max(last_min[last_min <= t], na.rm = TRUE)
    if(observed_min[t]==-Inf){observed_min[t]=NA}
    last_observed_min[t]=x[observed_min[t]]
    
    observed_max[t] <- max(last_max[last_max <= t], na.rm = TRUE)
    if(observed_max[t]==-Inf){observed_max[t]=NA}
    last_observed_max[t]=x[observed_max[t]]
  }
  
  # Create the dataframe
  df <- data.frame(
    last_time_min = observed_min,
    last_time_max = observed_max,
    last_observed_min,
    last_observed_max
  )
  
  return(df)
}


last_max_value <- function(is_max, values) {
  # Initialize a vector to store the result
  result <- numeric(length(values))
  
  # Initialize the last observed max
  last_max <- NA
  
  # Iterate through the values
  for (i in seq_along(values)) {
    if (is_max[i] == 1) {
      # If the current position is a max, store the value
      last_max <- values[i]
    }
    # Store the last observed max
    result[i] <- last_max
  }
  
  return(result)
}

# Define the function for minima
last_min_value <- function(is_min, values) {
  # Initialize a vector to store the result
  result <- numeric(length(values))
  
  # Initialize the last observed min
  last_min <- NA
  
  # Iterate through the values
  for (i in seq_along(values)) {
    if (is_min[i] == 1) {
      # If the current position is a min, store the value
      last_min <- values[i]
    }
    # Store the last observed min
    result[i] <- last_min
  }
  
  return(result)
}

library(ggplot2)
library(splus2R)
library(runner)

theta_trans_plot <- function(data, data_name,l=5,l2=50,l3=70,tt_thres_maxmin=2.01,
                             tt_thres_diffmaxmin=0.5) {
  
  # Function to get the last observed min and max value
  
  # Arguments:
  # data: dataframe with the data (t theta and type columns are required)
  # data_name: name of the dataset
  # l: span for peaks function
  # l2: span for runner function (most important)
  # l3: span for I(a<1|a>1)
  # tt_thres: threshold for theta (peaks (and valleys) above (below) this value are not considered as peaks (valleys))
  # tt_thres_diffmaxmin: threshold for the difference between max and min
  
  t_orig=data$t
  data$t=seq_along(data$t)
  data$t_orig=t_orig
  data$type <- as.factor(data$type)
  
  # Unique types and total rows
  unique_types <- unique(data$type)
  total_rows <- dim(data)[1]
  
  # Custom colors based on unique 'type'
  custom_colors <- c("-1" = "blue", "0" = "red", "50" = "green", "100" = "gray50")
  
  # Theta transformation
  theta_trans <- data$theta
  theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
  data$theta_trans <- theta_trans
  
  # Plot 0: Scatter plot of theta
  P0=ggplot(data, aes(x = t, y = theta)) +
    geom_point(aes(color = type)) +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_classic()+
    labs(title = data_name)
  
  # Plot 1: Scatter plot of theta_trans
  P1=ggplot(data, aes(x = t, y = theta_trans)) +
    geom_point(aes(color = type)) +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_classic()+
    labs(title = data_name)
  
  # Detecting peaks and valleys
  #l <- 5
  maxs <- as.numeric(peaks(data$theta_trans, span = l))
  mins <- as.numeric(peaks(-data$theta_trans, span = l))
  
  # Above 2.5 and below -2.5 are not considered as max and min
  maxs[which(data$theta_trans>(tt_thres_maxmin))] <- 0
  mins[which(data$theta_trans<(-tt_thres_maxmin))] <- 0
  
  # Smoothing with runner function
  #l2 <- 50
  count_min <- runner::runner(mins, k = l2, f = sum, na_pad = TRUE)
  count_min <- c(count_min[-(1:floor(l2 / 2))], rep(NA, round(l2 / 2)))
  
  count_max <- runner::runner(maxs, k = l2, f = sum, na_pad = TRUE)
  count_max <- c(count_max[-(1:floor(l2 / 2))], rep(NA, round(l2 / 2)))

  # dat_min <- data.frame(t = data$t, count_min = count_min, type = data$type)
  # dat_max <- data.frame(t = data$t, count_max = count_max, type = data$type)
  # dat_max_min <- merge(dat_max, dat_min, by = c('t', 'type'))
  
  data$count_min <- count_min
  data$count_max <- count_max
  
  # Plot 2: Scatter plot of count_min
  P2=ggplot(data, aes(x = t, y = count_min, color = type)) +
    geom_point() +
    labs(x = "t", y = "count_min") +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_minimal() +
    theme(legend.position = "top")+
    labs(title = paste(data_name," lag = ", l2))
  
  # Plot 3: Scatter plot of count_max
  P3=ggplot(data, aes(x = t, y = count_max, color = type)) +
    geom_point() +
    labs(x = "t", y = "count_max") +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_minimal() +
    theme(legend.position = "top")+
    labs(title = paste(data_name," lag = ", l2))
  
  # Add value_min, value_max and diffmaxmin
  data$value_min <- last_min_value(mins, data$theta_trans)
  data$value_max <- last_max_value(maxs, data$theta_trans)
  data$diffmaxmin <- data$value_max - data$value_min
  
  data$I_diffmaxmin = as.numeric(data$diffmaxmin>tt_thres_diffmaxmin)
  
  # Plot 4: Scatter plot of value_min
  P4=ggplot(data, aes(x = t, y = value_min, color = type)) +
    geom_point() +
    labs(x = "t", y = "value_min") +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_minimal() +
    theme(legend.position = "top")+
    labs(title = paste(data_name," lag = ", l2))
  
  # Plot 5: Scatter plot of value_max
  P5=ggplot(data, aes(x = t, y = value_max, color = type)) +
    geom_point() +
    labs(x = "t", y = "value_max") +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_minimal() +
    theme(legend.position = "top")+
    labs(title = paste(data_name," lag = ", l2))
  
  # Plot 6: Scatter plot of diffmaxmin
  P6=ggplot(data, aes(x = t, y = diffmaxmin, color = type)) +
    geom_point() +
    labs(x = "t", y = "diffmaxmin") +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_minimal() +
    theme(legend.position = "top")+
    labs(title = paste(data_name," lag = ", l2))
  
  P6.1=ggplot(data, aes(x = t, y = I_diffmaxmin, color = type)) +
    geom_point() +
    labs(x = "t", y = paste("I(diffmaxmin>", tt_thres_diffmaxmin, ")")) +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_minimal() +
    theme(legend.position = "top")+
    labs(title = paste(data_name," lag = ", l2))
  
  
  cust_fun=function(x){
    all(x==T)
  }
  
  ind_a=I(data$a<1)
  ind_a2=I(data$a>1)
  ind_a_mov <- runner::runner(ind_a, k = l3, f = cust_fun, na_pad = TRUE)
  ind_a2_mov <- runner::runner(ind_a2, k = l3, f = cust_fun, na_pad = TRUE)
  data$ind_a=as.factor(as.numeric(ind_a_mov+ind_a2_mov))
  
  #plot(data$a,col=ind_a_mov+1)
  
  # Plot 7: Scatter plot of a
  P7=ggplot(data, aes(x = t, y = a)) +
    geom_point(aes(color = type)) +
    scale_color_manual(values = custom_colors, name = "Type") +
    theme_classic()+
    labs(title = data_name)
  
  
  # Plot 8: Scatter plot of a with moving all(I(a<1))
  P8=ggplot(data, aes(x = t, y = a)) +
    geom_point(aes(color = ind_a)) +
    scale_color_manual(values = c("0" = "black", "1" = "magenta2"), name = "I(a<1|a>1)") +
    theme_classic()+
    labs(title = paste(data_name," lag = ", l3))
  
  # Return the dataset with the new features
  return(list(P_theta=P0,
              P_thetatrans=P1,
              P_count_min=P2,
              P_count_max=P3,
              P_value_min=P4,
              P_value_max=P5,
              P_diff_maxmin=P6,
              P_I_diff_maxmin=P6.1,
              P_a=P7,
              P_Ia=P8,
              data=data))
}

compute_feat_basic=function(data,wdn){
  # data is a dataframe with columns t, theta, a
  # wdn is a vector with the window sizes for the moving sd (two components)
  
  # Transform theta
  theta_trans <- data$theta
  theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
  data$theta_trans <- abs(theta_trans)
  data$theta=theta_trans
  
  # Add the first difference of theta and a
  data$dtheta=c(NA,diff(data$theta))
  data$da=c(NA,diff(data$a))
  
  # Moving sd
  data$sd_w1_dtheta=rollapply(data$dtheta, wdn[1], sd, fill=NA)
  data$sd_w2_dtheta=rollapply(data$dtheta, wdn[2], sd, fill=NA)
  
  data$sd_w1_da=rollapply(data$da, wdn[1], sd, fill=NA)
  data$sd_w2_da=rollapply(data$da, wdn[2], sd, fill=NA)
  
  data=data[complete.cases(data),]
  
  return(data)
  
}

comp_feat_theta=function(data,l2,l=5){
  
  # data is a dataframe with columns t and theta_trans
  
  # t_orig=data$t
  # data$t=seq_along(data$t)
  # data$t_orig=t_orig

  # Theta transformation
  maxs <- as.numeric(peaks(data$theta_trans, span = l))
  mins <- as.numeric(peaks(-data$theta_trans, span = l))
  
  count_min <- runner::runner(mins, k = l2, f = sum, na_pad = TRUE)
  count_min <- c(count_min[-(1:floor(l2 / 2))], rep(NA, round(l2 / 2)))
  
  count_max <- runner::runner(maxs, k = l2, f = sum, na_pad = TRUE)
  count_max <- c(count_max[-(1:floor(l2 / 2))], rep(NA, round(l2 / 2)))
  
  data$count_min <- count_min
  data$count_max <- count_max
  data$value_min <- last_min_value(mins, data$theta_trans)
  data$value_max <- last_max_value(maxs, data$theta_trans)
  data$diffmaxmin <- data$value_max - data$value_min
  
  data=subset(data,select=-c(theta_trans))
  
  #drop NA
  data=data[complete.cases(data),]
  
  colnames(data)[-1]=paste0(colnames(data)[-1],l2)
  
  return(data)
  
}

cust_fun=function(x){
  all(x==T)
}

comp_feat_a=function(data,l3){
  
  # data is a dataframe with columns t and a
  
  ind_a=I(data$a<1)
  ind_a2=I(data$a>1)
  ind_a_mov <- runner::runner(ind_a, k = l3, f = cust_fun, na_pad = TRUE)
  ind_a2_mov <- runner::runner(ind_a2, k = l3, f = cust_fun, na_pad = TRUE)
  data$ind_a=as.numeric(ind_a_mov+ind_a2_mov)
  
  data=subset(data,select=-a)
  
  colnames(data)[-1]=paste0(colnames(data)[-1],l3)
  data=data[complete.cases(data),]
  
  return(data)
}

comp_feat_theta_a=function(dat,wdn){
  dat=compute_feat_basic(dat[,c("t","a","theta")],wdn)
  dat_theta_short=comp_feat_theta(dat[,c("t","theta_trans")],wdn[1])
  dat_theta_long=comp_feat_theta(dat[,c("t","theta_trans")],wdn[2])
  #dat_a_short=comp_feat_a(dat[,c("t","a")],wdn[1])
  dat_a_long=comp_feat_a(dat[,c("t","a")],wdn[2])
  data=merge(dat,dat_theta_short,by="t")
  data=merge(data,dat_theta_long,by="t")
  #data=merge(data,dat_a_short,by="t")
  data=merge(data,dat_a_long,by="t")
  
  return(data)
  
}


feat_comput_theta_a <- function(data, data_name,
                                wdn=c(10,30),
                                l=5,
                                l2=c(10,30),
                                l3=c(10,30),
                                tt_thres_maxmin=2.7,
                                tt_thres_diffmaxmin=0.5) {
  
  # LAST UPDATE: 2024-10-22
  
  # Function to get the last observed min and max value, MA and SD of theta and dtheta, and moving indicator 
  # feature I(a<1|a>1)
  
  # Arguments:
  # data: dataframe with the data (t theta and type columns are required)
  # data_name: name of the dataset
  # wdn: window sizes for moving average and standard deviation of theta and dtheta (default 10,30)
  # l: span for peaks function (default 5)
  # l2: window sizes for runner function (default 10 and 30)
  # l3: window sizes for I(a<1|a>1) (default 10 and 30)
  # tt_thres: threshold for theta (peaks (and valleys) above (below) this value are not considered as peaks (valleys))
  # tt_thres_diffmaxmin: threshold for the difference between max and min
  
  # t_orig=data$t
  # data$t=seq_along(data$t)
  # data$t_orig=t_orig
  #data$type <- as.factor(data$type)
  
  # Unique types and total rows
  # unique_types <- unique(data$type)
  # total_rows <- dim(data)[1]
  # 
  # # Custom colors based on unique 'type'
  # custom_colors <- c("-1" = "blue", "0" = "red", "50" = "green", "100" = "gray50")
  
  # Theta transformation
  theta_trans <- data$theta
  theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
  data$theta <- theta_trans
  
  # # Plot 0: Scatter plot of theta
  # P0=ggplot(data, aes(x = t, y = theta)) +
  #   geom_point(aes(color = type)) +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_classic()+
  #   labs(title = data_name)
  # 
  # # Plot 1: Scatter plot of theta_trans
  # P1=ggplot(data, aes(x = t, y = theta_trans)) +
  #   geom_point(aes(color = type)) +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_classic()+
  #   labs(title = data_name)
  
  # Detecting peaks and valleys
  #l <- 5
  maxs <- as.numeric(peaks(data$theta, span = l))
  mins <- as.numeric(peaks(-data$theta, span = l))
  
  # Values above tt_thres_maxmin and below -tt_thres_maxmin are not considered as max and min
  maxs[which(data$theta>(tt_thres_maxmin))] <- 0
  mins[which(data$theta<(-tt_thres_maxmin))] <- 0
  
  # Smoothing with runner function
  #l2 <- 50
  count_min_short <- runner::runner(mins, k = l2[1], f = sum, na_pad = TRUE)
  count_min_short <- c(count_min_short[-(1:floor(l2[1] / 2))], rep(NA, round(l2[1] / 2)))
  
  count_min_long <- runner::runner(mins, k = l2[2], f = sum, na_pad = TRUE)
  count_min_long <- c(count_min_long[-(1:floor(l2[2] / 2))], rep(NA, round(l2[2] / 2)))
  
  count_max_short <- runner::runner(maxs, k = l2[1], f = sum, na_pad = TRUE)
  count_max_short <- c(count_max_short[-(1:floor(l2[1] / 2))], rep(NA, round(l2[1] / 2)))
  
  count_max_long <- runner::runner(maxs, k = l2[2], f = sum, na_pad = TRUE)
  count_max_long <- c(count_max_long[-(1:floor(l2[2] / 2))], rep(NA, round(l2[2] / 2)))
  
  
  data$count_min_short <- count_min_short
  data$count_max_short <- count_max_short
  
  data$count_min_long <- count_min_long
  data$count_max_long <- count_max_long
  
  # Plot 2: Scatter plot of count_min
  # P2=ggplot(data, aes(x = t, y = count_min, color = type)) +
  #   geom_point() +
  #   labs(x = "t", y = "count_min") +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_minimal() +
  #   theme(legend.position = "top")+
  #   labs(title = paste(data_name," lag = ", l2))
  # 
  # # Plot 3: Scatter plot of count_max
  # P3=ggplot(data, aes(x = t, y = count_max, color = type)) +
  #   geom_point() +
  #   labs(x = "t", y = "count_max") +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_minimal() +
  #   theme(legend.position = "top")+
  #   labs(title = paste(data_name," lag = ", l2))
  
  # Add value_min, value_max and diffmaxmin
  data$value_min <- last_min_value(mins, data$theta)
  data$value_max <- last_max_value(maxs, data$theta)
  data$diffmaxmin <- data$value_max - data$value_min
  
  data$I_diffmaxmin = as.numeric(data$diffmaxmin>tt_thres_diffmaxmin)
  
  data$dtheta=c(NA,diff(data$theta))
  
  data$ma_w1_theta=rollapply(data$theta, wdn[1], mean, fill=NA)
  data$sd_w1_theta=rollapply(data$theta, wdn[1], sd, fill=NA)
  data$sd_w1_dtheta=rollapply(data$dtheta, wdn[1], sd, fill=NA)
  
  if(length(wdn)>1){
    data$ma_w2_theta=rollapply(data$theta, wdn[2], mean, fill=NA)
    data$sd_w2_theta=rollapply(data$theta, wdn[2], sd, fill=NA)
    data$sd_w2_dtheta=rollapply(data$dtheta, wdn[2], sd, fill=NA)
  }
  
  # # Plot 4: Scatter plot of value_min
  # P4=ggplot(data, aes(x = t, y = value_min, color = type)) +
  #   geom_point() +
  #   labs(x = "t", y = "value_min") +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_minimal() +
  #   theme(legend.position = "top")+
  #   labs(title = paste(data_name," lag = ", l2))
  # 
  # # Plot 5: Scatter plot of value_max
  # P5=ggplot(data, aes(x = t, y = value_max, color = type)) +
  #   geom_point() +
  #   labs(x = "t", y = "value_max") +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_minimal() +
  #   theme(legend.position = "top")+
  #   labs(title = paste(data_name," lag = ", l2))
  # 
  # # Plot 6: Scatter plot of diffmaxmin
  # P6=ggplot(data, aes(x = t, y = diffmaxmin, color = type)) +
  #   geom_point() +
  #   labs(x = "t", y = "diffmaxmin") +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_minimal() +
  #   theme(legend.position = "top")+
  #   labs(title = paste(data_name," lag = ", l2))
  # 
  # P6.1=ggplot(data, aes(x = t, y = I_diffmaxmin, color = type)) +
  #   geom_point() +
  #   labs(x = "t", y = paste("I(diffmaxmin>", tt_thres_diffmaxmin, ")")) +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_minimal() +
  #   theme(legend.position = "top")+
  #   labs(title = paste(data_name," lag = ", l2))
  
  
  cust_fun=function(x){
    all(x==T)
  }
  
  ind_a=I(data$a<1)
  ind_a2=I(data$a>1)
  
  # Short
  ind_a_mov <- runner::runner(ind_a, k = l3[1], f = cust_fun, na_pad = TRUE)
  ind_a2_mov <- runner::runner(ind_a2, k = l3[1], f = cust_fun, na_pad = TRUE)
  #data$ind_a_short=as.factor(as.numeric(ind_a_mov+ind_a2_mov))
  data$ind_a_short=as.numeric(ind_a_mov+ind_a2_mov)
  
  # Long
  ind_a_mov <- runner::runner(ind_a, k = l3[2], f = cust_fun, na_pad = TRUE)
  ind_a2_mov <- runner::runner(ind_a2, k = l3[2], f = cust_fun, na_pad = TRUE)
  #data$ind_a_long=as.factor(as.numeric(ind_a_mov+ind_a2_mov))
  data$ind_a_long=as.numeric(ind_a_mov+ind_a2_mov)
  
  # TD indicator
  data$I_TD=as.numeric(data$value_max*data$value_min<0)
  
  #plot(data$a,col=ind_a_mov+1)
  
  # Plot 7: Scatter plot of a
  # P7=ggplot(data, aes(x = t, y = a)) +
  #   geom_point(aes(color = type)) +
  #   scale_color_manual(values = custom_colors, name = "Type") +
  #   theme_classic()+
  #   labs(title = data_name)
  # 
  # 
  # # Plot 8: Scatter plot of a with moving all(I(a<1))
  # P8=ggplot(data, aes(x = t, y = a)) +
  #   geom_point(aes(color = ind_a)) +
  #   scale_color_manual(values = c("0" = "black", "1" = "magenta2"), name = "I(a<1|a>1)") +
  #   theme_classic()+
  #   labs(title = paste(data_name," lag = ", l3))
  
  # Return the dataset with the new features
  return(
    # list(P_theta=P0,
    #           P_thetatrans=P1,
    #           P_count_min=P2,
    #           P_count_max=P3,
    #           P_value_min=P4,
    #           P_value_max=P5,
    #           P_diff_maxmin=P6,
    #           P_I_diff_maxmin=P6.1,
    #           P_a=P7,
    #           P_Ia=P8,
    data=data
    #)
  )
}
