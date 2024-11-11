library(ggplot2)
library(zoo)
library(splus2R)
library(dplyr)

trans_theta=function(data){
  theta_trans <- data$theta
  theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
  data$theta <- theta_trans
  return(data)
}

inv_trans_theta=function(data){
  theta_inv=data$theta
  theta_inv[which(theta_inv<0)]=theta_inv[which(theta_inv<0)]+2*pi
  data$theta <- theta_inv
  return(data)
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


# Simple cases ------------------------------------------------------------

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
  
  #datC=dat[complete.cases(dat),]
  
  return(dat)
  
}



# Complex cases -----------------------------------------------------------

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

max_min_feat=function(data,tt_thres_maxmin=3,
                      tt_thres_diffmaxmin=pi/4,
                      l=5
                      ,l2=c(5,30)
){
  
  # Theta transformation
  # theta_trans <- data$theta
  # theta_trans[which(theta_trans > pi)] <- theta_trans[which(theta_trans > pi)] - 2 * pi
  # data$theta <- theta_trans
  
  maxs <- as.numeric(peaks(data$theta, span = l))
  mins <- as.numeric(peaks(-data$theta, span = l))
  
  # Values above tt_thres_maxmin and below -tt_thres_maxmin are not considered as max and min
  maxs[which(data$theta>(tt_thres_maxmin))] <- 0
  mins[which(data$theta<(-tt_thres_maxmin))] <- 0
  
  #
  data$maxs=maxs
  data$mins=mins
  
  # Add value_min, value_max and diffmaxmin
  data$value_min <- last_min_value(mins, data$theta)
  data$value_max <- last_max_value(maxs, data$theta)
  
  data$value_min <- last_min_value(mins, data$theta)
  data$value_max <- last_max_value(maxs, data$theta)
  data$diffmaxmin <- data$value_max - data$value_min
  data$I_diffmaxmin = as.numeric(data$diffmaxmin<tt_thres_diffmaxmin)
  
  count_min_short <- runner::runner(mins, k = l2[1], f = sum, na_pad = TRUE)
  count_min_short <- c(count_min_short[-(1:floor(l2[1] / 2))], rep(NA, round(l2[1] / 2)))
  count_max_short <- runner::runner(maxs, k = l2[1], f = sum, na_pad = TRUE)
  count_max_short <- c(count_max_short[-(1:floor(l2[1] / 2))], rep(NA, round(l2[1] / 2)))
  data$count_min_short <- count_min_short
  data$count_max_short <- count_max_short
  
  if(length(l2)>1){
    count_min_long <- runner::runner(mins, k = l2[2], f = sum, na_pad = TRUE)
    count_min_long <- c(count_min_long[-(1:floor(l2[2] / 2))], rep(NA, round(l2[2] / 2)))
    
    count_max_long <- runner::runner(maxs, k = l2[2], f = sum, na_pad = TRUE)
    count_max_long <- c(count_max_long[-(1:floor(l2[2] / 2))], rep(NA, round(l2[2] / 2)))
    
    data$count_min_long <- count_min_long
    data$count_max_long <- count_max_long
  }
  
  
  # TD indicator
  data$I_TD=as.numeric(data$value_max*data$value_min>0)
  
  # HS indicator
  data$I_HS=as.numeric(data$value_max*data$value_min<0&data$value_max<0)
  
  # QS indicator
  data$I_QS=as.numeric(data$value_max*data$value_min<0&data$value_max>0)
  
  # Check if distance between max and min is small to distinguish between NR and QS
  # data$I_QS=data$I_QS*data$I_diffmaxmin
  
  
  
  data$mean_osc=data$I_HS*(data$value_min+data$value_max+2*pi)/2+
    data$I_QS*(data$value_min+data$value_max)/2+
    data$I_TD*(data$value_min+data$value_max)/2
  
  data=data[,c("t","maxs","mins",
               "value_min","value_max",
               "diffmaxmin",
               "I_diffmaxmin",
               "I_TD","I_HS","I_QS",
               "mean_osc")]
  
  return(data)
}

a_feat=function(data,l3){
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
  if(length(l3)>1){
    ind_a_mov <- runner::runner(ind_a, k = l3[2], f = cust_fun, na_pad = TRUE)
    ind_a2_mov <- runner::runner(ind_a2, k = l3[2], f = cust_fun, na_pad = TRUE)
    #data$ind_a_long=as.factor(as.numeric(ind_a_mov+ind_a2_mov))
    data$ind_a_long=as.numeric(ind_a_mov+ind_a2_mov)
  }
  
  return(data[,c("t","ind_a_short","ind_a_long")])
  
}

thetavol_feat=function(data,wdn=c(5,30)){
  
  data$dtheta=c(NA,diff(data$theta))
  
  data$ma_theta_short=rollapply(data$theta, wdn[1], mean, fill=NA)
  data$sd_theta_short=rollapply(data$theta, wdn[1], sd, fill=NA)
  data$sd_dtheta_short=rollapply(data$dtheta, wdn[1], sd, fill=NA)
  
  if(length(wdn)>1){
    data$ma_theta_long=rollapply(data$theta, wdn[2], mean, fill=NA)
    data$sd_theta_long=rollapply(data$theta, wdn[2], sd, fill=NA)
    data$sd_dtheta_long=rollapply(data$dtheta, wdn[2], sd, fill=NA)
  }
  
  return(data[,c("t","ma_theta_short","sd_theta_short",
                 "ma_theta_long","sd_theta_long"
                 ,"sd_dtheta_long","sd_dtheta_short"
  )])
}

BAC=function(obj1,obj2,levs=3){
  A=as.character(obj1)
  B=as.character(obj2)
  
  A=factor(A,levels=1:levs)
  B=factor(B,levels=1:levs)
  
  return(caret::confusionMatrix(A,B)$overall[1])
}
