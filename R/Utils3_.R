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

SJM_lambdakappa=function(lambda,kappa,K=2,df,Lnsat,Ksat=6,alpha0=NULL,K0=NULL,pers0=.95,true_states=NULL){
  
  #  df is a data.frame WITHOUT time column
  
  if(is.null(K0)){
    K0=K
  }
  
  # library(reticulate)
  # import("scipy")
  # source_python('SJ.py')
  
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
    true_states=order_states(true_states)
    est_states=order_states(est_states)
    overlap=sum(true_states==est_states)/N
    ARI=pdfCluster::adj.rand.index(true_states,est_states)
    return(list(FTIC=FTIC,
                BIC=BIC,
                AIC=AIC,
                ARI=ARI,
                overlap=overlap,
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
compute_feat=function(dat,wdn=c(10,75),am1=F,wdn_decomp=10){
  
  # Add first differences for each variable to dat2
  dat$da=c(NA,diff(dat$a))
  dat$de=c(NA,diff(dat$e))
  dat$dtheta=c(NA,diff(dat$theta))
  dat$domega=c(NA,diff(dat$omega))
  
  ###
  # Add moving average
  
  if(length(wdn)>1){
    dat$ma_long_a=rollapply(dat$a, wdn[2], mean, fill=NA)
    dat$ma_long_e=rollapply(dat$e, wdn[2], mean, fill=NA)
    dat$ma_long_theta=rollapply(dat$theta, wdn[2], mean, fill=NA)
    dat$ma_long_omega=rollapply(dat$omega, wdn[2], mean, fill=NA)
    
    dat$sd_long_a=rollapply(dat$a, wdn[2], sd, fill=NA)
    dat$sd_long_e=rollapply(dat$e, wdn[2], sd, fill=NA)
    dat$sd_long_theta=rollapply(dat$theta, wdn[2], sd, fill=NA)
    dat$sd_long_omega=rollapply(dat$omega, wdn[2], sd, fill=NA)
    
    dat$sd_long_da=rollapply(dat$da, wdn[2], sd, fill=NA)
    dat$sd_long_de=rollapply(dat$de, wdn[2], sd, fill=NA)
    dat$sd_long_dtheta=rollapply(dat$dtheta, wdn[2], sd, fill=NA)
    dat$sd_long_domega=rollapply(dat$domega, wdn[2], sd, fill=NA)
    
    dat$corr_long_da_de=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","de")],
                                                     width=wdn[2], function(x) cor(x[,1],x[,2]),
                                                     by.column=FALSE))
    dat$corr_long_da_dtheta=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","dtheta")],
                                                         width=wdn[2], function(x) cor(x[,1],x[,2]),
                                                         by.column=FALSE))
    dat$corr_long_da_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("da","domega")],
                                                         width=wdn[2], function(x) cor(x[,1],x[,2]),
                                                         by.column=FALSE))
    dat$corr_long_de_dtheta=c(rep(NA,wdn[2]-1),rollapply(dat[,c("de","dtheta")],
                                                         width=wdn[2], function(x) cor(x[,1],x[,2]), 
                                                         by.column=FALSE))
    dat$corr_long_de_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("de","domega")],
                                                         width=wdn[2], function(x) cor(x[,1],x[,2]), 
                                                         by.column=FALSE))
    dat$corr_long_dtheta_domega=c(rep(NA,wdn[2]-1),rollapply(dat[,c("dtheta","domega")],
                                                             width=wdn[2], function(x) cor(x[,1],x[,2]), 
                                                             by.column=FALSE))
  }
  
  dat$ma_short_a=rollapply(dat$a, wdn[1], mean, fill=NA)
  dat$ma_short_e=rollapply(dat$e, wdn[1], mean, fill=NA)
  dat$ma_short_theta=rollapply(dat$theta, wdn[1], mean, fill=NA)
  dat$ma_short_omega=rollapply(dat$omega, wdn[1], mean, fill=NA)
  
  # Add moving standard deviation for each variable to dat
  dat$sd_short_a=rollapply(dat$a, wdn[1], sd, fill=NA)
  dat$sd_short_e=rollapply(dat$e, wdn[1], sd, fill=NA)
  dat$sd_short_theta=rollapply(dat$theta, wdn[1], sd, fill=NA)
  dat$sd_short_omega=rollapply(dat$omega, wdn[1], sd, fill=NA)
  
  # Add moving st. dev. of first differences
  
  # Add moving correlation 
  library(zoo)
  
  # Add moving correlations between first differences
  dat$corr_short_da_de=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","de")],
                                                    width=wdn[1], function(x) cor(x[,1],x[,2]),
                                                    by.column=FALSE))
  dat$corr_short_da_dtheta=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","dtheta")],
                                                        width=wdn[1], function(x) cor(x[,1],x[,2]),
                                                        by.column=FALSE))
  dat$corr_short_da_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("da","domega")],
                                                        width=wdn[1], function(x) cor(x[,1],x[,2]),
                                                        by.column=FALSE))
  dat$corr_short_de_dtheta=c(rep(NA,wdn[1]-1),rollapply(dat[,c("de","dtheta")],
                                                        width=wdn[1], function(x) cor(x[,1],x[,2]), 
                                                        by.column=FALSE))
  dat$corr_short_de_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("de","domega")],
                                                        width=wdn[1], function(x) cor(x[,1],x[,2]), 
                                                        by.column=FALSE))
  dat$corr_short_dtheta_domega=c(rep(NA,wdn[1]-1),rollapply(dat[,c("dtheta","domega")],
                                                            width=wdn[1], function(x) cor(x[,1],x[,2]), 
                                                            by.column=FALSE))
  
  if(am1){
    dat$am1=dat$a-1
    dat$belowabove1=as.numeric(
      rollapply(dat$a, wdn[2], function(x){all(x<1)|all(x>1)}, fill=NA)
    )
  }
  
  tsa=ts(dat$a,frequency = wdn_decomp)
  tsa.stl <- stl(tsa, 
                 s.window=wdn_decomp,
                 t.window=NULL,
                 na.action = na.approx,
                 robust=T)
  dat$seas_a=tsa.stl$time.series[,1]
  dat$trend_a=tsa.stl$time.series[,2]
  dat$remainder_a=tsa.stl$time.series[,3]
  
  tse=ts(dat$e,frequency = wdn_decomp)
  tse.stl <- stl(tse, 
                 s.window=wdn_decomp,
                 t.window=NULL,
                 na.action = na.approx,
                 robust=T)
  
  dat$seas_e=tsa.stl$time.series[,1]
  dat$trend_e=tsa.stl$time.series[,2]
  dat$remainder_e=tsa.stl$time.series[,3]

  tstheta=ts(dat$theta,frequency = wdn_decomp)
  tstheta.stl <- stl(tstheta, 
                     s.window=wdn_decomp,
                     t.window=NULL,
                     na.action = na.approx,
                     robust=T)
  
  dat$seas_theta=tsa.stl$time.series[,1]
  dat$trend_theta=tsa.stl$time.series[,2]
  dat$remainder_theta=tsa.stl$time.series[,3]
  
  tsomega=ts(dat$omega,frequency = wdn_decomp)
  tsomega.stl <- stl(tsomega, 
                     s.window=wdn_decomp,
                     t.window=NULL,
                     na.action = na.approx,
                     robust=T)
  
  dat$seas_omega=tsa.stl$time.series[,1]
  dat$trend_omega=tsa.stl$time.series[,2]
  dat$remainder_omega=tsa.stl$time.series[,3]
  
  # Add moving standard deviation for trend, seas, and reminder of each variable
  dat$sd_short_seas_a=rollapply(dat$seas_a, wdn[1], sd, fill=NA)
  dat$sd_short_trend_a=rollapply(dat$trend_a, wdn[1], sd, fill=NA)
  dat$sd_short_remainder_a=rollapply(dat$remainder_a, wdn[1], sd, fill=NA)
  
  dat$sd_short_seas_e=rollapply(dat$seas_e, wdn[1], sd, fill=NA)
  dat$sd_short_trend_e=rollapply(dat$trend_e, wdn[1], sd, fill=NA)
  dat$sd_short_remainder_e=rollapply(dat$remainder_e, wdn[1], sd, fill=NA)
  
  dat$sd_short_seas_theta=rollapply(dat$seas_theta, wdn[1], sd, fill=NA)
  dat$sd_short_trend_theta=rollapply(dat$trend_theta, wdn[1], sd, fill=NA)
  dat$sd_short_remainder_theta=rollapply(dat$remainder_theta, wdn[1], sd, fill=NA)

  dat$sd_short_seas_omega=rollapply(dat$seas_omega, wdn[1], sd, fill=NA)
  dat$sd_short_trend_omega=rollapply(dat$trend_omega, wdn[1], sd, fill=NA)
  dat$sd_short_remainder_omega=rollapply(dat$remainder_omega, wdn[1], sd, fill=NA)
  
  if(length(wdn)>1){
    dat$sd_long_seas_a=rollapply(dat$seas_a, wdn[2], sd, fill=NA)
    dat$sd_long_trend_a=rollapply(dat$trend_a, wdn[2], sd, fill=NA)
    dat$sd_long_remainder_a=rollapply(dat$remainder_a, wdn[2], sd, fill=NA)
    
    dat$sd_long_seas_e=rollapply(dat$seas_e, wdn[2], sd, fill=NA)
    dat$sd_long_trend_e=rollapply(dat$trend_e, wdn[2], sd, fill=NA)
    dat$sd_long_remainder_e=rollapply(dat$remainder_e, wdn[2], sd, fill=NA)
    
    dat$sd_long_seas_theta=rollapply(dat$seas_theta, wdn[2], sd, fill=NA)
    dat$sd_long_trend_theta=rollapply(dat$trend_theta, wdn[2], sd, fill=NA)
    dat$sd_long_remainder_theta=rollapply(dat$remainder_theta, wdn[2], sd, fill=NA)
    
    dat$sd_long_seas_omega=rollapply(dat$seas_omega, wdn[2], sd, fill=NA)
    dat$sd_long_trend_omega=rollapply(dat$trend_omega, wdn[2], sd, fill=NA)
    dat$sd_long_remainder_omega=rollapply(dat$remainder_omega, wdn[2], sd, fill=NA)
  }
  
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

