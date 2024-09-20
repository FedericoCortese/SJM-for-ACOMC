source("Utils3_.R")
df2015XX169=read.table("propagation_2015XX169_new_v2.txt",header=T)

dat=df2015XX169
dat$type=as.factor(dat$type)
TT=dim(dat)[1]
custom_colors <- c("-1" = "blue", "0" = "red", "50" = "green")

# a --------------------------------------------------------------

windows()
ggplot(dat, aes(x = t, y = a)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()

# da interesting only in terms of volatility
dat$da=c(NA,diff(dat$a))
windows()
ggplot(dat, aes(x = t, y = da)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()

dat$d2a=c(NA,diff(dat$da))
windows()
ggplot(dat, aes(x = t, y = d2a)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()

# DD2a=DD(dat$a,2,out.plot = T)
# windows()
# DD2a$plot
# 
# res_dfa=last_min_max(last_min=DD2a$tmin,last_max=DD2a$tmax,TT=TT,dat$a)
# windows()
# plot(res_dfa$last_observed_max-res_dfa$last_observed_min,col=dat$type,
#      ylab='last_observed_max-last_observed_min',xlab='t',main='a')


# Local min max at lag 10
l=5
times=seq(1,TT)
library(splus2R)
maxs=as.numeric(peaks(dat$a,span=l))
mins=as.numeric(peaks(-dat$a,span=l))

res_dfa_short=last_min_max(last_min=times[which(mins==1)],last_max=times[which(maxs==1)],TT=TT,dat$a)
res_dfa_short$diff <- (res_dfa_short$last_observed_max - res_dfa_short$last_observed_min)/2

custom_colors <- c("-1" = "blue", "0" = "red", "50" = "green")
windows()
ggplot(res_dfa_short, aes(x = 1:nrow(res_dfa_short), y = abs(diff), color = factor(dat$type))) +
  geom_point() +
  labs(x = "t", y = "|M(a) - m(a)|", title = paste("lag =", l)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_color_manual(values = custom_colors, name = "Type") +  # Apply custom colors
  guides(color = guide_legend(override.aes = list(shape = 19)))

# Local min max at 20
l=20
maxs=as.numeric(peaks(dat$a,span=l))
mins=as.numeric(peaks(-dat$a,span=l))
res_dfa_long=last_min_max(last_min=times[which(mins==1)],last_max=times[which(maxs==1)],TT=TT,dat$a)
res_dfa_long$diff <- (res_dfa_long$last_observed_max - res_dfa_long$last_observed_min)/2

windows()
ggplot(res_dfa_long, aes(x = 1:nrow(res_dfa_long), y = abs(diff), color = factor(dat$type))) +
  geom_point() +
  labs(x = "t", y = "|M(a) - m(a)|", title = paste("lag =", l)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_color_manual(values = custom_colors, name = "Type") +  # Apply custom colors
  guides(color = guide_legend(override.aes = list(shape = 19)))




# theta -------------------------------------------------------------------
windows()
ggplot(dat, aes(x = t, y = theta)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()

dat$dtheta=c(NA,diff(dat$theta))
windows()
ggplot(dat, aes(x = t, y = dtheta)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()

d=5
dat$ddtheta=c(rep(NA,d),diff(dat$theta,differences = d))
windows()
ggplot(dat, aes(x = t, y = ddtheta)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()+
  labs(title=paste("d=",d))

# DD2theta=DD(dat$theta,2,out.plot = T)
# windows()
# DD2theta$plot
# last_min=DD2theta$tmin
# last_max=DD2theta$tmax
# 
# res_dftheta=last_min_max(last_min=DD2theta$tmin,last_max=DD2theta$tmax,TT=TT,dat$theta)
# 
# plot(res_dftheta$last_observed_max-res_dftheta$last_observed_min,col=dat$type,
#      ylab='last_observed_max-last_observed_min',xlab='t',main='theta')

l = 3
maxs = as.numeric(peaks(dat$theta, span = l))   # Use theta instead of a
mins = as.numeric(peaks(-dat$theta, span = l))  # Use theta instead of a

# Compute last observed min and max for short span (lag=10)
res_dftheta_short = last_min_max(last_min = times[which(mins == 1)], last_max = times[which(maxs == 1)], TT = TT, dat$theta)
res_dftheta_short$diff <- (res_dftheta_short$last_observed_max - res_dftheta_short$last_observed_min)/2

# Create the plot for lag=10
windows()
ggplot(res_dftheta_short, aes(x = 1:nrow(res_dftheta_short), y = diff, color = factor(dat$type))) +
  geom_point() +
  labs(x = "t", y = "(M(theta) - m(theta))/2", title = paste("lag =", l)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_color_manual(values = custom_colors, name = "Type") +  # Apply custom colors
  guides(color = guide_legend(override.aes = list(shape = 19)))  # Matching the pch 19

# Local min max at lag=20 for theta
l = 20
maxs = as.numeric(peaks(dat$theta, span = l))   # Use theta instead of a
mins = as.numeric(peaks(-dat$theta, span = l))  # Use theta instead of a

# Compute last observed min and max for long span (lag=20)
res_dftheta_long = last_min_max(last_min = times[which(mins == 1)], last_max = times[which(maxs == 1)], TT = TT, dat$theta)
res_dftheta_long$diff <- res_dftheta_long$last_observed_max - res_dftheta_long$last_observed_min

# Create the plot for lag=20
windows()
ggplot(res_dftheta_long, aes(x = 1:nrow(res_dftheta_long), y = abs(diff), color = factor(dat$type))) +
  geom_point() +
  labs(x = "t", y = "|M(theta) - m(theta)|", title = paste("lag =", l)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_color_manual(values = custom_colors, name = "Type") +  # Apply custom colors
  guides(color = guide_legend(override.aes = list(shape = 19)))  # Matching the pch 19


# theta transformed -------------------------------------------------------

theta_trans=dat$theta
theta_trans[which(theta_trans>pi)]=theta_trans[which(theta_trans>pi)]-2*pi
dat$theta_trans=abs(theta_trans)


G1=ggplot(dat[1:250,], aes(x = t, y = theta_trans)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()
G2=ggplot(dat[1700:1999,], aes(x = t, y = theta_trans)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()
windows()
ggpubr::ggarrange(G1,G2,nrow=2)

dat$dtheta_trans=c(NA,diff(dat$theta_trans))
windows()
ggplot(dat, aes(x = t, y = dtheta_trans)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()


l = 5
times=seq(1,TT)
library(splus2R)
maxs = as.numeric(peaks(dat$theta_trans, span = l))   
mins = as.numeric(peaks(-dat$theta_trans, span = l)) 

plot(dat$theta_trans,col=maxs+1,pch=19)
plot(dat$theta_trans,col=mins+1,pch=19)

l2=25
#
count_min=runner::runner(
  mins,
  k = l2,
  f = mean,
  na_pad=T
)
count_min=count_min[-(1:floor(l2/2))]
count_min=c(count_min,rep(NA,round(l2/2)))
plot(count_min,col=dat$type,ylab='count_min',xlab='t')
legend('topright',legend=levels(dat$type),col=1:3,pch=19)
# 
count_max=runner::runner(
  maxs,
  k = l2,
  f = mean,
  na_pad=T
)
count_max=count_max[-(1:floor(l2/2))]
count_max=c(count_max,rep(NA,round(l2/2)))
plot(count_max,col=dat$type,ylab='count_max',xlab='t')


# Compute last observed min and max for short span (lag=10)
# res_dftheta_trans = last_min_max(last_min = times[which(mins == 1)], 
#                                  last_max = times[which(maxs == 1)], TT = TT, dat$theta_trans)
# res_dftheta_trans$diff <- (res_dftheta_short$last_observed_max - res_dftheta_short$last_observed_min)/2
# 
# plot(res_dftheta_trans$diff,col=dat$type,ylab='(M(theta) - m(theta))/2',xlab='t')

# Count min and max in the neighborhood
# library(zoo)
# count_min=rollapply(mins, width = l2, sum, align = "left")
# count_max=rollapply(maxs, width = l2, sum, align = "center")


# e -----------------------------------------------------------------------

windows()
ggplot(dat, aes(x = t, y = e)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()

# omega -------------------------------------------------------------------

windows()
ggplot(dat, aes(x = t, y = omega)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()

# DD1omega=DD(dat$omega,1,out.plot = T)
# windows()
# DD1omega$plot

times = seq(1, TT)
library(splus2R)

# Local min max at lag=2 for omega
l = 5
maxs = as.numeric(peaks(dat$omega, span = l))   # Use omega instead of theta
mins = as.numeric(peaks(-dat$omega, span = l))  # Use omega instead of theta

# Compute last observed min and max for short span (lag=2)
res_dfomega_short = last_min_max(last_min = times[which(mins == 1)], last_max = times[which(maxs == 1)], 
                                 TT = TT, dat$omega)
res_dfomega_short$diff <- res_dfomega_short$last_observed_max - res_dfomega_short$last_observed_min

# Create the plot for lag=2
windows()
ggplot(res_dfomega_short, aes(x = 1:nrow(res_dfomega_short), y = diff, color = factor(dat$type))) +
  geom_point() +
  labs(x = "t", y = "M(omega) - m(omega)", title = paste("lag =", l)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_color_manual(values = custom_colors, name = "Type") +  # Apply custom colors
  guides(color = guide_legend(override.aes = list(shape = 19)))  # Matching the pch 19

# Local min max at lag=20 for omega
l = 20
maxs = as.numeric(peaks(dat$omega, span = l))   # Use omega instead of theta
mins = as.numeric(peaks(-dat$omega, span = l))  # Use omega instead of theta

# Compute last observed min and max for long span (lag=20)
res_dfomega_long = last_min_max(last_min = times[which(mins == 1)], last_max = times[which(maxs == 1)], TT = TT, dat$omega)
res_dfomega_long$diff <- res_dfomega_long$last_observed_max - res_dfomega_long$last_observed_min

# Create the plot for lag=20
windows()
ggplot(res_dfomega_long, aes(x = 1:nrow(res_dfomega_long), y = abs(diff), color = factor(dat$type))) +
  geom_point() +
  labs(x = "t", y = "|M(omega) - m(omega)|", title = paste("lag =", l)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_color_manual(values = custom_colors, name = "Type") +  # Apply custom colors
  guides(color = guide_legend(override.aes = list(shape = 19)))  # Matching the pch 19


l=5
precgreat=as.numeric(I(dat$omega[1:(TT-(l-1))]>dat$omega[l:TT]))

plot(precgreat,col=dat$type,ylab='precgreat',xlab='t',main='omega')


