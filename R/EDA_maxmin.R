source("Utils3_.R")
df2015XX169=read.table("propagation_2015XX169_new_v2.txt",header=T)

dat=df2015XX169
dat$type=as.factor(dat$type)
TT=dim(dat)[1]

# a features --------------------------------------------------------------
windows()
ggplot(dat, aes(x = t, y = a)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = c("blue", "red","green")) +
  theme_classic()

dat$da=c(NA,diff(dat$a))
windows()
ggplot(dat, aes(x = t, y = da)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = c("blue", "red","green")) +
  theme_classic()

dat$d2a=c(NA,diff(dat$da))
windows()
ggplot(dat, aes(x = t, y = d2a)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = c("blue", "red","green")) +
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
l=10
times=seq(1,TT)
library(splus2R)
maxs=as.numeric(peaks(dat$a,span=l))
mins=as.numeric(peaks(-dat$a,span=l))

res_dfa_short=last_min_max(last_min=times[which(mins==1)],last_max=times[which(maxs==1)],TT=TT,dat$a)
res_dfa_short$diff <- (res_dfa_short$last_observed_max - res_dfa_short$last_observed_min)/2

custom_colors <- c("-1" = "black", "0" = "blue", "50" = "red")
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

custom_colors <- c("-1" = "black", "0" = "blue", "50" = "red")
windows()
ggplot(res_dfa_long, aes(x = 1:nrow(res_dfa_long), y = abs(diff), color = factor(dat$type))) +
  geom_point() +
  labs(x = "t", y = "|M(a) - m(a)|", title = paste("lag =", l)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
  scale_color_manual(values = custom_colors, name = "Type") +  # Apply custom colors
  guides(color = guide_legend(override.aes = list(shape = 19)))




# theta features -------------------------------------------------------------------
windows()
ggplot(dat, aes(x = t, y = theta)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = c("blue", "red","green")) +
  theme_classic()

dat$dtheta=c(NA,diff(dat$theta))
windows()
ggplot(dat, aes(x = t, y = dtheta)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = c("blue", "red","green")) +
  theme_classic()

d=5
dat$ddtheta=c(rep(NA,d),diff(dat$theta,differences = d))
windows()
ggplot(dat, aes(x = t, y = ddtheta)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = c("blue", "red","green")) +
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

l = 2
maxs = as.numeric(peaks(dat$theta, span = l))   # Use theta instead of a
mins = as.numeric(peaks(-dat$theta, span = l))  # Use theta instead of a

# Compute last observed min and max for short span (lag=10)
res_dftheta_short = last_min_max(last_min = times[which(mins == 1)], last_max = times[which(maxs == 1)], TT = TT, dat$theta)
res_dftheta_short$diff <- (res_dftheta_short$last_observed_max - res_dftheta_short$last_observed_min)/2

# Define custom colors
custom_colors <- c("-1" = "black", "0" = "blue", "50" = "red")

# Create the plot for lag=10
windows()
ggplot(res_dftheta_short, aes(x = 1:nrow(res_dftheta_short), y = diff, color = factor(dat$type))) +
  geom_point() +
  labs(x = "t", y = "M(theta) - m(theta)", title = paste("lag =", l)) +
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


# e -----------------------------------------------------------------------

windows()
ggplot(dat, aes(x = t, y = e)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = c("blue", "red","green")) +
  theme_classic()

# omega -------------------------------------------------------------------

windows()
ggplot(dat, aes(x = t, y = omega)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = c("blue", "red","green")) +
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

# Define custom colors
custom_colors <- c("-1" = "black", "0" = "blue", "50" = "red")

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


precgreat=as.numeric(I(dat$omega[1:(TT-1)]>dat$omega[2:TT]))

preclow=as.numeric(I(dat$omega[1:(TT-1)]<dat$omega[2:TT]))
