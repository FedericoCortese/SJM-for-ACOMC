source("Utils3_.R")
df164207=read.table("propagation_164207_new_v2.txt",header=T)

dat=df164207
dat$type=as.factor(dat$type)
ut=unique(dat$type)
TT=dim(dat)[1]
#custom_colors <- c("-1" = "blue", "0" = "red", "50" = "green")
# custom_colors depending on unique(type)
custom_colors <- c("-1" = "blue", "0" = "red", "50" = "green","100"="gray50")

theta_trans=dat$theta
theta_trans[which(theta_trans>pi)]=theta_trans[which(theta_trans>pi)]-2*pi
dat$theta_trans=abs(theta_trans)

windows()
ggplot(dat, aes(x = t, y = theta_trans)) +
  geom_point(aes(color = type)) +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_classic()

l = 5
times=seq(1,TT)
library(splus2R)
maxs = as.numeric(peaks(dat$theta_trans, span = l))   
mins = as.numeric(peaks(-dat$theta_trans, span = l)) 

# plot(dat$theta_trans,col=maxs+1,pch=19)
# plot(dat$theta_trans,col=mins+1,pch=19)

l2=50
#
count_min=runner::runner(
  mins,
  k = l2,
  f = mean,
  na_pad=T
)
count_min=count_min[-(1:floor(l2/2))]
count_min=c(count_min,rep(NA,round(l2/2)))
# plot(count_min,col=dat$type,ylab='count_min',xlab='t')
# legend('topright',legend=levels(dat$type),col=1:3,pch=19)
dat_min <- data.frame(
  t = seq_along(count_min),
  count_min = count_min,
  type = dat$type
)
count_max=runner::runner(
  maxs,
  k = l2,
  f = mean,
  na_pad=T
)
count_max=count_max[-(1:floor(l2/2))]
count_max=c(count_max,rep(NA,round(l2/2)))
dat_max <- data.frame(
  t = seq_along(count_max),
  count_max = count_max,
  type = dat$type
)

dat_max_min=merge(dat_max,dat_min,by=c('t','type'))

# Scatter plot using ggplot2 for count_min
ggplot(dat_max_min, aes(x = t, y = count_min, color = type)) +
  geom_point() +
  labs(x = "t", y = "count_min") +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_minimal() +
  theme(legend.position = "top")

# Scatter plot using ggplot2 for count_max
ggplot(dat_max_min, aes(x = t, y = count_max, color = type)) +
  geom_point() +
  labs(x = "t", y = "count_max") +
  scale_color_manual(values = custom_colors,name="Type") +
  theme_minimal() +
  theme(legend.position = "top")
# 

value_min=last_min_value(mins,dat$theta_trans)
dat_max_min$value_min=value_min
value_max=last_max_value(maxs,dat$theta_trans)
dat_max_min$value_max=value_max
dat_max_min$diffmaxmin=dat_max_min$value_max-dat_max_min$value_min
# 
ggplot(dat_max_min, aes(x = t, y = value_min, color = dat$type)) +
  geom_point() +
  labs(x = "t", y = "value_min") +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_minimal() +
  theme(legend.position = "top")

ggplot(dat_max_min, aes(x = t, y = value_max, color = dat$type)) +
  geom_point() +
  labs(x = "t", y = "value_max") +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_minimal() +
  theme(legend.position = "top")

ggplot(dat_max_min, aes(x = t, y = diffmaxmin, color = dat$type)) +
  geom_point() +
  labs(x = "t", y = "diffmaxmin") +
  scale_color_manual(values = custom_colors, name = "Type") +
  theme_minimal() +
  theme(legend.position = "top")
