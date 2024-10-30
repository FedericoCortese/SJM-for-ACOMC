df100000915=read.table("100000915.txt",header = T)

# Only HS and QS

start=which.min(abs(df100000915$t-193792105.70 ))
end=which.min(abs(df100000915$t-258552105.70))

df100000915=df100000915[start:end,]

summary(diff(df100000915$t))
#plot(diff(df100000915$t),type='p')

#df100000915$type=-1
data=df100000915
data_name="100000915"

library(ggplot2)
library(dplyr)
library(plotly)

source("Utils4_.R")

max_lag=10
min_lag=5
# %FC importante che questa soglia sia abbastanza grande da distinguere QS dal resto
tt_thres_diffmaxmin=0.07

data_thetavol=thetavol_feat(data,
                            wdn=c(min_lag,max_lag))
data_a=a_feat(data,l3=c(min_lag,max_lag))

data_maxmin=max_min_feat(data,
                         tt_thres_maxmin = 2.4,
                         l2=c(min_lag,max_lag),
                         tt_thres_diffmaxmin=tt_thres_diffmaxmin)
names(data_maxmin)

#data_maxmin=data_maxmin[,c("t","I_TD","I_HS","I_QS","mean_osc")]

# data_maxmin=data_maxmin[,c("t","mean_osc")]
# %FC provo a reintrodurre value max min, e diff max min
# data_maxmin=data_maxmin[,c("t","mean_osc","value_min" ,"value_max","diffmaxmin","I_diffmaxmin")]

data_fin=merge(data_thetavol,data_a,by="t")
data_fin=merge(data_fin,data_maxmin,by="t")

names(data_fin)

Y=data_fin[complete.cases(data_fin),]
Y=Y[,-1]

data2=merge(data_fin,df100000915,by="t")
data2 <- data2 %>%
  mutate(Segment_max = cumsum(maxs != lag(maxs, default = first(maxs))),
         Segment_min = cumsum(mins != lag(mins, default = first(mins))),
         maxs_mins=mins+maxs*2)

data2 <- data2 %>%
  group_by(Segment_max) %>%
  mutate(next_t_max = dplyr::lead(t), next_theta_max = dplyr::lead(theta))

data2 <- data2 %>%
  group_by(Segment_min) %>%
  mutate(next_t_min = dplyr::lead(t), next_theta_min = dplyr::lead(theta))


p_I <- ggplot(data2) + 
  labs(title = "100000915", x = "t", y = "Theta") +
  geom_line(aes(x=t,y=theta),col='grey')+
  geom_point(aes(x=t,y=theta),col='grey40',size=1)+
  geom_point(aes(x = t, y = theta, color = factor(maxs_mins)), size = 1) +
  scale_color_manual(values = c("0" = "grey50", "1"='blue' ,"2" = "red")) +
  geom_line(aes(x=t,y=I_TD),col='green4')+
  geom_line(aes(x=t,y=I_HS),col='cyan3')+
  geom_line(aes(x=t,y=I_QS),col="violet")+
  #geom_line(aes(x=t,y=I_diffmaxmin*2),col='black')+
  theme_minimal() +
  theme(legend.position = "none") 
#p_I
# Convert the ggplot object to a plotly interactive plot
ggplotly(p_I)
