source("Utils3_.R")

df164207=read.table("propagation_164207_new_v2.txt",header=T)
length(unique(df164207$type))
tt_164207=theta_trans_plot(df164207,'164207',l2=20)
tt_164207[[1]]
tt_164207[[2]]
tt_164207[[3]]
tt_164207[[4]]
tt_164207[[5]]
tt_164207[[6]]
tt_164207[[7]]
tt_164207[[8]]
tt_164207[[9]]
tt_164207[[10]]

df2001GO2=read.table("propagation_2001GO2_new_v2.txt",header=T)
length(unique(df2001GO2$type))
tt_2001GO2=theta_trans_plot(df2001GO2,'2001GO2',l2=50)
tt_2001GO2[[1]]
tt_2001GO2[[2]]
tt_2001GO2[[3]]
tt_2001GO2[[4]]
tt_2001GO2[[5]]
tt_2001GO2[[6]]


df2002AA29=read.table("propagation_2002AA29_new_v2.txt",header=T)
length(unique(df2002AA29$type))
tt_2002AA29=theta_trans_plot(df2002AA29,'2002AA29',l2=20)
tt_2002AA29[[1]]
tt_2002AA29[[2]]
tt_2002AA29[[3]]
tt_2002AA29[[4]]
tt_2002AA29[[5]]
tt_2002AA29[[6]]

# Il seguente ha una riga in meno (no ground truth)
df2005UH6=read.table("propagation_2005UH6_new_v2.txt",header=T)
df2005UH6$type=0
length(unique(df2005UH6$type))
tt_2005UH6=theta_trans_plot(df2005UH6,'2005UH6',l2=50,l3=50)
tt_2005UH6[[1]]
tt_2005UH6[[2]]
tt_2005UH6[[3]]
tt_2005UH6[[4]]
tt_2005UH6[[5]]
tt_2005UH6[[6]]
tt_2005UH6[[7]]
tt_2005UH6[[8]]
tt_2005UH6[[9]]

df2014OL339=read.table("propagation_2014OL339_new_v2.txt",header=T)
length(unique(df2014OL339$type))
tt_2014OL339=theta_trans_plot(df2014OL339,'2014OL339',l2=20,l3=150)
tt_2014OL339[[1]]
tt_2014OL339[[2]]
tt_2014OL339[[3]]
tt_2014OL339[[4]]
tt_2014OL339[[5]]
tt_2014OL339[[8]]
tt_2014OL339[[9]]

df2015SO2=read.table("propagation_2015SO2_new_v2.txt",header=T)
length(unique(df2015SO2$type))
tt_2015SO2=theta_trans_plot(df2015SO2,'2015SO2',l2=50)
tt_2015SO2[[1]]
tt_2015SO2[[2]]
tt_2015SO2[[3]]
tt_2015SO2[[4]]
tt_2015SO2[[5]]
tt_2015SO2[[6]]


df2015XX169=read.table("propagation_2015XX169_new_v2.txt",header=T)
length(unique(df2015XX169$type))
tt_2015XX169=theta_trans_plot(df2015XX169,'2015XX169',l2=20,l3=1000)
tt_2015XX169[[1]]
tt_2015XX169[[2]]
tt_2015XX169[[3]]
tt_2015XX169[[4]]
tt_2015XX169[[7]]
tt_2015XX169[[9]]

# MOST DIFFICULT
df2016CA138=read.table("propagation_2016CA138_new_v2.txt",header=T)
length(unique(df2016CA138$type))
tt_2016CA138=theta_trans_plot(df2016CA138,'2016CA138',l2=20)
tt_2016CA138[[1]]
tt_2016CA138[[2]]
tt_2016CA138[[3]]
tt_2016CA138[[4]]
tt_2016CA138[[7]]
tt_2016CA138[[8]]
tt_2016CA138[[9]]

df2016CO246=read.table("propagation_2016CO246_new_v2.txt",header=T)
length(unique(df2016CO246$type))
tt_2016CO246=theta_trans_plot(df2016CO246,'2016CO246',l2=20)
tt_2016CO246[[1]]
tt_2016CO246[[2]]
tt_2016CO246[[3]]
tt_2016CO246[[4]]
tt_2016CO246[[7]]

df2016HO3=read.table("propagation_2016HO3_new_v2.txt",header=T)
length(unique(df2016HO3$type))
tt_2016HO3=theta_trans_plot(df2016HO3,'2016HO3',l2=20)
tt_2016HO3[[1]]
tt_2016HO3[[2]]
tt_2016HO3[[3]]
tt_2016HO3[[4]]
tt_2016HO3[[7]]
tt_2016HO3[[9]]

df2019GM1=read.table("propagation_2019GM1_new_v2.txt",header=T)
length(unique(df2019GM1$type))
tt_2019GM1=theta_trans_plot(df2019GM1,'2019GM1',l2=20)
tt_2019GM1[[1]]
tt_2019GM1[[2]]
tt_2019GM1[[3]]
tt_2019GM1[[4]]
tt_2019GM1[[7]]
tt_2019GM1[[9]]

df2020PN1=read.table("propagation_2020PN1_new_v2.txt",header=T)
length(unique(df2020PN1$type))
tt_2020PN1=theta_trans_plot(df2020PN1,'2020PN1',l2=20,l3=50)
tt_2020PN1[[1]]
tt_2020PN1[[2]]
tt_2020PN1[[3]]
tt_2020PN1[[4]]
tt_2020PN1[[7]]
tt_2020PN1[[8]]
tt_2020PN1[[9]]


df2020PP1=read.table("propagation_2020PP1_new_v2.txt",header=T)
length(unique(df2020PP1$type))
tt_2020PP1=theta_trans_plot(df2020PP1,'2020PP1',l2=20)
tt_2020PP1[[1]]
tt_2020PP1[[2]]
tt_2020PP1[[3]]
tt_2020PP1[[4]]
tt_2020PP1[[7]]
tt_2020PP1[[9]]
