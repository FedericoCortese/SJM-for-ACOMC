library(ggplot2)
library(ggpubr)
library(dplyr)
library(zoo)
source("Utils3_.R")
# 2015XX169 ---------------------------------------------------------------

df2015XX169=read.table("propagation_2015XX169_new_v2.txt",header=T)

p2015XX169=plot_real4(df2015XX169)

Ptrue2015XX169=ggarrange(p2015XX169$Pa,
                         p2015XX169$Pe,
                         p2015XX169$Ptheta,
                         p2015XX169$Pomega,
                         nrow=4
                         #,main="164207"
)
png(filename = "2015XX169TRUE.png",width = 1600, height = 900)
annotate_figure(Ptrue2015XX169,
                top = text_grob("2015XX169 - true", color = "black", face = "bold", size = 14))
dev.off()
true_states=df2015XX169$type
# TRUE STATES


df2015XX169=df2015XX169[,c("a","e","theta","omega")]
#df2015XX169=compute_feat(df2015XX169,wdn=10,am1=T)
df2015XX169=compute_feat(df2015XX169,am1=T)


lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2015XX169)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2015XX169,function(x)x$FTIC))
# unlist(lapply(est2015XX169,function(x)x$overlap))

load("est2015XX169.RData")

modsel2015XX169=data.frame(hp,FTIC=unlist(lapply(est2015XX169,function(x)x$FTIC)),
                           overlap=unlist(lapply(est2015XX169,function(x)x$overlap)),
                           ARI=unlist(lapply(est2015XX169,function(x)x$ARI))
)

modsel2015XX169

estw2015XX169=data.frame(weight=est2015XX169[[1]]$est_weights,
                         feat=colnames(df2015XX169))

# Sort by weight
estw2015XX169=estw2015XX169[order(estw2015XX169$weight,decreasing = T),]
estw2015XX169

# windows()
# par(mfrow=c(4,1))
# plot(df2015XX169$sddtheta,type='l',main="2015XX169",ylab="sddtheta")
# plot(df2015XX169$sdwdn_omega,type='l',ylab="sdwdn_omega")
# plot(df2015XX169$sdwdn_a,type='l',ylab="sdwdn_a")
# plot(df2015XX169$sdwdn_theta,type='l',ylab="sdwdn_theta")

df2015XX169_b=data.frame(df2015XX169,type=est2015XX169[[12]]$est_states)
p2015XX169_b=plot_real4(df2015XX169_b)

Pest2015XX169_b=ggarrange(p2015XX169_b$Pa,
                          p2015XX169_b$Pe,
                          p2015XX169_b$Ptheta,
                          p2015XX169_b$Pomega,
                          nrow=4
                          #,main="2015XX169"
)

png(filename = "2015XX169EST.png",width = 1600, height = 900)
annotate_figure(Pest2015XX169_b,
                top = text_grob("2015XX169 - est", color = "black", face = "bold", size = 14))
dev.off()



# 2016CA138 ---------------------------------------------------------------
df2016CA138=read.table("propagation_2016CA138_new_v2.txt",header=T)

p2016CA138=plot_real4(df2016CA138)

Ptrue2016CA138=ggarrange(p2016CA138$Pa,
                      p2016CA138$Pe,
                      p2016CA138$Ptheta,
                      p2016CA138$Pomega,
                      nrow=4
                      #,main="164207"
)

png(filename = "2016CA138TRUE.png",width = 1600, height = 900)
annotate_figure(Ptrue2016CA138,
                top = text_grob("2016CA138 - true", color = "black", face = "bold", size = 14))
dev.off()
true_states=df2016CA138$type
# TRUE STATES


df2016CA138=df2016CA138[,c("a","e","theta","omega")]
df2016CA138=compute_feat(df2016CA138,am1=T)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2016CA138)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2016CA138,function(x)x$FTIC))
# unlist(lapply(est2016CA138,function(x)x$overlap))

load("est2016CA138.RData")

modsel2016CA138=data.frame(hp,FTIC=unlist(lapply(est2016CA138,function(x)x$FTIC)),
                           overlap=unlist(lapply(est2016CA138,function(x)x$overlap)),
                           ARI=unlist(lapply(est2016CA138,function(x)x$ARI))
)

modsel2016CA138
sel=4
estw2016CA138=data.frame(weight=est2016CA138[[sel]]$est_weights,
                         feat=colnames(df2016CA138))

# Sort by weight
estw2016CA138=estw2016CA138[order(estw2016CA138$weight,decreasing = T),]
estw2016CA138

# windows()
# par(mfrow=c(4,1))
# plot(df2016CA138$sddtheta,type='l',main="2016CA138",ylab="sddtheta")
# plot(df2016CA138$sdwdn_omega,type='l',ylab="sdwdn_omega")
# plot(df2016CA138$sdwdn_a,type='l',ylab="sdwdn_a")
# plot(df2016CA138$sdwdn_theta,type='l',ylab="sdwdn_theta")

df2016CA138_b=data.frame(df2016CA138,type=est2016CA138[[sel]]$est_states)
p2016CA138_b=plot_real4(df2016CA138_b)

Pest2016CA138_b=ggarrange(p2016CA138_b$Pa,
                          p2016CA138_b$Pe,
                          p2016CA138_b$Ptheta,
                          p2016CA138_b$Pomega,
                          nrow=4
                          #,main="2016CA138"
)
windows()
annotate_figure(Pest2016CA138_b,
                top = text_grob("2016CA138 - est", color = "black", face = "bold", size = 14))


# 2016CO246 ---------------------------------------------------------------

df2016CO246=read.table("propagation_2016CO246_new_v2.txt",header=T)

p2016CO246=plot_real4(df2016CO246)

Ptrue2016CO246=ggarrange(p2016CO246$Pa,
                      p2016CO246$Pe,
                      p2016CO246$Ptheta,
                      p2016CO246$Pomega,
                      nrow=4
                      #,main="164207"
)

png(filename = "2016CO246TRUE.png",width = 1600, height = 900)
annotate_figure(Ptrue2016CO246,
                top = text_grob("2016CO246 - true", color = "black", face = "bold", size = 14))
dev.off()
true_states=df2016CO246$type
# TRUE STATES


df2016CO246=df2016CO246[,c("a","e","theta","omega")]
df2016CO246=compute_feat(df2016CO246,am1=T)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2016CO246)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2016CO246,function(x)x$FTIC))
# unlist(lapply(est2016CO246,function(x)x$overlap))

load("est2016CO246.RData")

modsel2016CO246=data.frame(hp,FTIC=unlist(lapply(est2016CO246,function(x)x$FTIC)),
                           overlap=unlist(lapply(est2016CO246,function(x)x$overlap)),
                           ARI=unlist(lapply(est2016CO246,function(x)x$ARI))
)

modsel2016CO246

sel=13

estw2016CO246=data.frame(weight=est2016CO246[[sel]]$est_weights,
                         feat=colnames(df2016CO246))

# Sort by weight
estw2016CO246=estw2016CO246[order(estw2016CO246$weight,decreasing = T),]
estw2016CO246

# windows()
# par(mfrow=c(4,1))
# plot(df2016CO246$sddtheta,type='l',main="2016CO246",ylab="sddtheta")
# plot(df2016CO246$sdwdn_omega,type='l',ylab="sdwdn_omega")
# plot(df2016CO246$sdwdn_a,type='l',ylab="sdwdn_a")
# plot(df2016CO246$sdwdn_theta,type='l',ylab="sdwdn_theta")

df2016CO246_b=data.frame(df2016CO246,type=est2016CO246[[sel]]$est_states)
p2016CO246_b=plot_real4(df2016CO246_b)

Pest2016CO246_b=ggarrange(p2016CO246_b$Pa,
                          p2016CO246_b$Pe,
                          p2016CO246_b$Ptheta,
                          p2016CO246_b$Pomega,
                          nrow=4
                          #,main="2016CO246"
)

png(filename = "2016CO246EST.png",width = 1600, height = 900)
annotate_figure(Pest2016CO246_b,
                top = text_grob("2016CO246 - est", color = "black", face = "bold", size = 14))
dev.off()


# 2014OL339 ---------------------------------------------------------------
df2014OL339=read.table("propagation_2014OL339_new_v2.txt",header=T)

p2014OL339=plot_real4(df2014OL339)

Ptrue2014OL339=ggarrange(p2014OL339$Pa,
                      p2014OL339$Pe,
                      p2014OL339$Ptheta,
                      p2014OL339$Pomega,
                      nrow=4
                      #,main="164207"
)

png(filename = "2014OL339TRUE.png",width = 1600, height = 900)
annotate_figure(Ptrue2014OL339,
                top = text_grob("2014OL339 - true", color = "black", face = "bold", size = 14))
dev.off()
true_states=df2014OL339$type
# TRUE STATES


df2014OL339=df2014OL339[,c("a","e","theta","omega")]
df2014OL339=compute_feat(df2014OL339,am1=T)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2014OL339)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2014OL339,function(x)x$FTIC))
# unlist(lapply(est2014OL339,function(x)x$overlap))

load("est2014OL339.RData")

modsel2014OL339=data.frame(hp,FTIC=unlist(lapply(est2014OL339,function(x)x$FTIC)),
                           overlap=unlist(lapply(est2014OL339,function(x)x$overlap)),
                           ARI=unlist(lapply(est2014OL339,function(x)x$ARI))
)

modsel2014OL339
sel=17
estw2014OL339=data.frame(weight=est2014OL339[[sel]]$est_weights,
                         feat=colnames(df2014OL339))

# Sort by weight
estw2014OL339=estw2014OL339[order(estw2014OL339$weight,decreasing = T),]
estw2014OL339

# windows()
# par(mfrow=c(4,1))
# plot(df2014OL339$sddtheta,type='l',main="2014OL339",ylab="sddtheta")
# plot(df2014OL339$sdwdn_omega,type='l',ylab="sdwdn_omega")
# plot(df2014OL339$sdwdn_a,type='l',ylab="sdwdn_a")
# plot(df2014OL339$sdwdn_theta,type='l',ylab="sdwdn_theta")

df2014OL339_b=data.frame(df2014OL339,type=est2014OL339[[sel]]$est_states)
p2014OL339_b=plot_real4(df2014OL339_b)

Pest2014OL339_b=ggarrange(p2014OL339_b$Pa,
                          p2014OL339_b$Pe,
                          p2014OL339_b$Ptheta,
                          p2014OL339_b$Pomega,
                          nrow=4
                          #,main="2014OL339"
)

png(filename = "2014OL339EST.png",width = 1600, height = 900)
annotate_figure(Pest2014OL339_b,
                top = text_grob("2014OL339 - est", color = "black", face = "bold", size = 14))
dev.off()


# 2020PN1 -----------------------------------------------------------------

df2020PN1=read.table("propagation_2020PN1_new_v2.txt",header=T)

p2020PN1=plot_real4(df2020PN1)

Ptrue2020PN1=ggarrange(p2020PN1$Pa,
                      p2020PN1$Pe,
                      p2020PN1$Ptheta,
                      p2020PN1$Pomega,
                      nrow=4
                      #,main="164207"
)

png(filename = "2020PN1TRUE.png",width = 1600, height = 900)
annotate_figure(Ptrue2020PN1,
                top = text_grob("2020PN1 - true", color = "black", face = "bold", size = 14))
dev.off()

true_states=df2020PN1$type
# TRUE STATES


df2020PN1=df2020PN1[,c("a","e","theta","omega")]
df2020PN1=compute_feat(df2020PN1,am1=T)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2020PN1)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2020PN1,function(x)x$FTIC))
# unlist(lapply(est2020PN1,function(x)x$overlap))

load("est2020PN1.RData")

modsel2020PN1=data.frame(hp,FTIC=unlist(lapply(est2020PN1,function(x)x$FTIC)),
                         overlap=unlist(lapply(est2020PN1,function(x)x$overlap)),
                         ARI=unlist(lapply(est2020PN1,function(x)x$ARI))
)

modsel2020PN1

sel=12
estw2020PN1=data.frame(weight=est2020PN1[[sel]]$est_weights,
                       feat=colnames(df2020PN1))

# Sort by weight
estw2020PN1=estw2020PN1[order(estw2020PN1$weight,decreasing = T),]
estw2020PN1

# windows()
# par(mfrow=c(4,1))
# plot(df2020PN1$sddtheta,type='l',main="2020PN1",ylab="sddtheta")
# plot(df2020PN1$sdwdn_omega,type='l',ylab="sdwdn_omega")
# plot(df2020PN1$sdwdn_a,type='l',ylab="sdwdn_a")
# plot(df2020PN1$sdwdn_theta,type='l',ylab="sdwdn_theta")


df2020PN1_b=data.frame(df2020PN1,type=est2020PN1[[sel]]$est_states)
p2020PN1_b=plot_real4(df2020PN1_b)

Pest2020PN1_b=ggarrange(p2020PN1_b$Pa,
                        p2020PN1_b$Pe,
                        p2020PN1_b$Ptheta,
                        p2020PN1_b$Pomega,
                        nrow=4
                        #,main="2020PN1"
)
png(filename = "2020PN1EST.png",width = 1600, height = 900)
annotate_figure(Pest2020PN1_b,
                top = text_grob("2020PN1 - est", color = "black", face = "bold", size = 14))
dev.off()
