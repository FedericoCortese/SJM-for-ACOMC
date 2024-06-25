source("Utils3_.R")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(zoo)

# 164207 ------------------------------------------------------------------

df164207=read.table("propagation_164207_new_v2.txt",header=T)

p164207=plot_real(df164207,"blue")

Ptrue164207=ggarrange(p164207$Pa,
                      p164207$Pe,
                      p164207$Ptheta,
                      p164207$Pomega,
                      nrow=4
                      #,main="164207"
)
windows()
annotate_figure(Ptrue164207,
                top = text_grob("164207 - true", color = "black", face = "bold", size = 14))


true_states=df164207$type
# TRUE STATES


df164207=df164207[,c("a","e","theta","omega")]
df164207=compute_feat(df164207)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df164207)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est164207,function(x)x$FTIC))
# unlist(lapply(est164207,function(x)x$overlap))

modsel164207=data.frame(hp,FTIC=unlist(lapply(est164207,function(x)x$FTIC)),
                        overlap=unlist(lapply(est164207,function(x)x$overlap)),
                        ARI=unlist(lapply(est164207,function(x)x$ARI))
)

modsel164207

estw164207=data.frame(weight=est164207[[6]]$est_weights,feat=colnames(df164207))

# Sort by weight
estw164207=estw164207[order(estw164207$weight,decreasing = T),]
estw164207

windows()
par(mfrow=c(4,1))
plot(df164207$sddtheta,type='l',main="164207",ylab="sddtheta")
plot(df164207$sdwdn_omega,type='l',ylab="sdwdn_omega")
plot(df164207$sdwdn_a,type='l',ylab="sdwdn_a")
plot(df164207$sdwdn_theta,type='l',ylab="sdwdn_theta")

df164207_b=data.frame(df164207,type=est164207[[6]]$est_states)
p164207_b=plot_real(df164207_b,"red")

Pest164207_b=ggarrange(p164207_b$Pa,
                      p164207_b$Pe,
                      p164207_b$Ptheta,
                      p164207_b$Pomega,
                      nrow=4
                      #,main="164207"
)
windows()
annotate_figure(Pest164207_b,
                top = text_grob("164207 - est", color = "black", face = "bold", size = 14))


# 2001GO2 -----------------------------------------------------------------

df2001GO2=read.table("propagation_2001GO2_new_v2.txt",header=T)

p2001GO2=plot_real(df2001GO2,"blue")

Ptrue2001GO2=ggarrange(p2001GO2$Pa,
                       p2001GO2$Pe,
                       p2001GO2$Ptheta,
                       p2001GO2$Pomega,
                       nrow=4
                       #,main="2001GO2"
)
windows()
annotate_figure(Ptrue2001GO2,
                top = text_grob("2001GO2 - true", color = "black", face = "bold", size = 14))


true_states=df2001GO2$type
# TRUE STATES


df2001GO2=df2001GO2[,c("a","e","theta","omega")]
df2001GO2=compute_feat(df2001GO2)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2001GO2)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2001GO2,function(x)x$FTIC))
# unlist(lapply(est2001GO2,function(x)x$overlap))

modsel2001GO2=data.frame(hp,FTIC=unlist(lapply(est2001GO2,function(x)x$FTIC)),
                         overlap=unlist(lapply(est2001GO2,function(x)x$overlap)),
                         ARI=unlist(lapply(est2001GO2,function(x)x$ARI))
)

modsel2001GO2

estw2001GO2=data.frame(weight=est2001GO2[[9]]$est_weights,feat=colnames(df2001GO2))

# Sort by weight
estw2001GO2=estw2001GO2[order(estw2001GO2$weight,decreasing = T),]
estw2001GO2

windows()
par(mfrow=c(4,1))
plot(df2001GO2$a,type='l',main="2001GO2",ylab="a")
plot(df2001GO2$sddtheta,type='l',ylab="sddtheta")
plot(df2001GO2$sdwdn_theta,type='l',ylab="sdwdn_theta")
plot(df2001GO2$sdde,type='l',ylab="sdde")

df2001GO2_b=data.frame(df2001GO2,type=est2001GO2[[9]]$est_states)
p2001GO2_b=plot_real(df2001GO2_b,"red")

Pest2001GO2_b=ggarrange(p2001GO2_b$Pa,
                        p2001GO2_b$Pe,
                        p2001GO2_b$Ptheta,
                        p2001GO2_b$Pomega,
                        nrow=4
                        #,main="2001GO2"
)
windows()
annotate_figure(Pest2001GO2_b,
                top = text_grob("2001GO2 - est", color = "black", face = "bold", size = 14))


# 2002AA29 --------------------------------------------------------------------

df2002AA29=read.table("propagation_2002AA29_new_v2.txt",header=T)

# 2016HO3 -----------------------------------------------------------------

df2016HO3=read.table("propagation_2016HO3_new_v2.txt",header=T)
unique(df2016HO3$type)

df2016HO3=read.table("propagation_2016HO3_new_v2.txt",header=T)

p2016HO3=plot_real(df2016HO3,"blue")

Ptrue2016HO3=ggarrange(p2016HO3$Pa,
                       p2016HO3$Pe,
                       p2016HO3$Ptheta,
                       p2016HO3$Pomega,
                       nrow=4
                       #,main="2016HO3"
)
windows()
annotate_figure(Ptrue2016HO3,
                top = text_grob("2016HO3 - true", color = "black", face = "bold", size = 14))


true_states=df2016HO3$type
# TRUE STATES


df2016HO3=df2016HO3[,c("a","e","theta","omega")]
df2016HO3=compute_feat(df2016HO3)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2016HO3)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2016HO3,function(x)x$FTIC))
# unlist(lapply(est2016HO3,function(x)x$overlap))

modsel2016HO3=data.frame(hp,FTIC=unlist(lapply(est2016HO3,function(x)x$FTIC)),
                         overlap=unlist(lapply(est2016HO3,function(x)x$overlap)),
                         ARI=unlist(lapply(est2016HO3,function(x)x$ARI))
)

modsel2016HO3

estw2016HO3=data.frame(weight=est2016HO3[[8]]$est_weights,feat=colnames(df2016HO3))

# Sort by weight
estw2016HO3=estw2016HO3[order(estw2016HO3$weight,decreasing = T),]
estw2016HO3

windows()
par(mfrow=c(4,1))
plot(df2016HO3$sddtheta,type='l',main="2016HO3",ylab="sddtheta")
plot(df2016HO3$sdwdn_omega,type='l',ylab="sdwdn_omega")
plot(df2016HO3$sdwdn_a,type='l',ylab="sdwdn_a")
plot(df2016HO3$sdwdn_theta,type='l',ylab="sdwdn_theta")

df2016HO3_b=data.frame(df2016HO3,type=est2016HO3[[8]]$est_states)
p2016HO3_b=plot_real(df2016HO3_b,"red")

Pest2016HO3_b=ggarrange(p2016HO3_b$Pa,
                        p2016HO3_b$Pe,
                        p2016HO3_b$Ptheta,
                        p2016HO3_b$Pomega,
                        nrow=4
                        #,main="2016HO3"
)
windows()
annotate_figure(Pest2016HO3_b,
                top = text_grob("2016HO3 - est", color = "black", face = "bold", size = 14))



# 2019GM1 -----------------------------------------------------------------

df2019GM1=read.table("propagation_2019GM1_new_v2.txt",header=T)
unique(df2019GM1$type)

df2019GM1=read.table("propagation_2019GM1_new_v2.txt",header=T)

p2019GM1=plot_real(df2019GM1,"blue")

Ptrue2019GM1=ggarrange(p2019GM1$Pa,
                       p2019GM1$Pe,
                       p2019GM1$Ptheta,
                       p2019GM1$Pomega,
                       nrow=4
                       #,main="2019GM1"
)
windows()
annotate_figure(Ptrue2019GM1,
                top = text_grob("2019GM1 - true", color = "black", face = "bold", size = 14))


true_states=df2019GM1$type
# TRUE STATES


df2019GM1=df2019GM1[,c("a","e","theta","omega")]
df2019GM1=compute_feat(df2019GM1)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2019GM1)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2019GM1,function(x)x$FTIC))
# unlist(lapply(est2019GM1,function(x)x$overlap))

modsel2019GM1=data.frame(hp,FTIC=unlist(lapply(est2019GM1,function(x)x$FTIC)),
                         overlap=unlist(lapply(est2019GM1,function(x)x$overlap)),
                         ARI=unlist(lapply(est2019GM1,function(x)x$ARI))
)

modsel2019GM1

estw2019GM1=data.frame(weight=est2019GM1[[8]]$est_weights,feat=colnames(df2019GM1))

# Sort by weight
estw2019GM1=estw2019GM1[order(estw2019GM1$weight,decreasing = T),]
estw2019GM1

windows()
par(mfrow=c(4,1))
plot(df2019GM1$sddtheta,type='l',main="2019GM1",ylab="sddtheta")
plot(df2019GM1$sdwdn_omega,type='l',ylab="sdwdn_omega")
plot(df2019GM1$sdwdn_a,type='l',ylab="sdwdn_a")
plot(df2019GM1$sdwdn_theta,type='l',ylab="sdwdn_theta")

df2019GM1_b=data.frame(df2019GM1,type=est2019GM1[[8]]$est_states)
p2019GM1_b=plot_real(df2019GM1_b,"red")

Pest2019GM1_b=ggarrange(p2019GM1_b$Pa,
                        p2019GM1_b$Pe,
                        p2019GM1_b$Ptheta,
                        p2019GM1_b$Pomega,
                        nrow=4
                        #,main="2019GM1"
)
windows()
annotate_figure(Pest2019GM1_b,
                top = text_grob("2019GM1 - est", color = "black", face = "bold", size = 14))


# 2020PP1 -----------------------------------------------------------------

df2020PP1=read.table("propagation_2020PP1_new_v2.txt",header=T)
unique(df2020PP1$type)

p2020PP1=plot_real(df2020PP1,"blue")

Ptrue2020PP1=ggarrange(p2020PP1$Pa,
                       p2020PP1$Pe,
                       p2020PP1$Ptheta,
                       p2020PP1$Pomega,
                       nrow=4
                       #,main="2020PP1"
)
windows()
annotate_figure(Ptrue2020PP1,
                top = text_grob("2020PP1 - true", color = "black", face = "bold", size = 14))


true_states=df2020PP1$type
# TRUE STATES


df2020PP1=df2020PP1[,c("a","e","theta","omega")]
df2020PP1=compute_feat(df2020PP1)

lambda=c(0,5,10,15,20,30)
kappa=seq(1,ceiling(sqrt(dim(df2020PP1)[2])),by=1)
hp=expand.grid(lambda=lambda,kappa=kappa)
# unlist(lapply(est2020PP1,function(x)x$FTIC))
# unlist(lapply(est2020PP1,function(x)x$overlap))

modsel2020PP1=data.frame(hp,FTIC=unlist(lapply(est2020PP1,function(x)x$FTIC)),
                         overlap=unlist(lapply(est2020PP1,function(x)x$overlap)),
                         ARI=unlist(lapply(est2020PP1,function(x)x$ARI))
)

modsel2020PP1

estw2020PP1=data.frame(weight=est2020PP1[[8]]$est_weights,feat=colnames(df2020PP1))

# Sort by weight
estw2020PP1=estw2020PP1[order(estw2020PP1$weight,decreasing = T),]
estw2020PP1

windows()
par(mfrow=c(4,1))
plot(df2020PP1$sddtheta,type='l',main="2020PP1",ylab="sddtheta")
plot(df2020PP1$sdwdn_omega,type='l',ylab="sdwdn_omega")
plot(df2020PP1$sdwdn_a,type='l',ylab="sdwdn_a")
plot(df2020PP1$sdwdn_theta,type='l',ylab="sdwdn_theta")

df2020PP1_b=data.frame(df2020PP1,type=est2020PP1[[8]]$est_states)
p2020PP1_b=plot_real(df2020PP1_b,"red")

Pest2020PP1_b=ggarrange(p2020PP1_b$Pa,
                        p2020PP1_b$Pe,
                        p2020PP1_b$Ptheta,
                        p2020PP1_b$Pomega,
                        nrow=4
                        #,main="2020PP1"
)
windows()
annotate_figure(Pest2020PP1_b,
                top = text_grob("2020PP1 - est", color = "black", face = "bold", size = 14))
