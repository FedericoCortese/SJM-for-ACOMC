source("Utils3_.R")


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
windows()
annotate_figure(Ptrue2016CA138,
                top = text_grob("2016CA138 - true", color = "black", face = "bold", size = 14))

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
windows()
annotate_figure(Ptrue2016CO246,
                top = text_grob("2016CO246 - true", color = "black", face = "bold", size = 14))

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
windows()
annotate_figure(Ptrue2014OL339,
                top = text_grob("2014OL339 - true", color = "black", face = "bold", size = 14))


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
windows()
annotate_figure(Ptrue2015XX169,
                top = text_grob("2015XX169 - true", color = "black", face = "bold", size = 14))


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
windows()
annotate_figure(Ptrue2020PN1,
                top = text_grob("2020PN1 - true", color = "black", face = "bold", size = 14))
