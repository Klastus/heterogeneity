#### sprawdzenie czy rozdzial na small normal i big jest prawidłowy ####

spr.G1 <- area.G1
spr.G2 <- area.G2
spr.G1 <- subset(spr.G1, ! spr.G1$phase %in% c("outliers", "outliers1"))
spr.G2 <- subset(spr.G2, ! spr.G2$phase %in% c("outliers", "outliers1"))

spr <- rbind(spr.G1, spr.G2)
spr$size <- factor(spr$size, levels = c("small", "normal", "big"))

experiment <- "PT31-image"
a <- ggplot(dane.multip, aes(y= Intensity_MeanIntensity_Alexa,
                 x= factor(time.1.1)))+
  geom_boxplot()+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "pS1, mean, cells"))+
  ylim(0,200)
a

b <- ggplot(PT31image, aes(y= Intensity_IntegratedIntensity_Alexa,
                      x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size))+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "pS1, integrated, cells"))+
  ylim(0,2e6)
b

c <- ggplot(PT31image, aes(y= Intensity_MeanIntensity_DAPI,
                           x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size))+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "pS1, mean, cells"))
  #ylim(0,200)
c

d <- ggplot(PT31image, aes(y= Intensity_IntegratedIntensity_DAPI,
                           x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size))+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "pS1, integrated, cells"))
  #ylim(0,2e6)
d
b <- ggplot(spr, aes(y= Intensity_IntegratedIntensity_Alexa,
                     x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size))+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "size, integrted, cells"))+
  ylim(0,1e6)
b

c <- ggplot(pS1.means, aes(y= average_mean,
                     x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size))+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "size, mean, wells"))+
  ylim(0,200)
c

experiment <- "PT31"
d <- ggplot(pS1.means, aes(y= average_integrated,
                           x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size))+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "size, integrated, wells"))+
  ylim(0,0.5e6)
d

e <- ggplot(red, aes(y= Intensity_MeanIntensity_Alexa,
                     x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size, group=well.name))+
  facet_grid(.~size)+
  theme_jetka()+
  ggtitle(paste(experiment, "size, mean, cells"))+
  ylim(0,250)
e


pdf(paste(getwd(), "/R/model/output/", experiment, " area DAPI.pdf",sep=''), width=12,height=6)
a
b
c
d
dev.off()

experiment <- "PT31"
f <- ggplot(spr, aes(y= Intensity_MedianIntensity_Alexa,
                     x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size))+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "size, median, cells"))+
  ylim(0,200)
f
experiment <- "PT31"
g <- ggplot(spr, aes(y= (Intensity_UpperQuartileIntensity_Alexa-Intensity_LowerQuartileIntensity_Alexa),
                     x= factor(time.1.1)))+
  geom_boxplot(aes(fill=size))+
  facet_grid(phase~stimulation.1.1)+
  theme_jetka()+
  ggtitle(paste(experiment, "size, quartile, cells"))+
  ylim(0,200)
g

saveRDS(data.plot.31, paste(getwd(), "/R/model/output/", experiment, " wells.RDS",sep=''))
saveRDS(pS1.means, paste(getwd(), "/R/model/output/", experiment, " area wells.RDS",sep=''))
#
 blue <- subset(pS1.means, pS1.means$time.1.1==30 & pS1.means$stimulation.1.1==1 &
                 pS1.means$size=="big" & pS1.means$phase=="G2/M")

red <- subset(spr, spr$time.1.1==15 & spr$stimulation.1.1==1 &
                                 spr$phase=="G2/M")

#### 3D ####
install.packages("plotly")
library(plotly)
packageVersion('plotly')

spr$Intensity_MeanIntensity_Alexa
p <- plot_ly(spr, x= ~Intensity_MeanIntensity_Alexa, y=~Intensity_IntegratedIntensity_DAPI,
             z=~spr$AreaShape_Area) %>%
  add_marker(size=0.0001)
p
