#PT.Win:
setwd('Y:/PiotrT/cellcyclemodel')
#PT.MAC
setwd("/Users/piotrt/Documents/IPPT PT/JAK-STAT model/cellcycle")

prepare_data_area <- function (input.folder="PT31", which.phase.1=1,
                               which.phase.2=3, bord.1=2000, bord.2=2800,
                               bord.1.G2=2800, bord.2.G2=4200, to.wells=0, 
                               phase.list, multiply=1){
  ####data loading ####
  dane <- read.csv(paste(getwd(), "R", input.folder, 'ShrinkedNuclei.csv', sep = '/'),
                   sep = ',',header = 1)
  dane <- normalize_data(dane)$data
  

  if ("antibody.1.3" %in% colnames(dane)) {
    dane <- subset(dane, dane$antibody.1.3 ==  "pS1")
  }
  dane.cc <- dane
  dane.cc$phase <- ifelse(dane.cc$Intensity_IntegratedIntensity_DAPI < phase.list[[1]], "outliers1",
                          ifelse(dane.cc$Intensity_IntegratedIntensity_DAPI < phase.list[[2]], "G1",
                                 ifelse(dane.cc$Intensity_IntegratedIntensity_DAPI < phase.list[[3]], "S",
                                        ifelse(dane.cc$Intensity_IntegratedIntensity_DAPI < phase.list[[4]],
                                               "G2/M", "outliers"))))
  
  # multiplying data of time=0 to get stim=0 in every time
  if (multiply==0){
  dane.multip <- dane.cc} else {
    stim0<-dane.cc[dane.cc$stimulation.1.1 == 0, ]
    stim0.time<-dane.cc[dane.cc$stimulation.1.1 == 0, ]
    list.stim0 <- list(dane.cc, stim0, stim0, stim0.time, stim0.time)
    
    list.stim0[[2]]$stimulation.1.1 <- 1
    list.stim0[[3]]$stimulation.1.1 <- 0.1
    list.stim0[[4]]$time.1.1 <- 15
    list.stim0[[5]]$time.1.1 <- 90
    
    dane.multip <<- do.call("rbind", list.stim0)
  }  
  
  #dane.multip <- dane.cc
  ###preparing data to treat it as one well=one observation###
  pS1.wells <- dane.multip[, c("well.name", "Intensity_IntegratedIntensity_Alexa",
                               "Intensity_MeanIntensity_DAPI",
                               "Intensity_IntegratedIntensity_DAPI",
                               "Intensity_MeanIntensity_Alexa", "phase", "time.1.1", 
                               "stimulation.1.1", "AreaShape_Area")]

  phases <- c("G1", "S", "G2/M")
  area.G1 <- subset(pS1.wells, pS1.wells$phase == phases[[which.phase.1]])
  area.G1$size <- ifelse(area.G1$AreaShape_Area < bord.1, "small",
                         ifelse(area.G1$AreaShape_Area > bord.2 , "big", "normal"))
  
  area.G2 <- subset(pS1.wells, pS1.wells$phase == phases[[which.phase.2]])
  area.G2$size <- ifelse(area.G2$AreaShape_Area < bord.1.G2, "small",
                         ifelse(area.G2$AreaShape_Area > bord.2.G2 , "big", "normal"))
  
  area <- rbind(area.G1, area.G2)
  wells <- unique(area$well.name)
  size <- c("small", "big", "normal")
  pS1.means <- data.frame()
  if (to.wells==0){
    area$size <- factor(area$size, levels = c("small", "normal", "big"))
    area } else{
  for(ph in c("G1", "G2/M")){
  for(well in wells){
    for(si in size){
      pS1.subset <- area[area$well.name == well &
                              area$size == si & area$phase == ph,]
      for (stim in unique(pS1.subset$stimulation.1.1)){
        pS1.rbind <- pS1.subset[pS1.subset$stimulation.1.1 == stim, ]
        pS1.rbind$average_integrated <- mean(pS1.rbind$Intensity_IntegratedIntensity_Alexa)
        pS1.rbind$average_mean <- mean(pS1.rbind$Intensity_MeanIntensity_Alexa)
        pS1.means <- rbind(pS1.means, pS1.rbind[1, ])
      }
    }
  } 
  }  
  pS1.means$size <- factor(pS1.means$size, levels = c("small", "normal", "big"))
  pS1.means
   
  }
  }
##### ###### 
#plotting
experiment <- "PT29 c/nc"
dane.multip.G1 <- dane.multip[dane.multip$phase=="G2/M",]
t <- ggplot(dane.multip.G1, aes(x = AreaShape_Area))+
  geom_histogram(bins=100)+
  xlim(0,6e3)+
  geom_vline(xintercept=2800, col="brown")+
  geom_vline(xintercept=4050, col="green")+
  theme_jetka(base_size = 7)+
  ggtitle(paste(experiment, "G2/M cells"))
t

test <- prepare_data_area(input.folder="PT31", which.phase.1=1, which.phase.2=3, 
                          bord.1=2000, bord.2=2800, bord.1.G2=2800, bord.2.G2=4200,
                          phase.list=list("out1"=0.9e6, "G1"=1.5e6, "S"=1.95e6, "G2/M"=2.8e6),
                          to.wells=0)

PT31image <- prepare_data_area(input.folder="PT31-image", which.phase.1=1, which.phase.2=3,
                          bord.1=1900, bord.2=2800, bord.1.G2=2800, bord.2.G2=4200,
                          phase.list=list("out1"=1.6e6, "G1"=2.6e6, "S"=3.3e6, "G2/M"=4.3e6),
                          to.wells=0, multiply=0)
                                  

# dla PT29:
# G2 bord.1=3900, bord.2=6000 ?
# G1 bord.1=3300, bord.2=4600 ?
# dla PT31:
# G2 bord.1=2800, bord.2=4200 
# G1 bord.1=2000, bord.2=2800
# dla PT31-image:
# G2 bord.1=2800, bord.2=4200 
# G1 bord.1=1900, bord.2=2800
# PT29 c/nc
# G2 bord.1=2800, bord.2=4050 
# G1 bord.1=1950, bord.2=2800

