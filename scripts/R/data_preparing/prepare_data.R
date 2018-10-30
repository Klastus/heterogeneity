#Data preparing script:
prepare_data <- function (input.folder="Y:/PiotrT/conferences/Japan/poster/input",
                          exp.name="PT31",
                          plot.hist=0, 
                          phase.list, 
                          to.wells,
                          which.phase.1=0, which.phase.2=which.phase.1){
  # prepare_data function prepares .csv experimental data for plotting, that is:
  # assigns each cell to a cell phase, multiplies  control data, 
  # averages each well and subsets the data according to cell phase
  #
  # Args:
  # input.folder- folder where the .csv data is located
  # plot.hist- logical, should the histogram be drawn (helpful in 
  # deciding the borders of cc)
  # phase.list- have to be a LIST!!! it lists all the borders of cc
  # which.phase.1- according to which phase the data should be subsetted. 
  # add which.phase.2 argument if you want to subset according to two phases
  
  ####data loading ####
  dane <- read.csv(paste(input.folder, 'ShrinkedNuclei.csv', sep = '/'),
                   sep = ',',header = 1)
  dane <- normalize_data(dane)$data
  
  #### dane data processing ####
  
  #cell cycle phases identification
  if ("antibody.1.3" %in% colnames(dane)) {
    dane <- subset(dane, dane$antibody.1.3 ==  "pS1")
  }
  
  
  if(plot.hist == 1){
    hist.cycle <- paste(exp.name, "cycle", sep = ' ')
    hist <- ggplot(dane, aes(Intensity_IntegratedIntensity_DAPI))+
      geom_histogram(bins=100)+
      geom_vline(xintercept=phase.list[[1]], col="brown")+
      geom_vline(xintercept=phase.list[[2]], col="red")+
      geom_vline(xintercept=phase.list[[3]], col="green")+
      geom_vline(xintercept=phase.list[[4]], col="blue")+
      xlim(0,7e+6)+
      ggtitle(hist.cycle)+
      theme_jetka(base_size = 10)
    print(hist)
  }
  dane.cc <- dane

  colnames(dane.cc)[which(colnames(dane.cc)=="Intensity_IntegratedIntensity_DAPI")] <- 
    "DNA_signal"
  
  dane.cc$phase <- cc.assign(dane.cc, 
                             cc.borders = phase.list, 
                             DAPI.column = "DNA_signal")
  

  # in case you need a boxplots comparing separate cycle phases AND all phases use:
  # list.dane <- list(cc=dane.cc, all=dane.cc)
  # list.dane[[2]]$phase <- "all"
  # dane.cc.rbind <- rbind(list.dane[[1]], list.dane[[2]])
  
  # multiplying data of time=0 to get stim=0 in every time
  stim0 <- dane.cc[dane.cc$stimulation.1.1 == 0, ]
  list.stim0 <- list(dane.cc, stim0, stim0)
  
  list.stim0[[2]]$stimulation.1.1 <- 1
  list.stim0[[3]]$stimulation.1.1 <- 0.1
  
  dane.cells <- do.call("rbind", list.stim0)
  
  ### preparing data to treat it as one well=one observation ###
  dane.cells <- dane.cells[, c("well.name", "Intensity_IntegratedIntensity_Alexa",
                               "DNA_signal",
                               "Intensity_MeanIntensity_Alexa", "phase", "time.1.1", 
                               "stimulation.1.1")]
  colnames(dane.cells) <- c("well.name", "Integrated_Alexa",
                            "DNA_signal",
                            "Mean_Alexa", "phase", "time.1.1", 
                            "stimulation.1.1")
  phases <- factor(c("G1", "S", "G2/M"), levels = c("G1", "S", "G2/M"))
  wells <- unique(dane.cells$well.name)
  dane.means <- data.frame()
  
  for(well in wells){
    for(ph in phases){
      dane.subset <- dane.cells[dane.cells$well.name == well & 
                                  dane.cells$phase == ph, ]
      for (stim in unique(dane.subset$stimulation.1.1)){
        dane.rbind <- dane.subset[dane.subset$stimulation.1.1 == stim, ]
        dane.rbind$Average_Integrated <- mean(dane.rbind$Integrated_Alexa)
        dane.rbind$Average_Mean <- mean(dane.rbind$Mean_Alexa)
        dane.means <- rbind(dane.means, dane.rbind[1, ])
      }
    }
  }
  cell.names <- c("well.name", 
                  "stimulation.1.1", 
                  "time.1.1", 
                  "Mean_Alexa",
                  "Integrated_Alexa",
                  "DNA_signal",
                  "phase") 
  well.names <- c(cell.names[1:(length(cell.names)-4)], 
                  "Average_Integrated", 
                  "Average_Mean",
                  cell.names[length(cell.names)]) 
  
  if (to.wells==0){
    dane.phase <- dane.cells[, cell.names]
    
  } else {
    dane.phase <- dane.means[, well.names]
    
  }
  
  #data subsetting in regards to CC
  dane.phase$phase <- factor(dane.phase$phase, levels = c("G1","S","G2/M") )
  if (which.phase.1 == 0){
    return(dane.phase) } else { 
      dane.data <- subset(dane.phase, dane.phase$phase == phases[[which.phase.1]] | 
                            dane.phase$phase == phases[[which.phase.2]])
      return(dane.data)
    }
}
# 
# test:
# data.plot.31 <- prepare_data(input.folder = "PT31/R/Image",
#                              phase.list=list("out1"=1.6e6, "G1"=2.6e6, 
#                                              "S"=3.3e6, "G2/M"=4.3e6),
#                              which.phase.1 = 2, which.phase.2 = 2,
#                              to.wells = 0)
# phase.list for PT29: 
# phase.list <- list("out1"=1.5e6, "G1"=2.4e6, "S"=3.1e6, "G2/M"=4.2e6)
# phase.list for PT31: 
# phase.list <- list("out1"=0.9e6, "G1"=1.5e6, "S"=1.95e6, "G2/M"=2.8e6)
# phase.list for PT31-image:
# phase.list=list("out1"=1.6e6, "G1"=2.6e6, "S"=3.3e6, "G2/M"=4.3e6),
# phase.list for PT29-cnc, version 1:
# phase.list <- list("out1"=0.85e6, "G1"=1.50e6, "S"=1.85e6, "G2/M"=2.6e6)
# phase.list for PT29-cnc, version 2:
# phase.list <- list("out1"=0.85e6, "G1"=1.45e6, "S"=1.95e6, "G2/M"=2.4e6)