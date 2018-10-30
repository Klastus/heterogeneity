## Script for plotting for Japan poster ##

if (Sys.info()["sysname"]=="Darwin"){
  # PT.MAC
  setwd("/Users/piotrt/Documents/IPPT PT/JAK-STAT model/cellcycle")} else 
    if (Sys.info()["sysname"]=="Linux"){
    # PT.Linux:
    setwd('/home/piotrek/Documents/heterogeneity_model/heterogeneity') } else {
      # PT.Win:
      setwd('Y:/PiotrT/cellcyclemodel')
    }
#### loading packages ####
library(grid)
#### loading requires scripts ####
core.path <- "/home/piotrek/Documents/heterogeneity_model/heterogeneity/scripts"
source(paste(core.path, 
             "/R/data_preparing/prepare_data.R", sep=''))
source(paste(core.path, 
             "/R/data_preparing/prepare_model.R", sep=''))

#####  loading experimental data #####
# for poster, only data from wells was used, so below 7 lines is not nessesary

PT31.all.cells <- prepare_data(input.folder = "/home/piotrek/Documents/PT/microscope_input/PT31-image",
                             phase.list=list("out1"=1.6e6, "G1"=2.6e6, "S"=3.3e6, "G2/M"=4.3e6),
                             which.phase.1 = 1, which.phase.2 = 1,
                             to.wells = 0)

PT31.all.cells <- PT31.all.cells[!grepl("A|B", PT31.all.cells$well.name), ]

#### loading model data ####
factor.for.model <- 1/280
CV.list <- list("real", 0.1, 0.3, 0.7, 0.9)
list.abundance <- list(1414)

start <- Sys.time()
for (abundance in list.abundance){
  if (abundance > 1414){
    replace <- TRUE
  } else {replace <- FALSE}
  PT31.sampled.histo <- rbind(sample_n(PT31.all.cells[PT31.all.cells$phase=="G1" & 
                                            PT31.all.cells$stimulation.1.1 == 0.1 &
                                            PT31.all.cells$time.1.1 == 15, ], abundance,
                                       replace = replace),
                        sample_n(PT31.all.cells[PT31.all.cells$phase=="G1" & 
                                                  PT31.all.cells$stimulation.1.1 == 1 &
                                                  PT31.all.cells$time.1.1 == 15, ], abundance, 
                                 replace = replace))
  
  plots <- list()
  for (CV in CV.list){
    file.to.G1.model <- paste("/home/piotrek/Documents/heterogeneity_model/input/sampled_traj/", 
                          abundance, "/",
                          CV,
                          "/sampled_traj.csv",
                          sep='')

model.G1.cells <- prepare_model_par(input.file = file.to.G1.model,
                                dane=PT31.all.cells[PT31.all.cells$phase=="G1", ], 
                                factor=factor.for.model,
                                no.cores = 8,
                                phase= "G1")
# head(model.G1.cells)
model.all.cells <- model.G1.cells
#### exact plotting cells ####

# subsetting of data, to merge and combine values from repetitional wells
# maybe boxplots instead of mean and sd ?:
title <- paste("PT31 cells, CV=", CV, "_", abundance, sep= '')
plots[["box_lines"]][[toString(CV)]] <- ggplot()+
  geom_line(data=model.all.cells[model.all.cells$stimulation.1.1 !=0, ], 
            aes(x = time,
                y = log10(value),
                group=sigma), 
            size = 0.8,
            alpha = 0.05)+
  geom_violin(data = PT31.all.cells[PT31.all.cells$phase=="G1" & 
                                       PT31.all.cells$stimulation.1.1 != 0 , ],
               aes(x = time.1.1, y = log10(Mean_Alexa), group = time.1.1),
               fill = "darkgoldenrod1", color= "goldenrod4")+
  theme_jetka(base_size =7)+
  facet_grid(phase ~ stimulation.1.1)+
  ylim(1.2, 3)+
  ggtitle(title)



# histograms somehow:

title <- paste("PT31 cells, CV=", CV, "_", abundance,  "_15 min", sep= '')
plots[["histo"]][[toString(CV)]] <- ggplot()+
  geom_histogram(data = PT31.sampled.histo,
                 aes(x = log10(Mean_Alexa)), 
                 colour = "goldenrod4", 
                 bins = 50,
                 fill = "darkgoldenrod1")+
  geom_histogram(data = model.all.cells[model.all.cells$time==15 &
                                          model.all.cells$stimulation.1.1 != 0, ],
                 aes(x = log10(value)), 
                 colour = "blue", 
                 fill = "darkblue", 
                 bins = 50, 
                 alpha = 0.5)+
  theme_jetka(base_size =7)+
  facet_grid(phase ~ stimulation.1.1)+
  xlim(1.5, 3.5)+
  ggtitle(title)
plots[["histo"]][[toString(CV)]]
  }
  plot.path <- paste(getwd(), "/plots/", abundance, "/", sep = '')
  if(!dir.exists(plot.path)){
    dir.create(plot.path)
  }
  pdf(file = paste(plot.path, "box.pdf", sep =''), width = 15, height = 7.5)
  print(plots[["box_lines"]])
  dev.off()
  
  pdf(file = paste(plot.path, "histo.pdf", sep =''), width = 15, height = 7.5)
  print(plots[["histo"]])
  dev.off()
}
end <- Sys.time()
end-start

# CV=1_sampled_1.5e3
