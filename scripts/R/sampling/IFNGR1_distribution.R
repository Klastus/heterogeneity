#### preparation of the ifngr distribution 
#### data preparation- loading and normalization ####

ifngr.fluo <- function(){
experiment <- "input/IFNGR/PT58"
path <- paste("/home/piotrek/Documents/heterogeneity_model/",
              experiment, sep='')
normaliz <- "/raw/"

IR.all <- normalize_data(read.table(paste(path, normaliz, 
                                              "CellsFiltered488.csv", sep=''),
                                        sep = ',', header = TRUE))$data
return(IR.all[IR.all$antibody.1.1!=0, ]$Intensity_IntegratedIntensity_Alexa555)
}

#### checking, whether a model distribution is well fitted
# IR.integrated <- 
#   IR.all[IR.all$antibody.1.1!=0, ]$Intensity_IntegratedIntensity_Alexa555
# 
# IR <- data.frame(integrated = IR.integrated)
# no <- length(IR[, 1])
# IR.distrib <- fitdist(IR$integrated, distr = "lnorm")
# 
# fitted.IR <- data.frame(x=rlnorm(n = no, 
#                                  meanlog = IR.distrib$estimate[1], 
#                                  sdlog = IR.distrib$estimate[2]))
# 
# 
# ggplot(data = fitted.IR) + geom_histogram(aes((x)), bins=30) + xlim(0, 1e7)
# ggplot(data = IR) + geom_histogram(aes((integrated)), bins=30) + xlim(0, 1e7)
# 
