#### preparation of the STAT1 distribution 
## microscopic data ##
#### data preparation- loading and normalization ####

s1.fluo.long <- function(){
experiment <- "input/STAT1/PT55"
path <- paste("/home/piotrek/Documents/heterogeneity_model/",
              experiment, sep='')
normaliz <- "/raw/"

merge.compartments <- function(dye, path, normaliz, dye.rename=dye){
  df.nuc <- read.table(paste(path, normaliz, "ShrinkedNucleiMasked", 
                             '.csv', sep=''), 
                       header=TRUE, sep=",")
  df.cells <- read.table(paste(path, normaliz, "CellsFiltered", 
                               dye, '.csv', sep=''), 
                         header=TRUE, sep=",")
  df.plasm <- read.table(paste(path, normaliz, "Cytoplasm", 
                               dye, '.csv', sep=''), 
                         header=TRUE, sep=",")
  df.nuc <- normalize_data(df.nuc)$data
  df.cells <- normalize_data(df.cells)$data
  df.plasm <- normalize_data(df.plasm)$data
  
  
  columns <- c("Intensity_IntegratedIntensity_Alexa488",
               "Intensity_MeanIntensity_Alexa488",
               "Intensity_IntegratedIntensity_Alexa555",
               "Intensity_MeanIntensity_Alexa555",
               "Intensity_IntegratedIntensity_Alexa647",
               "Intensity_MeanIntensity_Alexa647",
               "Intensity_IntegratedIntensity_DAPI_nonconfo",
               "Intensity_MeanIntensity_DAPI_nonconfo")
  
  
  new.column.core <- c("Integrated_Alexa488",
                       "Integrated_Alexa555",
                       "Integrated_Alexa647",
                       "Integrated_DAPI_nonconfo",
                       "Mean_Alexa488",
                       "Mean_Alexa555",
                       "Mean_Alexa647",
                       "Mean_DAPI_nonconfo")
  df.nuc2 <- variable.subset(df.nuc, columns, 
                             new.columns = paste0(new.column.core, 
                                                  "_nucleus", dye.rename))
  df.cells2 <- variable.subset(df.cells, columns, 
                               new.columns = paste0(new.column.core, 
                                                    "_cells", dye.rename))
  df.plasm2 <- variable.subset(df.plasm, columns, 
                               new.columns = paste0(new.column.core, 
                                                    "_plasm", dye.rename))
  
  return(cbind(df.cells[, head(colnames(df.cells), 5)],
               "AreaShape_Area_nuc"=df.nuc[["AreaShape_Area"]],
               df.nuc2, 
               df.cells2,
               df.plasm2,
               df.cells[, tail(colnames(df.cells), 12)]))
}
s1.all <- merge.compartments(dye="647", 
                             path=path, 
                             normaliz = normaliz, 
                             dye.rename = '')

s1.only <- s1.all[s1.all$antibody.1.1!=0 & s1.all$antibody.1.2 == 0 &
                    s1.all$cells.1.1 == "MEF", ]
return(s1.only$Integrated_Alexa555_cells)
}

s1.fluo <- function(){
  experiment <- "input/STAT1/PT55"
  path <- paste("/home/piotrek/Documents/heterogeneity_model/",
                experiment, sep='')
  normaliz <- "/raw/"
  
  merge.compartments <- function(dye, path, normaliz, dye.rename=dye){
    df.cells <- normalize_data(read.table(paste(path, normaliz, "CellsFiltered", 
                                 dye, '.csv', sep=''), 
                           header=TRUE, sep=","))$data

    return(df.cells)
  }
  s1.all <- merge.compartments(dye="647", 
                               path=path, 
                               normaliz = normaliz, 
                               dye.rename = '')
  
  s1.only <- s1.all[s1.all$antibody.1.1!=0 & s1.all$antibody.1.2 == 0 &
                      s1.all$cells.1.1 == "MEF", ]
  return(s1.only$Intensity_IntegratedIntensity_Alexa555)
}


#### checking, whether a model distribution is well fitted
# s1 <- data.frame(integrated = s1.fluo())
# no <- length(s1[, 1])
# s1.distrib <- fitdist(s1$integrated, distr = "lnorm")
# 
# fitted.s1 <- data.frame(x=rlnorm(n = no, 
#                                  meanlog = s1.distrib$estimate[1], 
#                                  sdlog = s1.distrib$estimate[2]))
# 
# 
# ggplot(data = fitted.s1) + geom_histogram(aes((x)), bins=30) + xlim(0, 1e7)
# ggplot(data = s1) + geom_histogram(aes((integrated)), bins=30) + xlim(0, 1e7)

