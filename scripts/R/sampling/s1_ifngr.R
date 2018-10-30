#### aim: to have n pairs of s1 and ifngr number. 
#### Later you will have to rescale it somehow :/
# core.path <- "/Users/piotrt/Documents/IPPT_PT/"
core.path <- "/home/piotrek/Documents/heterogeneity_model/heterogeneity/scripts"
# setwd(paste(core.path, "/scripts/", sep = ''))

source(paste(core.path, "/R/basic_scripts/normalize_data.R", sep=''))
source(paste(core.path, "/R/basic_scripts/theme_jetka.R", sep=''))
source(paste(core.path, 
             "/R/sampling/STAT1_distribution.R", sep=''))
source(paste(core.path, 
             "/R/sampling/IFNGR1_distribution.R", sep=''))
source(paste(core.path, 
             "/R/sampling/sample_protein.R", sep=''))

#### changing the number of sampled numbers:
origin.stat1 <- s1.fluo()
origin.ifngr <- ifngr.fluo()

factor.stat1 <- 152490/mean(origin.stat1)
factor.ifngr <- 140/mean(origin.ifngr)


list.CV <- list("real",  0.1, 0.3, 0.7, 0.9)

CV.ifngr <- TRUE
CV.stat1 <- TRUE

list.abundance <- list(1414)

for(CV in list.CV) {
  for (abundance in list.abundance){

    if(!CV.stat1){
      stat1 <- sample.protein(origin.stat1, 
                              CV.given = CV.stat1,
                              n = abundance)
    } else {
      stat1 <- sample.protein(origin.stat1, 
                                    CV.given = CV,
                                    n = abundance)
    }
    
    if(!CV.ifngr){
      ifngr <- sample.protein(origin.ifngr, 
                              CV.given = CV.ifngr,
                              n = abundance)
    } else {
      ifngr <- sample.protein(origin.ifngr, 
                                    CV.given = CV,
                                    n = abundance)
    }
    
sampled <- data.frame("stat1" = stat1$sample,
                     "ifngr" = ifngr$sample)
## ##
### what if those factors will be as mean of sampled distributions? 
## it would give same results as taking mean of origin distributions?
## ##

sampled$factor.stat1 <- sampled$stat1*factor.stat1
sampled$factor.ifngr <- sampled$ifngr*factor.ifngr

protein.folder.name <- names(which(list("STAT1" = CV.stat1,
                            "IFNGR" = CV.ifngr) ==FALSE))

if(identical(protein.folder.name, character(0))){
  protein.folder.name <- "IFNGR_STAT1"
}

path.to.save <- paste(core.path, 
                      "/MATLAB/models/JAKSTATdensityver1/sampled_protein/", 
                      protein.folder.name, "/",
                      abundance, "/",
                      CV,
                      sep = '')

if(!dir.exists(path.to.save)){
  dir.create(path.to.save, recursive = TRUE)
}
write.csv(sampled, file = paste(path.to.save, "/sampled_s1_ifngr.csv", sep =''))

  }
}
# helping plotting
# ggplot(data = data.frame(integrated=origin.stat1)) + geom_histogram(aes(integrated), bins=50) + xlim(-1.5e6, 1.5e7)
# ggplot(data = sampled) + geom_histogram(aes(stat1), bins=50) + xlim(-1.5e6, 1.5e7)
# ggplot(data = sampled) + geom_histogram(aes(factor.stat1), bins=50) 
# ggplot(data = data.frame(integrated=origin.stat1)) + geom_histogram(aes(log(integrated)), bins=50) + xlim(10, 18)
# ggplot(data = sampled) + geom_histogram(aes(log(stat1)), bins=50) + xlim(10, 18)
