### Script for preparing model data with 5 sigma points for plotting each sigma point case ###
# #
if (Sys.info()["sysname"]=="Darwin"){
  # PT.MAC
  setwd("/Users/piotrt/Documents/IPPT PT/JAK-STAT model/cellcycle")} else {
    # PT.Win:
    setwd('Y:/PiotrT/cellcyclemodel')
  }

#Model preparing script:
prepare_model_sigma <- function(input.file = 'compenG2_data_model_sigmapoints.csv', 
                                stim.level = c(0, 0.1, 1), sigma = 5,
                                data=pS1.means, factor=1/400){
  # prepare_model function prepares .csv model data for plotting, that is:
  # inverts the table and names columns as variable and rows as time observations, 
  # melts the data, what is required for ggplotting
  # adjusts the values to be comparable to experimental data
  #
  # Arg:
  # stim.level- all levels of stimulations, has to be provided according to MatLab
  # data- the data according to which the background is calculated
  # factor- a factor by which the model is adjusted
  # sigma- number of sigma points used
  
  # model data loading #
  model<-read.csv(paste(getwd(), "MATLAB/models/JAKSTATdensityver0/Output", input.file, 
                        sep='/'),sep=',',header = 0)
  
  ####model data processing ####
  rep.times <- function (times){
    rep(paste("var", (1:(n.comb/sigma)), sep=''), times=times)
  }
  n.comb <- (ncol(JS) / l - 2)
  l <- length(stim.level)
  JS <- model
  JS <- data.frame(t(JS))
  

  colnames(JS) <- rep(c(rep.times(times=sigma), "stimulation.1.1","time"), times=l)
  row.names(JS) <- (1:nrow(JS))
  
  #melting model data
  n <- (n.comb+2)
  JS.melt <- data.frame()
  JS.a <- data.frame()
  
  for(x in (0:(l - 1))){
    JS.a <- JS[, c((x * n + 1):((x + 1) * n))]
    for (si in 0:(sigma-1)){
      n.var <- (n-2)/sigma
      JS.c <- JS.a[, (si * n.var + 1):(n.var*(si+1))]
      colnames(JS.c) <- rep.times(times=1)
      JS.c$stimulation.1.1 <- JS.a$stimulation.1.1
      JS.c$time <- JS.a$time
      JS.c$sigma <- si+1
      JS.b <- melt(JS.c, id.vars = c("time", "stimulation.1.1", "sigma"), 
                   value.name = "value", variable.name = "variable")
      JS.melt <- rbind(JS.melt, JS.b)
    }
  }
  
  #adjusting function
  adjust <- function (data=pS1.means, factor){
    background<-mean(data[data$time.1.1 ==  0, ]$average_mean)
    JS.14 <- JS.melt[JS.melt$variable == "var14", ]
    JS.14$value<-JS.14$value * factor + background
    return(JS.14)
  }
  JS.14n<-adjust(factor=factor)
}