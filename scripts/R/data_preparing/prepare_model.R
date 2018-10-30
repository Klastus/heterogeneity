### Script for preparing model data both with and w/o sigma points 
#   and next for plotting each sigma point case separately ###
# #

#Model preparing script:
prepare_model <- function(input.file,
                          stim.level = c(0, 0.1, 1),
                          dane,
                          phase,
                          factor,
                          bck=0,
                          no.cores,
                          sigma,
                          to.log=1){
  # prepare_model function prepares .csv model data for plotting, that is:
  # -inverts the table and names columns as variable and rows as time observations, 
  # -melts the data, what is required for ggplotting
  # -adjusts the values to be comparable to experimental data
  #
  # Arg:
  # stim.level- all levels of stimulations, has to be provided according to MatLab
  # data- the data according to which the background is calculated
  # factor- a factor by which the model is adjusted
  # sigma- number of sigma points used
  
  # model data loading #
  model <- read.csv(input.file, sep=',', header = 0)
  
  #### model data processing ####
  JS <- model
  JS <- data.frame(t(JS))
  l <- length(stim.level)
  m.size <- dim(JS)[2]
  
  if(sigma==-1){
    sigma <- ((m.size/3)-2)/19
  }
  n.comb <- (ncol(JS) / l - 2)
  
  rep.times <- function (times){
    rep(paste("var", (1:(n.comb/sigma)), sep=''), times=times)
  }
  
  
  if (sigma==0){
    colnames(JS) <- rep(c(paste("var", (1:n.comb), sep=''), 
                          "stimulation.1.1", "time"), times=l)
  } else {
    colnames(JS) <- rep(c(rep.times(times=sigma), "stimulation.1.1","time"), times=l)
  }  
  
  row.names(JS) <- (1:nrow(JS))
  
  #melting model data
  n <- (n.comb+2)
  JS.melt <- data.frame()
  JS.a <- data.frame()
  if (sigma==0){
    for(x in (0:(l - 1))){
      JS.a <- JS[, (x * n + 1):((x + 1) * n)]
      JS.b <- melt(JS.a, id.vars = c("time", "stimulation.1.1"), 
                   value.name = "value", variable.name = "variable")
      JS.melt <- rbind(JS.melt, JS.b)
    }
  } else {
    for(x in (0:(l - 1))){
      JS.a <- JS[, c((x * n + 1):((x + 1) * n))]
      for (si in 0:(sigma-1)){
        n.var <- (n-2)/sigma
        JS.c <- JS.a[, (si * n.var + 1):(n.var*(si+1))]
        colnames(JS.c) <- rep.times(times=1)
        JS.c$stimulation.1.1 <- JS.a$stimulation.1.1
        JS.c$time <- JS.a$time
        # JS.c$sigma <- si+1
        JS.c$sigma <- paste(JS.c$"var18"[1], "xS1\n",
                            JS.c$"var19"[1], "xIFNGR", sep = '')
        JS.b <- melt(JS.c, id.vars = c("time", "stimulation.1.1", "sigma"), 
                     value.name = "value", variable.name = "variable")
        JS.melt <- rbind(JS.melt, JS.b)
      }
    }
  }
  JS.melt$phase <- phase
  #adjusting function
  adjust <- function (data=dane, factor.tmp, bck.tmp=bck){
    if(bck.tmp==0){
      column <- colnames(data)[grepl("Mean", colnames(data))]
      if(to.log==1){
      background <- 10^(mean(log10(data[data$time.1.1 == 0, ][[column]])))
      } else {
        background <- mean(data[data$time.1.1 == 0, ][[column]])
      }
    } else {background <- bck}
    JS.14 <- JS.melt[JS.melt$variable == "var14", ]
    JS.14$value <- JS.14$value * factor.tmp + background
    return(JS.14)
  }
  JS.14n <- adjust(factor.tmp=factor)
  return(JS.14n)
}
 # testing
# model.1 <- prepare_model(input.file = 'compenG2_data_model_sigmapoints.csv',
#                          stim.level = c(0, 0.1, 1),
#                          sigma = 5, dane=data.plot.31, factor=1/400)

