#### function to sample data from given distribution but changed CV 
#### aim: to have n pairs of s1 and ifngr numbers, with given CV

#### ####

sample.protein <- function(data.fluo, CV.given, n){
  m <- mean(data.fluo)
  v <- var(data.fluo)
  sd <- v^0.5
  CV.real <- sd/m 
  
  if (CV.given == FALSE){
    return(list("sample" = rep(m, times = n),
                "CV.real" = 0,
                "CV.given" = 0, 
                "CV.sampled" = 0))
  }


  if(CV.given == "real"){
    CV.given <- CV.real
  }
  miu <- log(m^2/(v+m^2)^0.5)
  sigma <- log(CV.given^2+1)^0.5 ## most important formula!
  # lognormal CV:
  protein.sample <- rlnorm(n = n, 
                           meanlog = miu, 
                           sdlog = sigma)
  CV.sampled <- sd(protein.sample)/mean(protein.sample)
  return(list("sample" = protein.sample,
              "CV.real" = CV.real,
              "CV.given" = CV.given, 
              "CV.sampled" = CV.sampled))
}

# will be used in the future?:
# we have miu and log.CV, want to check the raw.CV. This is the formula:
# raw.CV <- (exp(miu*(2+log.CV))*exp((log.CV*miu)^2-1))^0.5/exp((log.CV^2*miu+2)*miu/2)
# library(fitdistrplus)
# s1.distrib <- fitdist(s1$integrated, distr = "lnorm")