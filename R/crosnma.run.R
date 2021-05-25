#!!! Think of what to monitor
#' Run Generic NMA or NMR model for dichotomous outcomes
#' @description Takes jags code from an object produced by \code{crosnma.model} and runs model using \code{jags}.
#'
#' @param model A \code{crosnmaModel} object produced by running \code{crosnma.model}.
#' @param n.adapt Number of adaptations for the mcmc chains.
#' @param n.burnin Number of burnin iterations for the mcmc chains.
#' @param n.iter Number of iterations for the mcmc chains.
#' @param thin Thinning factor for the mcmc chains. Default is 1.
#' @param n.chains Number of mcmc chains. Default is 3.
#' See \code{\link{jags.model}} for more info.
#'
#' @return \code{crosnma.run} returns an object of class \code{crosnmaRun} which is a list containing the following components:
#' @return \code{samples} - The MCMC samples produced by running the BUGS model.
#' @return \code{model} - The \code{BUGSnetModel} object obtained from \code{nma.model} which was used to run \code{jags}.
#' @return \code{trt.key} - Treatments mapped to numbers, used to run BUGS code.
#' @export
#' @seealso \code{\link{crosnma.model}},\code{\link{jags.model}}
#' @examples


crosnma.run <- function(model,
                     n.adapt = 1000,
                     n.burnin = floor(n.iter / 2),
                     n.iter,
                     thin=1,
                     n.chains=3,
                     quiet=TRUE
                     # inits = "DEFAULT"
){


  if(class(model) != "crosnmaModel")
    stop("\'model\' must be a valid crosnmaModel object created using the crosnma.model function.")


  seeds <- sample(.Machine$integer.max, n.chains, replace = FALSE)
  inits <- list()
  for (i in 1:n.chains)
    inits[[i]] <- list(.RNG.seed = seeds[i], .RNG.name = "base::Mersenne-Twister")


  jagsfit <- jags.model(textConnection(model$jags),        #Create a connection so JAGS can access the variables
                          model$data,
                          n.chains=n.chains,
                          n.adapt=n.adapt,
                          inits = inits,
                          quiet=quiet)

  # runjags.options(silent.jags=TRUE, silent.runjags=TRUE)
  if(n.burnin!=0) jagsburnin <- update(jagsfit, n.iter=n.burnin)

  # ---- monitor
  # basics
  make.monitor <- "d"
  if(model$trt.effect=="random") make.monitor <- c(make.monitor, "tau")
  # meta-regression
  if(!is.null(model$covariate)){
    if(model$split.regcoef){
      make.monitor.reg <- c()
      for (i in 1:length(model$covariate[[1]])) {
        make.monitor.reg0 <- c(paste0("bb_",i),paste0("bw_",i))
        make.monitor.reg <- c(make.monitor.reg,make.monitor.reg0)
      }
      if(model$regb.effect=='random'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg1 <- paste0("tau.bb_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg1)
        }
      }
      if(model$regw.effect=='random'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg2 <- paste0("tau.bw_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg2)
        }
      }
    }else{
      make.monitor.reg <- c()
      for (i in 1:length(model$covariate[[1]])) {
        make.monitor.reg0 <- paste0("b_",i)
        make.monitor.reg <- c(make.monitor.reg,make.monitor.reg0)
      }
      if(model$regb.effect=='random'|model$regw.effect=='random'){
        for (i in 1:length(model$covariate[[1]])) {
          make.monitor.reg1 <- paste0("tau.b_",i)
          make.monitor.reg <- c(make.monitor.reg,make.monitor.reg1)
        }
      }
    }
    make.monitor <- c(make.monitor, make.monitor.reg)
  }

  if(!is.null(model$method.bias)){
    make.monitor.bias <- c()
    if(model$method.bias=='adjust1'){
      if(model$bias.type=='both'){
        make.monitor.bias <- c("g1","g2")
        if(model$bias.effect=='random') make.monitor.bias <- c(make.monitor.bias,"tau.gamma1","tau.gamma2")
      }else if(model$bias.type%in%c('add','mult')){
        make.monitor.bias <- c("g")
        if(model$bias.effect=='random') make.monitor.bias <- c(make.monitor.bias,"tau.gamma")
      }
    }
    if(model$method.bias=='adjust2'){
      make.monitor.bias <- c("g")
      if(model$bias.effect=='random') make.monitor.bias <- c(make.monitor.bias,"tau.gamma")

    }
    make.monitor <- c(make.monitor, make.monitor.bias)
  }


  jagssamples <- coda.samples(jagsfit,
                              variable.names=make.monitor,
                              n.iter=n.iter,
                              thin=thin)

  crosrun <- structure(list(samples=jagssamples,
                         model=model,
                         "trt.key"=model$trt.key),
                    class = "crosnma")
  return(crosrun)
}
