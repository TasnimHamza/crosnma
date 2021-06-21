#' Run NMA or NMR model to cross-synthesis in NMA and NMR for dichotomous outcomes
#' @description Takes jags model from an object produced by \code{crosnma.model} and runs model using \code{jags}.
#'
#' @param model A \code{crosnmaModel} object produced by running \code{crosnma.model}.
#' @param n.adapt Number of adaptations for the mcmc chains.
#' @param n.burnin Number of burnin iterations for the mcmc chains.
#' @param n.iter Number of iterations for the mcmc chains.
#' @param thin Number of thinning for the mcmc chains. Default is 1.
#' @param n.chains Number of mcmc chains. Default is 2.
#' @param quiet A logical. If TRUE, the warning message will not be displayed
#' See \code{\link{jags.model}} for more info.
#'
#' @return \code{crosnma.run} returns an object of class \code{crosrun} which is a list containing the following components:
#' @return \code{samples}  The MCMC samples produced by running the BUGS model.
#' @return \code{model}  The \code{crosnmaModel} object obtained from \code{crosnma.model} which was used to run \code{jags}.
#' @return \code{trt.key}  A Table of the treatments and its mapped integer number (used in JAGS model).
#' @examples
#' # An example from participant-level data and study-level data.
#' # data
#' data(prt.data)
#' data(std.data)
#'  #=========================#
#'   # Create a jags model  #
#'  #=========================#
#'  # We conduct a network meta-analysis assuming a random effect model.
#'  # The data comes from randomised-controlled trials and non-randomised studies. They will be combined naively.
#'  # The data has 2 different formats: individual participant data (prt.data) and study-level data (std.data).
#' mod <- crosnma.model(prt.data=prt.data,
#'                   std.data=std.data,
#'                   trt=c('trt','trt'),
#'                   study=c('study','study'),
#'                   outcome=c('outcome','outcome'),
#'                   n='n',
#'                   design=c('design','design'),
#'                   reference='A',
#'                   trt.effect='random',
#'                   covariate = NULL,
#'                   method.bias='naive'
#'                    )
#'  #=========================#
#'     # Fit jags model  #
#'  #=========================#
#' fit <- crosnma.run(model=mod,
#'                 n.adapt = 20,
#'                 n.iter=50,
#'                 thin=1,
#'                 n.chains=3)
#'
#'  #=========================#
#'    # Display the output   #
#'  #=========================#
#' summary(fit)
#' plot(fit)
#'
#'
#' @seealso \code{\link{crosnma.model}},\code{\link{jags.model}}
#' @export

crosnma.run <- function(model,
                     n.adapt = 1000,
                     n.burnin = floor(n.iter / 2),
                     n.iter,
                     thin=1,
                     n.chains=2,
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
