devtools::install_github("TasnimHamza/crosnma",force = TRUE)
library(crosnma)
library(rjags)
load.module('mix') # needed for adjust2 method, library(rjags) might be needed

#-------- data --------#
myprt.data <- read.csv("~/Google Drive/GenericModelNMA/GenericNMAanalysis/IPD for analysis") # participant data
mystd.data <- read.csv("~/Google Drive/GenericModelNMA/GenericNMAanalysis/AD for analysis") # aggregate data

#-------- MCMC settings --------#
n.adapt = 2000
n.iter=10000
n.burnin = 4000
thin=1
n.chains=2
#-------- 3. Naive NMA --------#
# jags model: code+data
mod3 <- crosnma.model(prt.data=myprt.data,
                      std.data=mystd.data,
                      trt=c('TRT01A','treat'),
                      study=c('STUDYID','study'),
                      outcome=c('RELAPSE2year','r'),
                      n='n',
                      design=c('design','design'),
                      reference='Placebo',
                      trt.effect='random',
                      method.bias = 'naive',
                      covariate=list('AGE','age')
)
source("R/crosnma.codeCLEAN.R")
mod3$jagsmodel<-crosnma.codeCLEAN(ipd = F,
                                   ad = T,
                                   trt.effect='random',
                                   prior.tau.trt="dnorm(0,0.01)T(0,)",
                                   # -------- meta-regression
                                   split.regcoef =F,
                                   covariate=list(c('AGE','age')),

                                   reg0.effect='random',
                                   regb.effect='random',
                                   regw.effect='random',

                                   prior.tau.reg0="dnorm(0,0.01)T(0,)",
                                   prior.tau.regb=NULL,
                                   prior.tau.regw=NULL,
                                   # --------  bias adjustment
                                   bias.effect=NULL,
                                   bias.type=NULL,
                                   bias.covariate=NULL,

                                   prior.tau.gamma=NULL,

                                   prior.pi.high.rct=NULL,
                                   prior.pi.low.rct=NULL,
                                   prior.pi.high.nrs=NULL,
                                   prior.pi.low.nrs=NULL,

                                   method.bias = NULL,
                                   d.prior.nrs=NULL # required when method.bias='prior'
                                   )
# mod3$data$t.ad <- matrix(c(3,4,3,4),2,2,byrow = TRUE)

cat(mod3$jagsmodel)
# run jags
jagsfit3 <- crosnma.run(model=mod3,
                        n.adapt = n.adapt,
                        n.iter=n.iter,
                        n.burnin = n.burnin,
                        thin=thin,
                        n.chains=n.chains,
                        monitor=c('LOR'))

# output
summary(jagsfit3,expo=T)
coda::traceplot(jagsfit3$samples)


