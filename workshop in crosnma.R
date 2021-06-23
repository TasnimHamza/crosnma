#  Workshop in crosnma
#  Authors: Tasnim Hamza and Georgia Salanti

#-------- load the libray --------#
# install.packages("devtools") # install the package if you didn't before
devtools::install_github("TasnimHamza/crosnma",force = TRUE)
library(crosnma)
load.module('mix') # needed for adjust2 method, library(rjags) might be needed

#-------- data --------#
head(prt.data) # participant data
head(std.data) # aggregate data

#-------- 1. Naive NMA --------#
# jags model: code+data
mod1 <- crosnma.model(prt.data=prt.data,
                      std.data=std.data,
                      trt=c('trt','trt'),
                      study=c('study','study'),
                      outcome=c('outcome','outcome'),
                      n='n',
                      design=c('design','design'),
                      trt.effect='random',
                      reference='A',
                      method.bias = 'naive'
)

# run jags
jagsfit1 <- crosnma.run(model=mod1,
                        n.adapt = 20,
                        n.iter=100,
                        n.burnin = 40,
                        thin=1,
                        n.chains=2)

# output
summary(jagsfit1,expo=T)
coda::traceplot(jagsfit1$samples)

mean(prt.data$sex[prt.data$study==3])
#-------- 2. Naive NMR --------#
# jags model: code+data
mod2 <- crosnma.model(prt.data=prt.data,
                      std.data=std.data,
                      trt=c('trt','trt'),
                      study=c('study','study'),
                      outcome=c('outcome','outcome'),
                      n='n',
                      design=c('design','design'),
                      reference='A',
                      trt.effect='random',
                      #----------  meta-regression ----------
                      covariate = list(c('age'),c('age')),
                      split.regcoef = T,
                      #---------- bias adjustment ----------
                      method.bias='naive'
)
# run jags
jagsfit2 <- crosnma.run(model=mod2,
                        n.adapt = 500,
                        n.iter=10000,
                        n.burnin = 4000,
                        thin=1,
                        n.chains=3)
# output
summary(jagsfit2,expo=T)
coda::traceplot(jagsfit2$samples)

#-------- 3. prior NMA --------#
# jags model: code+data
mod3 <- crosnma.model(prt.data=prt.data,
                      std.data=std.data,
                      trt=c('trt','trt'),
                      study=c('study','study'),
                      outcome=c('outcome','outcome'),
                      n='n',
                      design=c('design','design'),
                      reference='A',
                      trt.effect='random',
                      #----------  meta-regression ----------
                      covariate = list(c('age'),c('age')),
                      split.regcoef = F,
                      #---------- bias adjustment ----------
                      method.bias='prior',
                      run.nrs=list(n.adapt = 10,
                                   n.iter=5000,
                                   n.burnin = 400,
                                   thin=1,
                                   n.chains=2))

cat(mod3$jagsmodel)# run jags
jagsfit3 <- crosnma.run(model=mod3,
                        n.adapt = 200,
                        n.iter=10000,
                        n.burnin = 4000,
                        thin=1,
                        n.chains=2)

# output
summary(jagsfit3,expo=T)
coda::traceplot(jagsfit3$samples)


#-------- 4. adjust1 NMA --------#
# jags model: code+data
mod4 <- crosnma.model(prt.data=prt.data,
                      std.data=std.data,
                      trt=c('trt','trt'),
                      study=c('study','study'),
                      outcome=c('outcome','outcome'),
                      n='n',
                      design=c('design','design'),
                      reference='A',
                      trt.effect='random',
                      #----------  meta-regression ----------
                      covariate = list(c('sex'),c('sex')),
                      split.regcoef = T,
                      #---------- bias adjustment ----------
                      method.bias='adjust1',
                      bias=c('bias','bias'),
                      bias.type='add',
                      bias.effect='common'
)
# run jags
jagsfit4 <- crosnma.run(model=mod4,
                        n.adapt = 500,
                        n.iter=10000,
                        n.burnin = 4000,
                        thin=1,
                        n.chains=3)
# output
summary(jagsfit4,expo=T)
coda::traceplot(jagsfit4$samples)
#-------- 5. adjust2 NMA --------#
# jags model: code+data
mod5 <- crosnma.model(prt.data=prt.data,
                      std.data=std.data,
                      trt=c('trt','trt'),
                      study=c('study','study'),
                      outcome=c('outcome','outcome'),
                      n='n',
                      design=c('design','design'),
                      reference='A',
                      trt.effect='random',
                      #---------- bias adjustment ----------
                      method.bias='adjust2',
                      bias=c('bias','bias'),
                      bias.type='add',
                      bias.effect='common',
                      bias.covariate = c('year','year'),#
                      #---------- assign a prior ----------
                      prior=list(tau.trt='dunif(0,3)',
                                 pi.high.rct='dbeta(5,1)',
                                 pi.low.rct='dbeta(1,20)',
                                 pi.high.nrs='dbeta(30,1)',
                                 pi.low.nrs='dbeta(1,2)'
                      )
)

# run jags
jagsfit5 <- crosnma.run(model=mod5,
                        n.adapt = 20,
                        n.iter=100,
                        n.burnin = 40,
                        thin=1,
                        n.chains=2)
# output
summary(jagsfit5,expo=T)
coda::traceplot(jagsfit5$samples)







