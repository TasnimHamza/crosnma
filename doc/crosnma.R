## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(crosnma)

## -----------------------------------------------------------------------------
head(prt.data)

## -----------------------------------------------------------------------------
head(std.data)

## -----------------------------------------------------------------------------
netplot(prt.data,std.data)

## -----------------------------------------------------------------------------
knitr::kable(ns.tab(prt.data,std.data))

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
# run jags
jagsfit1 <- crosnma.run(model=mod1,
              n.adapt = 20,
              n.iter=50,
              thin=1,
              n.chains=3)

## -----------------------------------------------------------------------------
knitr::kable(summary(jagsfit1,expo=T,))

## -----------------------------------------------------------------------------
# to be added

## -----------------------------------------------------------------------------
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
                   covariate = list(c('age','sex'),c('age','sex')),
                   #---------- bias adjustment ---------- 
                   method.bias='prior',
                   run.nrs=list(n.adapt = 10,
                                n.iter=10,
                                n.burnin = 5,
                                thin=1,
                                n.chains=1))

## -----------------------------------------------------------------------------
# run jags
jagsfit2 <- crosnma.run(model=mod2,
              n.adapt = 20,
              n.iter=50,
              thin=1,
              n.chains=3)

## -----------------------------------------------------------------------------
knitr::kable(summary(jagsfit2))

## -----------------------------------------------------------------------------
#plot(jagsfit1)

## -----------------------------------------------------------------------------
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
                   covariate = list(c('age','sex'),c('age','sex')),
                   #---------- bias adjustment ---------- 
                   method.bias='adjust1',
                   bias=c('bias','bias'), 
                   bias.type='add',
                   bias.effect='common',
                   #---------- assign a prior ---------- 
                   prior=list(tau.trt='dunif(0,3)',
                              pi.high.rct='dbeta(5,1)',
                              pi.low.rct='dbeta(1,20)',
                              pi.high.nrs='dbeta(30,1)',
                              pi.low.nrs='dbeta(1,2)'
                                  )
                  )

## -----------------------------------------------------------------------------
# run jags
jagsfit3 <- crosnma.run(model=mod3,
              n.adapt = 20,
              n.iter=50,
              thin=1,
              n.chains=3)

## -----------------------------------------------------------------------------
knitr::kable(summary(jagsfit3))

## -----------------------------------------------------------------------------
# to be added

## -----------------------------------------------------------------------------
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
                   covariate = list(c('age','sex'),c('age','sex')),
                   #---------- bias adjustment ----------
                   method.bias='adjust1',
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

## -----------------------------------------------------------------------------
# run jags
jagsfit4 <- crosnma.run(model=mod4,
              n.adapt = 20,
              n.iter=50,
              thin=1,
              n.chains=3)

## -----------------------------------------------------------------------------
knitr::kable(summary(jagsfit4))

## -----------------------------------------------------------------------------
# to be added

