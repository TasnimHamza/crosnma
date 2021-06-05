# 1. check if it works only for IPD: YES
devtools::install_github("TasnimHamza/crosnma")
library(crosnma)
# jags model: code+data
mod1 <- crosnma.model(prt.data=prt.data,
                   std.data=std.data,
                   trt=c('trt','trt'),
                   study=c('study','study'),
                   outcome=c('outcome','outcome'),
                   design=c('design','design'),
                   n='n',
                   covariate = list(c('age'),'age'),#list(c('AGE','SEX','EDSSBL'),c('age','sex','edss')),
                   reference='A',
                   method.bias = 'prior',
                   run.nrs=list(
                                n.adapt = 20,
                                n.iter=50,
                                n.burnin = 3,
                                thin=1,
                                n.chains=2)
)



# fit jags
fit1 <- crosnma.run(model=mod1,
                 n.adapt = 20,
                 n.iter=1000,
                 thin=1,
                 n.chains=3,
                 n.burnin = 400)
summary(fit1,expo = TRUE)

autocorr.plot(fit1$samples)
cumuplot(fit1$samples)
densplot(fit1$samples)



geweke.plot(fit1$samples)
plot(fit1$samples)
traceplot(fit1$samples)


library(netmeta)

data(Senn2013)

# Generation of an object of class 'netmeta' with reference
# treatment 'plac'
#

mod1$data$t.ipd
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = Senn2013, sm = "MD")

# Network graph with default settings
#
netgraph(net1)


