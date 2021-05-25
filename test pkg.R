# 1. check if it works only for IPD: YES
devtools::load_all()
devtools::install_github("TasnimHamza/crosnma")
library(crosnma)
# jags model: code+data
mod1 <- crosnma.model(prt.data=prt.data,
                   std.data=NULL,
                   trt=c('trt','trt'),
                   study=c('study','study'),
                   outcome=c('outcome','outcome'),
                   design=c('design','design'),
                   n='n',
                   covariate = list(c('age'),'age'),#list(c('AGE','SEX','EDSSBL'),c('age','sex','edss')),
                   reference='A',
                   method.bias = 'naive'
)



# fit jags
fit1 <- crosnma.run(model=mod1,
                 n.adapt = 20,
                 n.iter=100,
                 thin=1,
                 n.chains=3)
#library(coda)
summary(fit1,expo = TRUE)
# plot(fit1)
# gsub("d.\\",fit1$trt.key[,1])
# grepl("^d", rownames(summary(fit1)))
# rownames(summary(fit1))%in%match("d")
# match("d",rownames(summary(fit1)))
# sum.nam <- rownames(summary(fit1))
# sum.fit
# sum.nam[startsWith(sum.nam,"d")]
sum.fit <- summary(fit1)


library(rjags)
data(LINE)
LINE$recompile()
LINE.out <- coda.samples(LINE, c("alpha","beta","sigma"), n.iter=1000)
summary(LINE.out)
names(summary(LINE.out))
summary(LINE.out)$quantiles

library(coda)
autocorr.plot(fit1$samples)
cumuplot(fit1$samples)
densplot(fit1$samples)



geweke.diag(fit1$samples)
geweke.plot(fit1$samples)
plot(fit1$samples)
traceplot(fit1$samples)

#
J <- 8.0
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6)


jags.data <- list("y","sd","J")
jags.params <- c("mu","sigma","theta")
jags.inits <- function(){
  list("mu"=rnorm(1),"sigma"=runif(1),"theta"=rnorm(J))
}

## You can input data in 4 ways
## 1) data as list of character
jagsfit <- jags(data=list("y","sd","J"), inits=jags.inits, jags.params,
                n.iter=10, model.file=model.file)
print(jagsfit)
model.file <- system.file(package="R2jags", "model", "schools.txt")


data.frame(
  Parameters=c(paste("relative treatment effects ()"))
)


