
devtools::document()  # adopt the changes in the pkg

# create all Vignette related files
library(ryxogn2)
#usethis::use_vignette("crosnma")
devtools::build_vignettes()




success <- 0:20

plot(success,dbinom(success,size=20,prob=.3),
     type='h',
     main='Binomial Distribution (n=20, p=0.3)',
     ylab='Probability',
     xlab ='# events',
     lwd=3,
     las=1)
success <- 0:1
plot(0:1,dbinom(success,size=1,prob=.3),
     type='h',
     main='Binomial Distribution (n=1, p=0.3)',
     ylab='Probability',
     xlab ='# events',
     lwd=3,
     las=1)
