
devtools::document()  # adopt the changes in the pkg
devtools::build()
loadNamespace("crosnma")
# create all Vignette related files
library(ryxogn2)
#usethis::use_vignette("crosnma")
devtools::build_vignettes()





