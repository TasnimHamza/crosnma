#
devtools::document()  # adopt the changes in the pkg

# create all Vignette related files
library(ryxogn2)
#usethis::use_vignette("crosnma")
devtools::build_vignettes()

# add bib
knitr::write_bib(c(.packages(), "bookdown"), "packages.bib")
