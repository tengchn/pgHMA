library(devtools)
library(roxygen2)

load_all()
document()
check()
build()
build_manual()
build_site()
build_vignettes()
install("../pgHMAtools/")
library(pgHMAtools)

