
mypkg <- c("shiny", "magrittr", "dplyr", "tidyr", "ggplot2", "cowplot", "gtools")

if (length(setdiff(mypkg, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(mypkg, rownames(installed.packages())), repos = "https://cran.r-project.org/")
}
