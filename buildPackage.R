rm(list = ls())  # vymazu Envrironment
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


suppressPackageStartupMessages(library(devtools))

################################################################################
# 

devtools::document()
devtools::install()