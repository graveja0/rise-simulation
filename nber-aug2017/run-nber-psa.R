rm(list=ls())
vN_PSA = 2
pkg = list("tidyverse","deSolve")
invisible(lapply(pkg, require, character.only = TRUE))
run.id.base <- "vogi-numerical"
Scenarios <- list (c("None","Single"),c("None","None"),c("None","Panel"),c("Panel","None"))

mainDir <- "../run-data"
subDir <- run.id.base

if (!file.exists(file.path(mainDir,subDir))){
    dir.create(file.path(mainDir, subDir))
  }
