# Trims only the reproducible output from a GLRT execution, namely
# the first iteration only. The first iteration of the algorithm has
# no sampling, so the values should match.

library(assertthat)
library(tidyverse)

source(str_c(Sys.getenv('INVAR_HOME'), '/R/shared/common.R'))

args <- commandArgs(TRUE)

assert_that(length(args) == 2, msg = "Expect exactly two arguments: the generated file and the output file")
assert_that(file.exists(args[1]), msg = str_c(args[1], " does not exist."))

readRDS(args[1]) %>%
    filter(ITERATION == 1) %>%
    saveRDS(args[2])
