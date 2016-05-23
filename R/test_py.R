## This script is an R adaptation of a test script to the python interface for GeodesicLM
## See: https://sourceforge.net/p/geodesiclm/code/ci/master/tree/pythonInterface/test.py

# Clear Workspace 
rm(list = ls())

setwd("P:/RIAP/RCS/Marc_Fletcher/C-code/GeodesicLM-master/GeodesicLM-master/R")

# Attach Libraries
library(magrittr)
library(minpack.lm)

# Set file paths and attach DLL
baseDir <- getwd() %>% dirname
dllPath <- paste(baseDir, "src", "geodesicLM.dll", sep = "/")
dyn.load(dllPath)

# Source helper files 
source("nls_ga.R")
source("minpack_lm.R")


rosenbrock <- function(x0, x1, A){
    c(1 - x0, A*(x1-x0^2))
}

init <- list(x0 = -1, x1 = 1)

control_args <- nls.ga.control(maxiter = 1e4, print_level = 5, iaccel = 1, avmax = 2)
nlsLM(formula = y ~ rosenbrock(x0, x1, 1000), data = list(y = c(0,0)), start = init, 
       control = control_args)
