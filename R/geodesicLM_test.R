
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

# Generate problem
n  <- 2e1
x  <- runif(n, 0, 10)
y  <- sin(x) + rnorm(n, 0, 0.2)
data <- data.frame(x = x, y = y)

myform <- y ~ a * sin(b * x + c) + d
init   <- c(a = 2, b = 0.5, c = 1, d = -0.5)

# control parameters for geodesicLM
control_args <- nls.ga.control(maxiter = 200, print_level = 2, iaccel = 0)

# Fit both models
fit.minpack <- nlsLM(formula = myform, data = data, start = init,
                     control = nls.lm.control(maxiter = 200))
fit.geodesic <- nlsLM2(formula = myform, data = data, start = init,
                       control = control_args)

dyn.unload(dllPath)
