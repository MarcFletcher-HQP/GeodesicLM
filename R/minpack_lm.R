nlsLM2 <- nlsLM
body(nlsLM2)[[28]] <- quote(NLS <- nls.ga(par = start, fn = FCT, jac = jac, control = control,
                                           lower = lower, upper = upper, ...))