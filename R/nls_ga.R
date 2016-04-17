## nls.lm.control is called from nlsLM, used to find convergence parameters
## for the nlsLM algorithm.
nls.ga.control <- function(ftol  = sqrt(.Machine$double.eps),
                           frtol = sqrt(.Machine$double.eps),
                           ptol  = sqrt(.Machine$double.eps), 
                           prtol = sqrt(.Machine$double.eps),
                           ctol  = sqrt(.Machine$double.eps),
                           crtol = sqrt(.Machine$double.eps),
                           gtol = sqrt(.Machine$double.eps), 
                           diag = numeric(0), epsfcn1 = sqrt(.Machine$double.eps), 
                           epsfcn2 = sqrt(.Machine$double.eps), factor = 100, 
                           maxfev = integer(), maxjev = integer(),
                           maxaev = integer(),
                           maxiter = 50, nprint = 0){
    
    list(ftol = ftol, ptol = ptol, gtol = gtol, diag = diag, 
         epsfcn1 = epsfcn1, epsfcn2 = epsfcn2, factor = factor, 
         maxfev = maxfev, maxiter = maxiter, nprint = nprint)
}

## nls.lm is called from nlsLM and calls the optimisation routines through
## the external call to nls_lm.
nls.ga <- function(par, lower=NULL, upper=NULL, fn, jac = NULL, acc = NULL, 
                   control = nls.lm.control(), ...){
    
    # Make a copy of the function provided as input
    fn1  <- function(par){
        fn(par, ...)
    }
    
    # Make a copy of the jacobian function (if supplied)
    jac1 <- if(!is.null(jac)){
        function(par) jac(par, ...)
    }
    
    acc1 <- if(!is.null(acc)){
        function(par) acc(par, ...)
    }
    
    # Set lower parameter bounds (if none specified)
    if(is.null(lower)){
        lower <- rep(-Inf,length(par))
    }
    
    # Set upper parameter bounds (if none specified)
    if(is.null(upper)){
        upper <- rep(Inf,length(par))
    }
    
    # Error message for not specfying a lower bound for all parameters
    if(length(lower) != length(par)){
        stop("length(lower) must be equal to length(par)")
    }
    
    # Error message for not specfying an upper bound for all parameters
    if(length(upper) != length(par)){
        stop("length(upper) must be equal to length(par)")
    }
    
    # Get nls_lm convergence criteria
    ctrl <- nls.ga.control()
    
    # If additional control options were specified use them along side the
    # defaults.
    if(!missing(control)){
        control              <- as.list(control)
        ctrl[names(control)] <- control
    }
    
    # Specify max function evaluations
    if(length(ctrl[["maxfev"]]) == 0){
        ctrl[["maxfev"]] <- 100*(length(unlist(par)) + 1)
    }
    
    # Specify max Jacobian evaluations
    if(length(ctrl[["maxjev"]]) == 0){
        ctrl[["maxjev"]] <- 100*(length(unlist(par)) + 1)
    }
    
    # Specify max Acceleration evaluations
    if(length(ctrl[["maxaev"]]) == 0){
        ctrl[["maxaev"]] <- 100*(length(unlist(par)) + 1)
    }
    
    # Use nls_lm to generate optimised parameter values.
    sink(file = "geo_lm_debug.txt")
    out <- .Call("geo_lm", par, lower, upper, fn1, jac1, acc1, 
                 ctrl, new.env())
                 # PACKAGE = "geodesicLM")
    sink()
    # Convert hessian output into a matrix
    out$hessian <- matrix(out$hessian, nrow = length(unlist(par)))
    
    # This section presents an error when parameters are given to the
    # model as a vector.
    names(out$par)        <-
        rownames(out$hessian) <-
        colnames(out$hessian) <-
        names(out$diag)       <- names(par)
    out$deviance <- sum(out$fvec^2)
    class(out)   <- "nls.lm"
    out
}

print.nls.lm <- function(x, digits = max(3, getOption("digits") - 3), ...){
    
    cat("Nonlinear regression via the Levenberg-Marquardt algorithm\n")
    cat("parameter estimates:", toString(x$par), "\n")
    cat("residual sum-of-squares: ", format(x$deviance, digits = digits),
        "\n", sep = '')
    cat("reason terminated: ", x$message, "\n", sep='')
    invisible(x)
}

deviance.nls.lm    <- function(object, ...) object$deviance
coef.nls.lm        <- function(object, ...) unlist(object$par)
residuals.nls.lm   <- function(object, ...) object$fvec
df.residual.nls.lm <- function(object, ...){
    length(resid(object)) - length(coef(object))
}

summary.nls.lm <- function (object, ...){
    
    param     <- coef(object)
    pnames    <- names(param)
    ibb       <- chol(object$hessian)
    ih        <- chol2inv(ibb)
    p         <- length(param)
    rdf       <- length(object$fvec) - p 
    resvar    <- deviance(object) / rdf
    se        <- sqrt(diag(ih) * resvar)
    names(se) <- pnames
    tval      <- param/se
    
    param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <- list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans <- list(residuals = object$fvec, sigma = sqrt(object$deviance/rdf),
                df = c(p, rdf), cov.unscaled = ih,
                info = object$info, niter = object$niter,
                stopmess = object$message, 
                coefficients = param)
    class(ans) <- "summary.nls.lm"
    ans
}

print.summary.nls.lm <- function (x, digits = max(3, getOption("digits") - 3), ...){
    
    df <- x$df
    rdf <- df[2]
    cat("\nParameters:\n")
    printCoefmat(x$coefficients, digits = digits, ...)
    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
    cat("Number of iterations to termination:", x$niter, "\n")
    cat("Reason for termination:", x$stopmess, "\n")
    invisible(x)
    
}

vcov.nls.lm <- function(object,...){
    object$deviance/(length(object$fvec)-length(object$par))*solve(object$hessian)
}

confint.nls.lm <- function(object, parm, level = 0.95, ...){
    cc <- coef(object)
    
    if(missing(parm)){
        parm <- seq_along(cc)
    }
    
    levs         <- c((1-level)/2,0.5+level/2)
    dfval        <- (length(object$fvec)-length(object$par))
    tdist        <- qt(levs[2],dfval)
    m1           <- outer(sqrt(diag(vcov(object))),c(-1,1))*tdist
    m2           <- sweep(m1,cc,MARGIN=1,"+")
    colnames(m2) <- paste(100*levs,"%",sep="")
    m2[parm,]
    
}