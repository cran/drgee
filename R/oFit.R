oFit <-
    function(object, ...) {

        if(class(object) != "drgeeData") {
            stop("An object of class \"drgeeData\" is expected")
        }

        axv <- cbind(object$ax, object$v)

        if (object$cond) {
            fit <- geeFitCond(y = object$y, x = axv, link = object$olink, object$id,
                              ...)
        } else {
            fit <- geeFit(y = object$y, x = axv, link = object$olink)
        }

        u <- apply(fit$eq.x, 2, '*', fit$res)
        d.u <- crossprod( fit$eq.x , fit$d.res ) / nrow(u)

        coefficients <- fit$coefficients
        names(coefficients) <- colnames(axv)

        vcov <- robVcov(u, d.u, object$id)
        dimnames(vcov) <- list(names(coefficients), names(coefficients))


        result <- list(coefficients = coefficients, vcov = vcov, optim.object = fit$optim.object,
                       optim.object.o = fit$optim.object, optim.object.e = NULL)

        return(result)
    }
