obeFit <-
    function(object, ...){

        if(class(object)!="drgeeData"){
            stop("An object of class \"drgeeData\" is expected")
        }

        if(ncol(object$z)>1)  warning("\nIn obe estimation,  RHS of eformula is ignored\n")

        fit <- geeFit(y=object$y, x=cbind(object$ax, object$v), link=object$olink, ...)

        beta.hat <- fit$coefficients[colnames(object$ax)]

        U <- cbind(object$ax, object$v) * fit$res

        d.U <- crossprod( cbind(object$ax, object$v) , fit$d.res ) / nrow(U)

        vcov <- robVcov(U, d.U, object$id)[1:ncol(object$ax), 1:ncol(object$ax),drop=FALSE]

        dimnames(vcov) <- list(names(beta.hat), names(beta.hat))

        result <- list(coefficients=beta.hat, vcov=vcov)

        return(result)
}

