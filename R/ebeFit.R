ebeFit <-
function(object, rootFinder=findRoots, ...){

	if(class(object)!="drgeeData"){
		stop("An object of class \"drgeeData\" is expected")
	}

    if(ncol(object$v)>1)  warning("\nIn ebe estimation,  RHS of oformula is ignored\n")

    if(object$olink=="logit"){
    	if(object$elink!="logit"){
    		stop("In exposure based estimation with olink \"logit\", elink needs to be \"logit\"\n" )
    	}
    	if(length(unique(object$a))!=2){
    		stop("When elink is \"logit\" the exposure needs to be binary")
    	}

        fit <- geeFit(object$a, cbind(object$yx, object$z), "logit")
    	beta.hat <- fit$coefficients[colnames(object$yx)]
        U <- cbind(object$yx, object$z) * fit$res
        d.U <- crossprod( cbind(object$yx, object$z), fit$d.res)  / nrow(U)
        optim.object <- NULL

	# If the outcome link is identity or log
    }else{
	# Calculate the outcome nuisance parameters
    	e.fit <- geeFit(object$a,object$z,object$elink)

        if(object$olink=="identity"){
            w <- t(object$x * e.fit$res)
            beta.hat <- as.vector(solve(w%*%object$ax) %*% w %*% object$y)
            optim.object <- NULL
            S <- as.vector(object$y - object$ax %*% beta.hat)
            d.S <- -object$ax
            optim.object <- NULL
        }else if(object$olink=="log"){
            eq.func <- function(beta, arg.list){
                crossprod(arg.list$w, exp(-arg.list$ax%*%beta) )
            }

            d.eq.func <- function(beta, arg.list){
                crossprod(arg.list$w*as.vector(exp(-arg.list$ax%*%beta)), -arg.list$ax)
            }

	    # Create an initial guess for the values of beta
            o.fit1 <- geeFit(object$y, cbind(1,object$ax), "log")
            beta.init <- o.fit1$coefficients[-1]
            w <- object$x * e.fit$res * as.vector(object$y)
            all.args <- c(list(beta.init=beta.init, eq.func=eq.func, d.eq.func=d.eq.func,
                          arg.list=list(w=w, ax=object$ax)),
                          list(...))

    	    # call equation solver with beta.init as initial guess and eq.func as estimation function
            root.object <- do.call(rootFinder, all.args )
            beta.hat <- root.object$roots
            optim.object <- root.object$optim.object
            S <- as.vector(object$y * exp(-object$ax %*% beta.hat))
            d.S <- -object$ax * as.vector(object$y * exp(-object$ax%*%beta.hat))
    	}

        U <- cbind( object$x * e.fit$res * S, object$z * e.fit$res)

        U1.d.beta <- crossprod(object$x * e.fit$res, d.S)
        U1.d.alpha <- crossprod(object$x * S, e.fit$d.res)
        U2.d.beta <- matrix(rep(0,ncol(object$ax)*ncol(object$z)), nrow=ncol(object$z) )
        U2.d.alpha <- crossprod(object$z, e.fit$d.res)
        d.U <- rbind( cbind( U1.d.beta, U1.d.alpha), cbind(U2.d.beta, U2.d.alpha) ) / nrow(U)
    }

    vcov <- robVcov(U, d.U, object$id)[1:ncol(object$ax), 1:ncol(object$ax),drop=FALSE]

    names(beta.hat) <- colnames(object$ax)

    dimnames(vcov) <- list(names(beta.hat), names(beta.hat))

    result <- list(coefficients=beta.hat, vcov=vcov,
                optim.object=optim.object)

    return(result)
}
