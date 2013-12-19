drFit <-
function(object, rootFinder=findRoots, ...){

	if(class(object)!="drgeeData"){
		stop("An object of class \"drgeeData\" is expected")
	}

    if(object$olink=="logit"){

        o.fit1 <- geeFit(object$y, cbind(object$ax, object$v), "logit")

        e.fit1 <- geeFit(object$a, cbind(object$yx, object$z), "logit")

        beta1.hat <- o.fit1$coefficients[colnames(object$ax)]
        gamma.hat <- o.fit1$coefficients[colnames(object$v)]
        beta2.hat <- e.fit1$coefficients[colnames(object$yx)]
        alpha.hat <- e.fit1$coefficients[colnames(object$z)]

        exp.gamma.V <- as.vector(exp( object$v %*% gamma.hat ))
        exp.alpha.Z <- as.vector(exp( object$z %*% alpha.hat ))

        eq.func <- function(beta, arg.list){

            exp.beta.X <- exp( arg.list$x.f%*%beta )
            exp.beta.AX <- exp( arg.list$ax.f%*%beta )

            p <- exp.beta.X * arg.list$exp.alpha.Z * ( 1 + arg.list$exp.gamma.V )
            q <- p + 1 + exp.beta.X * arg.list$exp.gamma.V
            E.star <- as.vector( p / q )

            res.A.E.star <- as.vector(arg.list$a.f-E.star)
            res.Y <- as.vector( arg.list$y.f - 1 + 1/(1 + exp.beta.AX * arg.list$exp.gamma.V) )

            return( as.vector( crossprod( arg.list$x.f, res.A.E.star * res.Y ) ) )

		}

        d.eq.func <- function(beta, arg.list){

            exp.beta.X <- as.vector(exp( arg.list$x.f%*%beta ))
            exp.beta.AX <- as.vector(exp( arg.list$ax.f%*%beta ))

            p <- exp.beta.X * arg.list$exp.alpha.Z * (1 + arg.list$exp.gamma.V)
            q <- p + 1 + exp.beta.X * arg.list$exp.gamma.V
            E.star <- p / q

            res.A.E.star <- as.vector(arg.list$a.f-E.star)
            res.Y <- as.vector( arg.list$y.f - 1 + 1/(1 + exp.beta.AX * arg.list$exp.gamma.V) )

            d.res.A.E.star <- -arg.list$x.f * (E.star / q) * ( q - p - exp.beta.X*arg.list$exp.gamma.V)
            d.res.Y <- -arg.list$ax.f * exp.beta.AX * exp.gamma.V / (1 + exp.beta.AX * arg.list$exp.gamma.V)^2

            return( crossprod( arg.list$x.f , d.res.A.E.star*res.Y + d.res.Y*res.A.E.star ) )
        }

        all.args <- c(list(beta.init=beta1.hat, eq.func=eq.func,
                           d.eq.func=d.eq.func,
                           arg.list=list(y.f=object$y, ax.f=object$ax,
                               a.f=object$a, x.f=object$x,
                               exp.gamma.V=exp.gamma.V, exp.alpha.Z=exp.alpha.Z)),
                      list(...))

        # call equation solver with beta.init as initial guess and eq.func as estimation function
   	root.object <- do.call(rootFinder, all.args)
   	beta.hat <- root.object$roots
        optim.object <- root.object$optim.object

        exp.beta.AX <- as.vector(exp( object$ax%*%beta.hat ) )
        exp.beta.X <- as.vector(exp(  object$x%*%beta.hat ) )
        exp.beta1.AX <- as.vector(exp( object$ax%*%beta1.hat ) )
        exp.beta2.YX <- as.vector(exp( object$yx%*%beta2.hat ) )
        exp.gamma.V <- as.vector(exp( object$v%*%gamma.hat ) )
        exp.alpha.Z <- as.vector(exp( object$z%*%alpha.hat ) )

        p <- as.vector(exp.beta.X*exp.alpha.Z*(1+exp.gamma.V))
        q <- as.vector(p+1+exp.beta.X*exp.gamma.V)
        E.star <- as.vector(p/q)

        res.A.E.star <- as.vector(object$a-E.star)
        res.Y <- as.vector(object$y - 1 + 1/(1 + exp.beta.AX*exp.gamma.V))

        d.res.A.E.star.gamma <- -object$v * exp.gamma.V * (exp.beta.X/q) * (exp.alpha.Z+E.star*(exp.alpha.Z-1))
        d.res.Y.gamma <- -object$v * exp.beta.AX * exp.gamma.V / (1 + exp.beta.AX*exp.gamma.V)^2

        d.res.A.E.star.alpha <- -object$z * E.star * (1-E.star)

        U <- cbind(object$x * res.A.E.star * res.Y,
                   cbind(object$ax, object$v) * o.fit1$res,
                   cbind(object$yx, object$z) * e.fit1$res)

        d.U1.beta <- d.eq.func(beta.hat, all.args$arg.list)
        d.U1.beta1.beta2 <-
            matrix(rep(0, ncol(object$x)*(ncol(object$ax)+ncol(object$yx))), nrow=ncol(object$x))
        d.U1.gamma <- crossprod( object$x , d.res.A.E.star.gamma*res.Y + d.res.Y.gamma*res.A.E.star )
        d.U1.alpha <- crossprod( object$x , d.res.A.E.star.alpha*res.Y )

        d.U2.beta <- matrix(rep(0,(ncol(object$ax)+ncol(object$v))*ncol(object$ax)),ncol=ncol(object$ax))
        d.U2.beta1 <- crossprod( cbind(object$ax, object$v), o.fit1$d.res[,colnames(object$ax)])
        d.U2.beta2 <- matrix(rep(0,(ncol(object$yx)+ncol(object$v))*ncol(object$ax)),ncol=ncol(object$yx))
        d.U2.gamma <- crossprod( cbind(object$ax, object$v), o.fit1$d.res[,colnames(object$v)])
        d.U2.alpha <- matrix(rep(0,(ncol(object$ax)+ncol(object$v))*ncol(object$z)),ncol=ncol(object$z))

        d.U3.beta <- matrix(rep(0,(ncol(object$yx)+ncol(object$z))*ncol(object$ax)),ncol=ncol(object$ax))
        d.U3.beta1 <- matrix(rep(0,(ncol(object$yx)+ncol(object$z))*ncol(object$ax)),ncol=ncol(object$ax))
        d.U3.beta2 <- crossprod( cbind(object$yx, object$z), e.fit1$d.res[,colnames(object$yx)])
        d.U3.gamma <- matrix(rep(0,(ncol(object$yx)+ncol(object$z))*ncol(object$v)),ncol=ncol(object$v))
        d.U3.alpha <- crossprod( cbind(object$yx, object$z), e.fit1$d.res[,colnames(object$z)])

	d.U <- rbind( cbind(d.U1.beta, d.U1.beta1.beta2, d.U1.gamma, d.U1.alpha),
                     cbind(d.U2.beta, d.U2.beta1, d.U2.beta2, d.U2.gamma, d.U2.alpha),
                     cbind(d.U3.beta, d.U3.beta1, d.U3.beta2, d.U3.gamma, d.U3.alpha) ) / nrow(U)

    }else{

    o.fit1 <- geeFit(object$y, cbind(object$ax, object$v), object$olink)

    e.fit <- geeFit(object$a, object$z, object$elink)

    beta1.hat <- o.fit1$coefficients[colnames(object$ax)]
    gamma.hat <- o.fit1$coefficients[colnames(object$v)]
    alpha.hat <- e.fit$coefficients

    if(object$olink=="identity"){
    	w <- t(object$x * e.fit$res)
        beta.hat <- as.vector(solve(w%*%object$ax) %*% w %*% (object$y-object$v%*%gamma.hat))
        optim.object <- NULL
        res.Y <- as.vector(object$y - object$ax%*%beta.hat - object$v%*%gamma.hat)
        d.res.Y.beta <- -object$ax
        d.res.Y.gamma <- -object$v

    }else if(object$olink=="log"){

        eq.func <- function(beta, arg.list){
            crossprod(arg.list$w, arg.list$y*exp(-arg.list$ax%*%beta) - arg.list$exp.v.gamma)
        }

        d.eq.func <- function(beta, arg.list){
            crossprod(arg.list$w*as.vector(arg.list$y*exp(-arg.list$ax%*%beta)), -arg.list$ax)
        }

        w <- object$x * e.fit$res

        all.args <- c(list(beta.init=beta1.hat, eq.func=eq.func,
                           d.eq.func=d.eq.func,
                           arg.list=list(w = w,
                               exp.v.gamma = exp(object$v%*%gamma.hat),
                               ax=object$ax, y=object$y)),
                      list(...))

   	# call equation solver with beta.init as initial guess and eq.func as estimation function
        root.object <- do.call(rootFinder, all.args)
        beta.hat <- root.object$roots
        optim.object <- root.object$optim.object
    	res.Y <- as.vector(object$y * exp(-object$ax %*% beta.hat) - exp(object$v%*%gamma.hat))
    	d.res.Y.beta <- -object$ax * as.vector((object$y) * exp(-object$ax %*% beta.hat))
    	d.res.Y.gamma <- -object$v * as.vector(exp(object$v%*%gamma.hat))
    }

    U <- cbind( object$x * e.fit$res * res.Y, cbind(object$ax, object$v) * o.fit1$res, object$z * e.fit$res )

    d.U1.beta <- crossprod(object$x * e.fit$res, d.res.Y.beta)
    d.U1.beta1 <- matrix(rep(0,ncol(object$x)*ncol(object$ax)), nrow=ncol(object$x))
    d.U1.gamma <- crossprod(object$x * e.fit$res, d.res.Y.gamma)
    d.U1.alpha <- crossprod(object$x * res.Y, e.fit$d.res)

    d.U2.beta <- matrix(rep(0,(ncol(object$ax)+ncol(object$v))*ncol(object$ax)), ncol=ncol(object$ax))
    d.U2.beta1.gamma <- crossprod(cbind(object$ax, object$v), o.fit1$d.res)
    d.U2.alpha <- matrix(rep(0,(ncol(object$ax)+ncol(object$v))*ncol(object$z)), ncol=ncol(object$z))

    d.U3.beta.beta1.gamma <- matrix(rep(0,ncol(object$z)*(2*ncol(object$ax)+ncol(object$v))),nrow=ncol(object$z))
    d.U3.alpha <- crossprod(object$z, e.fit$d.res)

    d.U <- rbind( cbind(d.U1.beta, d.U1.beta1, d.U1.gamma, d.U1.alpha),
                     cbind(d.U2.beta, d.U2.beta1.gamma, d.U2.alpha),
                     cbind(d.U3.beta.beta1.gamma, d.U3.alpha) ) / nrow(U)

}

    vcov <- robVcov(U, d.U, object$id)[1:ncol(object$ax), 1:ncol(object$ax),drop=FALSE]

    names(beta.hat) <- colnames(object$ax)

    dimnames(vcov) <- list(names(beta.hat), names(beta.hat))

    result <- list(coefficients=beta.hat, vcov=vcov,
                optim.object=optim.object)

    return(result)

}
