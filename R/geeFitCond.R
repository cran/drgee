geeFitCond <-
    function(y, x, link = c("identity","log","logit"), id, rootFinder =
             findRoots, ...) {

        link <- match.arg(link)

        x.cent <- x - apply(x, 2, function(z) ave(z,id))

        if (link == "identity") {

            y.cent <- y - ave(y,id)

            ## Solve the equation     t(instr) %*% (  y.cent - x.cent %*% beta ) = 0
            ## t(instr) %*% y.cent = t(instr) %*% x.cent %*% beta
            ## beta.hat = solve( crossprod(instr, x.cent) ) %*% t(instr) %*% y.cent

            beta.hat <- as.vector( solve(crossprod(x.cent, x.cent)) %*% crossprod(x.cent, y.cent) )

            return( list(coefficients = beta.hat,
                         res = as.vector(y.cent - x.cent %*% beta.hat),
                         d.res = -x.cent,
                         eq.x = x.cent,
                         optim.object = NULL,
                         naive.var = NULL))

        } else if (link == "log") {

            ## Create an initial estimate of beta
            ## with intercept being the mean of the outcome
            intercept.init <- colMeans(y)
            beta.init <- geeFit(y = y, x = cbind(rep(intercept.init, nrow(x)), x), link = link)$coefficients[-1]

            u.func <- function(beta, arg.list) {
                y.star <- arg.list$y * exp(-arg.list$x %*% beta)
                y.star.cent <- y.star - ave(y.star, arg.list$id)
                as.vector( crossprod( arg.list$x.cent , y.star.cent))
            }

            all.args <- c(list(beta.init = beta.init, eq.func = u.func, d.eq.func = NULL,
                               arg.list = list(y = y, x = x, x.cent = x.cent, id = id)),
                          list(...))

            ## Call equation solver with beta.init as initial guess and eq.func as estimation function
            root.object <- do.call(rootFinder, all.args)
            beta.hat <- root.object$roots

            y.star.hat <- as.vector(y * exp(-x %*% beta.hat))
            y.star.hat.cmean <- ave(y.star.hat, id)

            d.y.star.hat.beta <- y.star.hat * -x
            d.y.star.hat.beta.cmean <- apply(d.y.star.hat.beta, 2, function(z) ave(z,id))

            res = y.star.hat - y.star.hat.cmean
            d.res = d.y.star.hat.beta -  d.y.star.hat.beta.cmean

            names(beta.hat) <- colnames(x)

            return( list(coefficients = beta.hat,
                         res = res,
                         d.res = d.res,
                         eq.x = x.cent,
                         optim.object = root.object$optim.object,
                         naive.var = NULL) )

        } else if (link == "logit") {

            fit <- try(survival::clogit(y~x.cent + strata(id), method = "exact") )
            root.object <- list()
            root.object$optim.object <- NULL

            if (class(fit)[1] == "try-error") {

                beta.hat <- rep(NA, ncol(x))
                names(beta.hat) <- colnames(x)

                eq.res <- rep(NA, nrow(x))

                d.eq.res <- matrix( rep(NA, length(x)), nrow = nrow(x) )
                colnames(d.eq.res) <- colnames(x)

                res.list <- list(eq.res = eq.res, d.eq.res = d.eq.res)

                naive.var = NULL
            } else {
                beta.hat <- coef(fit)

                names(beta.hat) <- colnames(x)

                res.list <- clogitRes(beta.hat, y, x.cent, id )

            }

            return( list(coefficients = beta.hat,
                         res = res.list$eq.res,
                         d.res = res.list$d.eq.res,
                         eq.x = x.cent,
                         optim.object = root.object$optim.object,
                         naive.var = fit$var) )


        }

    }
