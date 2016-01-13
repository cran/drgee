gee <-
    function (formula,
              link = c("identity","log","logit"),
              data,
              cond = FALSE,
              clusterid,
              rootFinder = findRoots,
              ...
              ){

        call <- match.call()

        link <- match.arg(link)

        m <- match(c("data", "clusterid", "cond"), names(call), 0L)

        dD <- call[c(1L, m)]
        dD$oformula <- formula
        dD$olink <- link
        dD$estimation.method <- "o"

        dD[[1L]] <- quote(drgeeData)

        gee.data <- eval(dD, parent.frame())

        if (cond) {
            if (is.null(gee.data$v)) {
                stop("No parameters to estimate\n\n")
            } else {
                fit <- geeFitCond(gee.data$y, gee.data$v, link = link, gee.data$id,
                                  rootFinder, ...)
                naive.var <- fit$naive.var
            }

        } else {
            fit <- geeFit(gee.data$y, gee.data$v, link = link)
            naive.var <- NULL
        }

        coefficients = fit$coefficients
        names(coefficients) <- gee.data$v.names

        u <- fit$eq.x * fit$res

        d.u <- crossprod( fit$eq.x , fit$d.res ) / nrow(u)

        vcov <- as.matrix( robVcov(u, d.u, gee.data$id) )

        dimnames(vcov) <- list(gee.data$v.names, gee.data$v.names)

        result <- list(coefficients = coefficients,
                       vcov = vcov,
                       call = call,
                       cond = cond,
                       gee.data = gee.data,
                       optim.object = fit$optim.object,
                       res = fit$res,
                       d.res = fit$d.res,
                       naive.var = naive.var)

        class(result) <- "gee"

        return(result)
    }

print.gee <-
    function(x, digits = max(3L, getOption("digits") - 3L), ...) {
        if (length(x$coefficients)) {
            cat("\nCoefficients:\n")
            print.default(format(coef(x), digits = digits),
                          print.gap = 2, quote = FALSE)
            cat("\n")
        } else {
            cat("No coefficients\n\n")
        }

    }

summary.gee <-
    function(object, digits = max(3L, getOption("digits") - 3L), ...) {

	s.err <- sqrt(diag(as.matrix(vcov(object))))
	zvalue <- coef(object) / s.err
	pvalue <- 2 * pnorm(-abs(zvalue))

	coef.table <- as.matrix(cbind(coef(object), s.err, zvalue, pvalue))

	dimnames(coef.table) <- list(names(coef(object)), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

        summ <- summary(object$gee.data)

	ans <- list(call = object$call,
                    coefficients = coef.table,
                    vcov = vcov(object),
                    model.formula = summ$outcome.nuisance.model,
                    link = summ$olink,
                    n.obs = summ$n.obs,
                    n.clust = summ$n.clust)

        class(ans) <- "summary.gee"
        return(ans)
    }

print.summary.gee <-
    function(x, digits = max(3L, getOption("digits") - 3L),
             signif.stars = getOption("show.signif.stars"), ...){
        cat("\nCall:  ",
            paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

        cat("\nModel: ",x$model.formula,"\n")

        cat("\nLink function: ", x$link,"\n")

        if (length(x$coefficients)) {
            cat("\n")
            printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                         na.print = "NA", ...)
            cat("\n", x$n.obs, " complete observations used\n")

            if (x$n.clust < x$n.obs){
                cat("\nCluster-robust Std. errors\n", x$n.clust, " clusters\n")
            }
        } else {
            cat("No coefficients estimated\n\n")
        }

    }

coef.gee <- function(object, ...) {
    return(object$coefficients)
}

vcov.gee <- function(object, ...) {
    return(object$vcov)
}

