eFit <-
    function(object, rootFinder = findRoots, ...){

	if (class(object) != "drgeeData") {
            stop("An object of class \"drgeeData\" is expected")
	}

        if (object$olink == "logit") {

            if (object$elink != "logit") {
                warning("\nAssuming the logit link for the exposure nuisance model\n\n")
            }

            if (length(unique(object$a)) != 2) {
    		stop("When elink is \"logit\" the exposure needs to be binary")
            }

            ## Let y and a switch place and replace v with z
            ## and run retrospective logistic regression
            outcome <- object$outname
            object$outname <- object$expname
            object$expname <- outcome
            object$y <- object$a
            colnames(object$y) <- colnames(object$a)

            exposures.all <- colnames(object$ax)
            object$ax <- object$yx
            colnames(object$ax) <- exposures.all

            object$v <- object$z
            colnames(object$v) <- colnames(object$z)
            object$oterms <- object$eterms

            return(oFit(object))

            ## If the outcome link is identity or log
        } else {
            if (object$cond) {
                return(dreFitCond(object, omodel = FALSE, rootFinder = findRoots, ...))
            } else {
                return(dreFit(object, omodel = FALSE, rootFinder = findRoots, ...))
            }
        }
    }
