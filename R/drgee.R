drgee <-
function (oformula,
          eformula,
          iaformula=formula(~1),
          olink=c("identity","log","logit"),
          elink=c("identity","log","logit"),
          estimationMethod=c("dr","obe","ebe"),
          data=NULL,
          rootFinder=findRoots,
          clusterid=NULL,
          ...
          )
{

    call <- match.call()

    olink <- match.arg(olink)

    elink <- match.arg(elink)

    estimationMethod <- match.arg(estimationMethod)

    if(estimationMethod!="obe" & olink=="logit" & elink!="logit"){
        warning("\nFor dr and ebe estimation, olink=\"logit\" can only be combined with elink=\"logit\"\nelink changed to \"logit\"\n")
        elink<-"logit"
    }else{
        elink <- match.arg(elink)
    }

    drgeeData <- drgeeData(oformula, eformula, iaformula, olink, elink,
                                data, clusterid)

    if(estimationMethod=="obe"){
        res <- obeFit(drgeeData)
        res$optim.object=NULL
    }else if(estimationMethod=="ebe"){
        res <- ebeFit(drgeeData, rootFinder, ...)
    }else if(estimationMethod=="dr"){
        res <- drFit(drgeeData, rootFinder, ...)
    }

    res$call <- call

    res$drgeeData<- drgeeData

    res$estimationMethod <- estimationMethod

    class(res) <- c("drgee")

    return(res)
}

print.drgee <-
function(x, digits = max(3L, getOption("digits") - 3L), ...){
    if(length(x$coefficients)) {
        cat("\nCoefficients for main effect:\n")
        print.default(format(coef(x), digits = digits),
                      print.gap = 2, quote = FALSE)
        cat("\n")
    } else cat("No coefficients\n\n")

}

summary.drgee <-
function(object, digits = max(3L, getOption("digits") - 3L), ...){

	s.err <- sqrt(diag(as.matrix(vcov(object))))
	zvalue <- coef(object)/s.err
	pvalue <- 2*pnorm(-abs(zvalue))

	coef.table <- as.matrix(cbind(coef(object), s.err, zvalue, pvalue))

	dimnames(coef.table) <- list(names(coef(object)), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

	model.formulas <- summary(object$drgeeData)

	ans <- list(call=object$call, coefficients=coef.table,
                    vcov=vcov(object), estimationMethod=object$estimationMethod,
                    model.formulas=model.formulas, n.obs=model.formulas$n.obs,
                    n.clust=model.formulas$n.clust)

    class(ans) <- "summary.drgee"
    return(ans)
}

print.summary.drgee <-
function(x, digits = max(3L, getOption("digits") - 3L),
             signif.stars = getOption("show.signif.stars"), ...){
        cat("\nCall:  ",
            paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

		cat("\nMain model: ", x$model.formulas$main.model,"\n")

		if(x$estimationMethod!="ebe"){
			cat("\nOutcome nuisance model: ", x$model.formulas$outcome.nuisance.model,"\n")
		}

		cat("\nOutcome link function: ", x$model.formulas$olink,"\n")

		if(x$estimationMethod!="obe"){
			cat("\nExposure nuisance model: ", x$model.formulas$exposure.nuisance.model,"\n")
			cat("\nExposure link function: ", x$model.formulas$elink,"\n")
		}

        if(length(x$coefficients)){
        	cat("\n")
            printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                         na.print = "NA", ...)
            cat("\n", x$n.obs, " complete observations used\n")

            if(x$n.clust < x$n.obs){
                cat("\nCluster-robust Std. errors\n", x$n.clust, " clusters\n")
            }
        }else{
            cat("No coefficients estimated\n\n")
        }

    }

coef.drgee <-
function(object, ...){
    return(object$coefficients)
}

vcov.drgee <-
function(object, ...){
    return(object$vcov)
}
