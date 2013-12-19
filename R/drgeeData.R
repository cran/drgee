drgeeData <-
    function(oformula,
             eformula,
             iaformula=formula(~1),
             olink=c("identity","log","logit"),
             elink=c("identity","log","logit"),
             data=NULL,
             clusterid=NULL){

    olink <- match.arg(olink)

    elink <- match.arg(elink)

    oterms <- terms(oformula)
    Y <- as.character(oformula[[2]])
    V <- attr(oterms, "term.labels")

    eterms <- terms(eformula)
    A <- as.character(eformula[[2]])
    Z <- attr(eterms, "term.labels")

    X <- attr(terms(iaformula),"term.labels")

    if(is.null(data)){
        data <- attr(oterms, ".Environment")
    }

    if(is.null(clusterid)){
        tmp.formula <- formula(paste("~", paste(c(Y, A, X, V, Z),
                                                  collapse="+")), sep="")
        # Create a dataset taken from dataset or environment of a formula with nas removed
        mf <- model.frame(formula=tmp.formula, data=data, subset=NULL,
                          na.action=na.omit)
        id <- as.factor(1:nrow(mf))

    }else{
        # Create a dataset taken from dataset or environment of a formula with nas removed
        tmp.formula <- formula(paste("~", paste(c(Y, A, X, V, Z, clusterid),
                                                  collapse="+")), sep="")
        mf <- model.frame(formula=tmp.formula, data=data, subset=NULL,  na.action=na.omit)
        id <- as.factor(unlist(mf[clusterid]))
    }

    # Check factor variables
    if(is.factor(mf[[Y]])){
        if(olink!="logit"){
            stop("Factor outcome is only possible when outcome link is logit\n\n")
        }else if(nlevels(mf[[Y]])>2){
            stop("The outcome cannot be a factor with more than two levels\n\n")
        }
    }else if(olink=="log" & min(mf[[Y]])<0){
        stop("When outcome link is log the outcome needs to be positive\n\n")
    }


    # Extract outcome and exposure
    y <- model.matrix(formula(paste("~",Y)), mf)[,2,drop=FALSE]
    a <- model.matrix(formula(paste("~",A)), mf)[,-1,drop=FALSE]

    # Create design matrices
    v <- model.matrix(oformula,mf)
    z <- model.matrix(eformula,mf)
    x <- model.matrix(iaformula,mf)

    n.x <- ncol(x)

    # Create new factors for y*x
    n.a <- ncol(a)
    ax <- matrix(rep(NA,nrow(a)*n.a*n.x),nrow=nrow(a))
    # The first n.a columns of ax will be a
    ax[,1:n.a] <- a
    colnames(ax) <- rep(colnames(a),n.x)
    # The subsequent columns will be multiplications of a with x
    if(n.x>1){
        for(i in 2:n.x){
            # Create new factors for a*x
            ax[,(i-1)*n.a+(1:n.a)] <- a*x[,i]
            colnames(ax)[(i-1)*n.a+(1:n.a)] <- paste(colnames(a),".",colnames(x)[i],sep="")
        }
    }

    # Create new factors for y*x
    yx <- as.matrix(apply(x,2,"*",y))
    colnames(yx)[1] <- colnames(y)
    colnames(yx)[-1] <- paste(colnames(y),".",colnames(x)[-1],sep="")

    drgeeData <- list(y=y, a=a, x=x, ax=ax, v=v, z=z, yx=yx, id=id,
                       olink=olink, elink=elink)

    class(drgeeData) <- "drgeeData"
    return(drgeeData)
}

summary.drgeeData <- function(object, ...){

    outcome <- colnames(object$y)

    exposure <- colnames(object$a)

    interactions <- colnames(object$x)

    main.model <- paste(colnames(object$y), "~",
                          paste(colnames(object$ax), collapse=" + "), sep=" ")

    outcome.nuisance.model <- paste(colnames(object$y),"~",
                                      paste(c(colnames(object$v)[-1]),
                                            collapse=" + "), sep=" ")

    outcome.model <- paste(colnames(object$y),"~",
                             paste(c(colnames(object$ax), colnames(object$v)[-1]),
                                   collapse=" + "), sep=" ")

    exposure.nuisance.model <- paste(colnames(object$a),"~",
                                           paste(colnames(object$z)[-1],
                                                 collapse=" + "), sep=" ")

    exposure.model <- paste(colnames(object$a),"~",
                                           paste(c(colnames(object$yx), colnames(object$z)[-1]),
                                                 collapse=" + "), sep=" ")

    ans <- list(outcome=outcome,
                exposure=exposure,
                interactions=interactions,
                main.model=main.model,
                outcome.nuisance.model=outcome.nuisance.model,
                outcome.model=outcome.model,
                exposure.nuisance.model=exposure.nuisance.model,
                exposure.model=exposure.model,
                n.obs=length(object$id),
                n.clust=nlevels(object$id),
                olink=object$olink,
                elink=object$elink)

    class(ans) <- "summary.drgeeData"

    return(ans)
}

print.summary.drgeeData <-
    function(x, digits = max(3L, getOption("digits") -
                    3L), ...){

    cat("\nOutcome: ", colnames(x$outcome), "\nExposure: ", colnames(x$exposure),
        "\nInteractions: ", paste(colnames(x$interactions), collapse=", "))

    cat("\n\nMain model: ", x$main.model, "\nwith link function: ", x$olink)

    cat("\n\nOutcome nuisance model: ", x$outcome.nuisance.model, "\nwith link function: ", x$olink)

    cat("\n\nExposure nuisance model: ", x$exposure.nuisance.model, "\nwith link function: ", x$elink, "\n")

    cat("\n\n", x$n.obs, "with ", x$n.clust, " clusters\n")
}

