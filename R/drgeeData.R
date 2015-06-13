drgeeData <-
    function(outcome,
             exposure,
             oformula,
             eformula,
             iaformula = formula(~1),
             olink = c("identity","log","logit"),
             elink = c("identity","log","logit"),
             data,
             estimation.method = c("dr", "o", "e"),
             cond = FALSE,
             clusterid
             ) {

        call <- match.call()

        olink <- match.arg(olink)
        elink <- match.arg(elink)
        estimation.method <- match.arg(estimation.method)

        response.ok <- function (x, link = "identity") {

            if (is.factor(x)) {

                if (link != "logit"){
                    stop("\nFactor response only possible for the logit link\n\n");
                    return(FALSE);
                }

                if(nlevels(x) != 2){
                    stop("\nThe response in logistic regression needs to be binary\n\n");
                    return(FALSE);
                } else {
                    return(TRUE);
                }

            } else if (is.vector(x) & is.numeric(x)) {

                if (link == "log"){

                    if (min(x, na.rm = TRUE) < 0) {
                        stop("\nThe response in log-linear regression needs to be non-negative\n\n");
                        return(FALSE);
                    } else {
                        return(TRUE);
                    }

                } else if (link == "logit") {

                    if (min(x, na.rm = TRUE) < 0 | max(x, na.rm = TRUE) > 1) {
                        stop("\nThe response in logistic regression needs to be between 0 and 1\n\n");
                        return(FALSE);
                    } else if (any(x != 0 & x != 1, na.rm = TRUE)) {
                        warning("\nnon-integer response in logistic regression\n\n");
                        return(TRUE);
                    } else {
                        return(TRUE);
                    }

                } else {
                    return(TRUE);
                }

            } else {

                stop("\nThe response the regression needs to be a factor or a numeric vector\n\n");
                return(FALSE);

            }

        }


        if (missing(data)) {
            data <- parent.frame()
        }

        if (!missing(oformula)) {
            oterms <- terms(oformula)
            omf <- model.frame(formula = oformula, data = data, na.action = na.pass)

            ## we define the outcome as the response in the oformula
            if (attr(oterms,"response")) {
                outname <- all.vars(oformula)[1]
            } else {
                outname <- NULL
            }
        } else {

            if(estimation.method %in% c("dr", "o")) {
                stop("\nAn outcome nuisance model is needed\n\n");
            } else {
                oterms <- NULL;
                outname <- NULL;
            }
        }

        ## Extract the outcome if it is given
        if (!missing(outcome)) {

            ## If the outcome is not given as a string
            ## get the name of the object that was given
            ## as input
            if (is.character(outcome)) {
                outname <- outcome
            } else {
                outname <- as.character(call$outcome)
            }

            if (!missing(oformula)) {
                warning(paste("\nDuplicate specifications of the outcome, using ",
                              outname, "from the outcome argument\n\n") )
            }

        }



        if (is.null(outname)) {
            stop("An outcome needs to be specified\n\n")
        } else {
            tempof <- formula(paste("~", outname))

            if (is.environment(data)) {
                environment(tempof) <- data
            }

            yf <- model.frame(tempof, data = data, na.action = na.pass)
            y.all <- model.matrix(tempof, yf)

            if (dim(y.all)[2] != 2) {
                stop("\nThe outcome needs to be of exactly one dimension\n
or a factor with two levels\n\n")
            }

            y <- y.all[, 2, drop = F]

            ## Update the column names if the outcome was not numeric
            outname <- colnames(y)

            if (!response.ok(y[, 1, drop=TRUE], link = olink)) {
                stop("\nBad outcome!\n\n")
            }

            nobs <- nrow(y)

            ## Identify complete observations
            compl.rows <- !is.na(y[, 1])

        }

        ## Check that we have the same number of observations for all given
        ## objects
        if (!is.null(oterms) & !length(attr(oterms, "term.labels")) == 0) {
            if (nrow(omf) != nobs) {
                stop("The should be as many outcome as outcome nuisance model observations\n\n")
            } else {
                compl.rows <- compl.rows & !rowSums(is.na(omf))
            }
        }


        if (!missing(eformula)) {

            eterms <- terms(eformula)
            emf <- model.frame(formula = eformula, data = data, na.action = na.pass)

            if (!length(attr(eterms, "term.labels")) == 0 & nrow(emf) != nobs) {
                stop("The should be as many outcome as exposure nuisance model observations\n\n")
            } else {
                compl.rows <- compl.rows & !rowSums(is.na(emf))
            }

            if (attr(eterms,"response")) {
                expname <- all.vars(eformula)[1]
                ## If there is no response in eformula
                ## we cannot determine what the exposure is
            } else {
                expname <- NULL
            }

        } else {

            if(estimation.method %in% c("dr", "e")) {
                stop("\nAn exposure nuisance model is needed\n\n");
            } else {
                eterms <- NULL
                expname <- NULL
            }

        }

        ## Extract the exposure if it is given
        if (!missing(exposure)) {


            ## If the exposure is not given as a string
            ## get the name of the object that was given
            ## as input
            if (is.character(exposure)) {
                expname <- exposure
            } else {
                expname <- as.character(call$exposure)
            }

            if (!missing(eformula)) {
                warning(paste("\nDuplicate specifications of the exposure, using ",
                              expname, "from the exposure argument\n\n") )
            }

        }

        if (!is.null(expname)) {

            tempef <- formula(paste("~", expname))

            if (is.environment(data)) {
                environment(tempef) <- data
            }

            af <- model.frame(tempef, data = data, na.action = na.pass)
            a.all <- model.matrix(tempef, af)

            if (dim(a.all)[2] != 2) {
                stop("\nThe exposure needs to be of exactly one dimension\n
or a factor with two levels\n\n")
            }

            a <- a.all[, 2, drop = F]

            ## Update the column names if the outcome was not numeric
            expname <- colnames(a)

            if ( nrow(a) != nobs) {
                stop("The should be as many outcome as exposure observations\n\n")
            } else {
                compl.rows <- compl.rows & !is.na(as.vector(a))
            }

            if (estimation.method %in% c("dr", "e")) {
                if (!response.ok(a[, 1, drop=TRUE], link = elink) &
                    estimation.method != "o") {
                    stop("\nBad exposuree!\n")
                }
            }

        } else {
            a <- NULL
        }

        if (!missing(iaformula)) {

            iaterms <- terms(iaformula)
            iamf <- model.frame(formula = iaformula, data = data,
                                na.action = na.pass)

            if (length(attr(iaterms, "term.labels")) == 0) {
                iaterms <- NULL
            } else {
                if (nrow(iamf) != nobs) {
                    stop("The should be as many outcome as interaction observations\n\n")
                } else {
                    compl.rows <- compl.rows & !rowSums(is.na(iamf))
                }
            }

        } else {
            iaterms <- NULL
        }

        ## Get the clusterid if it is given and create a clusterid otherwise
        if (!missing(clusterid)) {

            if (is.character(clusterid)) {
                clustname <- clusterid
                if (is.environment(data)) {
                    id <- get(clustname, envir = data)
                } else if (is.data.frame(data)) {
                    id <- data[clustname]
                } else if (is.matrix(data)) {
                    id <- data[, clustname]
                } else {
                    stop(paste("The clusterid ", clusterid, " could not be found\n\n"))
                }

            } else {
                clustname <- call[["clusterid"]]
                id <- as.vector(clusterid)
            }

            if (is.list(id)) {
                id <- unlist(id)
            } else {
                id <- as.vector(id)
            }

            if (length(id) != nobs) {
                stop("Each observation needs to have exactly one observation of the clusterid\n\n")
            }

            compl.rows <- compl.rows & !is.na(id)

        } else {
            if (cond) {
                stop("For conditional methods, clusterid is required\n\n")
            } else {
                id <- 1:nrow(y)
                clustname <- "id"
            }
        }

        ## Make sure that the data is sorted by the cluster variable
        if (is.unsorted(id, na.rm = TRUE)) {
            obs.order <- order(id)
            ## The permutation of complete rows (boolean)
            compl.rows <- compl.rows[obs.order]
        } else {
            obs.order <- 1:nobs
        }

        ## For conditional logistic methods, we exclude non-informative clusters
        if (cond & olink == "logit") {
            if (estimation.method == "o") {

                id.tmp <- as.factor(id[obs.order])
                y.tmp <- y[obs.order, 1]
                ## For each cluster, identify the number of different values for the outcome
                ni.vals <- ave(as.vector(y.tmp), id.tmp, FUN = function(y) {length(unique(y[which(!is.na(y)) ]))})

                ## Only use outcome-discordant clusters
                compl.rows <- compl.rows & (ni.vals > 1)

            } else if (estimation.method == "e") {

                id.tmp <- as.factor(id[obs.order])
                a.tmp <- a[obs.order, 1]
                ## For each cluster, identify the number of different values for the outcome
                ni.vals <- ave(as.vector(a.tmp), id.tmp, FUN = function(a) {length(unique(a[which(!is.na(a))]))})

                ## Only use outcome-discordant clusters
                compl.rows <- compl.rows & (ni.vals > 1)

            }
        }


        nobsnew <- sum(compl.rows)

        ## Get the row numbers that we will use
        rows.to.use <- obs.order[which(compl.rows)]

        ## Extract observations for complete observations
        id <- as.factor(id[rows.to.use])
        y <- y[rows.to.use,, drop = F]

        if (!is.null(a)) {
            a <- a[rows.to.use,, drop = F]
        }

        ## Extract covariates
        if (!is.null(oterms)) {
            if (length(attr(oterms, "term.labels")) > 0) {
                v <- model.matrix(oterms, omf)[rows.to.use,, drop = F]
            } else {
                v <- matrix( rep(1, nobsnew), ncol = 1)
            }
        } else {
            v <- NULL
        }

        if (!is.null(eterms)) {
            if (length(attr(eterms, "term.labels")) > 0) {
                z <- model.matrix(eterms, emf)[rows.to.use,, drop = F]
            } else {
                z <- matrix( rep(1, nobsnew), ncol = 1)
            }
        } else {
            z <- NULL
        }

        x <- matrix( rep(1, nobsnew), ncol = 1)

        if (!is.null(iaterms)) {
            if (!length(attr(iaterms, "term.labels")) == 0) {
                x <- model.matrix(iaterms, iamf)[rows.to.use,, drop = F]
            }
        }

        if (olink == "logit") {
            yx <- x * as.vector(y)
            colnames(yx)[1] <- colnames(y)
            colnames(yx)[-1] <- paste(colnames(y), colnames(x)[-1], sep = ":")
        } else {
            yx <- NULL
        }

        if (!is.null(a)) {
            ax <- x * as.vector(a)
            colnames(ax)[1] <- colnames(a)
            colnames(ax)[-1] <- paste(colnames(a), colnames(x)[-1], sep = ":")
        } else {
            ax <- NULL
        }

        ## Do not use an intercept for conditional methods
        if (cond) {

            if (!is.null(oterms)) {
                if (length(attr(oterms, "term.labels"))>0 & attr(oterms,"intercept")) {
                    v <- v[, -1, drop = F]
                } else {
                    v <- NULL
                }
            }

            if (!is.null(eterms)) {
                if (length(attr(eterms, "term.labels"))>0 & attr(eterms,"intercept")) {
                    z <- z[, -1, drop = F]
                } else {
                    z <- NULL
                }
            }


        }

        drgee.data <- list(y = y, a = a, x = x, ax = ax, v = v, z = z, yx = yx,
                          id = id, idname = clustname,
                          outname = outname, expname = expname,
                          clustname = clustname,
                          olink = olink, elink = elink,
                          oterms = oterms, eterms = eterms,
                          cond = cond)

        class(drgee.data) <- "drgeeData"
        return(drgee.data)
    }

summary.drgeeData <-
    function(object, ...) {

        outcome <- colnames(object$y)
        exposure <- colnames(object$a)
        exposures.all <- colnames(object$ax)
        covariates <- setdiff( union(colnames(object$v),colnames(object$z)),  "(Intercept)" )
        interactions <- colnames(object$x)

        main.model <- paste(colnames(object$y), "~",
                            paste(colnames(object$ax), collapse = " + "), sep = " ")

        if (object$cond) {
            outcome.nuisance.model <- paste(colnames(object$y),"~",
                                            paste(c(colnames(object$v)),
                                                  collapse = " + "), sep = " ")
            outcome.model <- paste(colnames(object$y),"~",
                                   paste(c(colnames(object$ax), colnames(object$v)),
                                         collapse = " + "), sep = " ")

        } else {
            outcome.nuisance.model <- paste(colnames(object$y),"~",
                                            paste(c(colnames(object$v)[-1]),
                                                  collapse = " + "), sep = " ")
            outcome.model <- paste(colnames(object$y),"~",
                                   paste(c(colnames(object$ax), colnames(object$v)[-1]),
                                         collapse = " + "), sep = " ")

        }

        if (object$cond) {
            exposure.nuisance.model <- paste(colnames(object$a),"~",
                                             paste(colnames(object$z),
                                                   collapse = " + "), sep = " ")

            exposure.model <- paste(colnames(object$a),"~",
                                    paste(c(colnames(object$yx), colnames(object$z)),
                                          collapse = " + "), sep = " ")

        } else {
            exposure.nuisance.model <- paste(colnames(object$a),"~",
                                             paste(colnames(object$z)[-1],
                                                   collapse = " + "), sep = " ")

            exposure.model <- paste(colnames(object$a),"~",
                                    paste(c(colnames(object$yx), colnames(object$z)[-1]),
                                          collapse = " + "), sep = " ")

        }

        ans <- list(outcome = outcome,
                    exposure = exposure,
                    covariates = covariates,
                    interactions = interactions,
                    main.model = main.model,
                    outcome.nuisance.model = outcome.nuisance.model,
                    outcome.model = outcome.model,
                    exposure.nuisance.model = exposure.nuisance.model,
                    exposure.model = exposure.model,
                    n.obs = length(object$id),
                    n.clust = nlevels(object$id),
                    clustname = object$clustname,
                    olink = object$olink,
                    elink = object$elink,
                    cond = object$cond)

        class(ans) <- "summary.drgeeData"

        return(ans)
    }

print.summary.drgeeData <-
    function(x, digits = max(3L, getOption("digits") -
                    3L), ...) {

        cat("\nOutcome: ", colnames(x$outcome), "\nExposure: ", colnames(x$exposure),
            "\nInteractions: ", paste(colnames(x$interactions), collapse = ", "))

        cat("\n\nMain model: ", x$main.model, "\nwith link function: ", x$olink)

        cat("\n\nOutcome nuisance model: ", x$outcome.nuisance.model, "\nwith link function: ", x$olink)

        cat("\n\nExposure nuisance model: ", x$exposure.nuisance.model, "\nwith link function: ", x$elink, "\n")

        cat("\n\n", x$n.obs, "with ", x$n.clust, " clusters\n")
    }

