condit <-
    function(y, x.cent, id) {

        ## Assuming that the observations are sorted by id
        
        ## Find the total number of observations
        n.obs <- length(id)

        ## Find the number of parameters
        n.params <- ncol(x.cent)

        ## Indices in the original matrix
        idx <- seq_len(n.obs)

        y.dt <- as.data.table(list(idx = idx, id = id, y = y))
        setkey(y.dt, id)

        ## Find cluster sizes
        clust.size <- y.dt[, .N, by = id][[2]]

        ## Find cluster sums
        ysum <- y.dt[, sum(y), by = id][[2]]

        ## Find cluster sums
        min.idx <- y.dt[, min(idx), by = id][[2]]

        ## Find number of clusters
        n.clust <- length(ysum)
        
        ## Find concordant clusters
        disc <- rep(TRUE, n.clust)
        disc[ysum == 0 | ysum == clust.size] <- FALSE

        ## Find cluster with more than half being cases
        inv <- integer(n.clust)
        inv[ysum > 0.5 * clust.size] <- 1L

        ysum.disc <- ysum[disc]
        clust.size.disc <- clust.size[disc]
        min.idx.disc <- min.idx[disc]
        inv.disc <- inv[disc]

        ##########################################################
        ## Fit using clogit
        ##########################################################

        fit.clogit <- try(survival::clogit(y~x.cent + strata(id),
                                           method = "exact",
                                           subset = id %in% id[min.idx.disc] ) )

        if (class(fit.clogit)[1] == "try-error") {

            beta.hat <- rep(NA, n.params)

            eq.res <- rep(NA, n.obs)

            d.eq.res <- matrix( rep(NA, n.obs*n.params), nrow = n.obs )

            naive.var = NULL

        } else {

            beta.hat <- coef(fit.clogit)

            resids <- .Call("conditRes", beta.hat, ysum.disc,
                            clust.size.disc, min.idx.disc, inv.disc,
                            y, x.cent, PACKAGE = "drgee")

            eq.res <- as.vector(resids$res)

            d.eq.res <- resids$dres

            naive.var = vcov(fit.clogit)

        }

        return( list(coefficients = beta.hat,
                     res = eq.res,
                     d.res = d.eq.res,
                     eq.x = x.cent,
                     optim.object = NULL,
                     naive.var = naive.var) )

    }
