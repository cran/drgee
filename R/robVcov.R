robVcov <-
    function(u, d.u, id = NULL) {

        inv.d.u <- try(solve(d.u))

        if (class(inv.d.u) == 'try-error') {

            vcov <- matrix(rep( NA, ncol(d.u)^2 ), ncol = ncol(d.u) )

        } else {
            ## If data is not clustered we divide by
            ## the number of observations
            if (is.null(id) | length(unique(id)) == length(id)) {
                correction.term <- 1 / nrow(u)
            } else {
                ## If data is clustered we use the clustersums of u
                ## and divide by the number of clusters
                u <- apply(u, 2, function(x) tapply(x, as.factor(id), sum))
                correction.term <- length(unique(id)) / length(id)^2
            }

            vcov <- inv.d.u %*% cov(u) %*% t(inv.d.u) * correction.term
        }

        return(vcov)
    }
