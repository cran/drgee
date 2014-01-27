robVcov <-
    function(U, d.U, id=NULL){

        inv.d.U <- try(solve(d.U))

        if(class(inv.d.U)=='try-error'){

            vcov <- matrix(rep( NA, ncol(d.U)^2 ), ncol=ncol(d.U) )

        }else{
            # If data is not clustered we divide by
            # the number of observations
            if(is.null(id) | length(unique(id))==length(id)){
                correction.term <- 1 / nrow(U)

            # If data is clustered we use the clustersums of U
            # and divide by the number of clusters
            }else{
                U <- apply(U, 2, function(x) tapply(x, as.factor(id), sum))
                correction.term <- length(unique(id)) / length(id)^2
            }

            vcov <- inv.d.U %*% cov(U) %*% t(inv.d.U) * correction.term
        }

        return(vcov)
    }
