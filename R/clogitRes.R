clogitRes <-
    function(theta, y, x.cent, id) {

        ## Center covariates for numerical stability
        id.idx <- split(1:nrow(y),id)
        resid.list <- unlist(lapply( id.idx,
                                    function(idxs){clogitResCluster(theta, y[idxs,,drop = F],
                                                                    x.cent[idxs,,drop = F])} ), recursive=FALSE)
        ## Note that this requires that both x and y are sorted by id
        ## Otherwise we get the wrong order here

        d.eq.res.idx <- 2 * (1:length(id.idx))
        d.eq.res <- matrix( unlist(resid.list[d.eq.res.idx]), ncol=ncol(x.cent), byrow=TRUE)
        eq.res <- unlist(resid.list[-d.eq.res.idx])

	return( list(eq.res=eq.res, d.eq.res=d.eq.res) )

    }

clogitResCluster <-
    function(theta, yi, xi) {

	si <- sum(yi)
	ni <- nrow(yi)

        ## If more than half are cases
        ## we calculate the probability
        ## of being a non-case.
        ## Otherwise we calculate the
        ## probability of being a case

        if (si > 0.5*ni){
            p.case <- FALSE
            si <- ni - si
            xi <- -xi
        }else{
            p.case <- TRUE
        }

        nbetas <- ncol(xi)

  	## If all y.i are 0 or all y.i are 1, then return a 0 vector
  	if (si == 0 | si == ni) {

            return( list(eq.res = rep(0,ni),
                         d.eq.res = t( matrix(rep(0, ni*nbetas), nrow=ni) )
                         ) )

            ## If cluster sum si is =1
  	} else if (si == 1) {
            exp.theta.xi <- as.vector( exp(xi %*% theta) )
            sum.exp.theta.xi <- sum( exp.theta.xi )
            yi.hat <- as.vector( exp.theta.xi / sum.exp.theta.xi )

            if (p.case) {
                eq.res <- as.vector( yi - yi.hat )
                d.eq.res <- t( - xi * yi.hat + yi.hat %o% colSums(xi * exp.theta.xi) / sum(exp.theta.xi) )
            } else {
                eq.res <- as.vector( yi - ( rep(1, ni) - yi.hat ) )
                d.eq.res <- -t( - xi * yi.hat + yi.hat %o% colSums(xi * exp.theta.xi) / sum(exp.theta.xi) )
            }



            return( list(eq.res = eq.res, d.eq.res = d.eq.res ) )

            ## If cluster sum si is >1
  	} else {

            exp.theta.xi <- exp(xi %*% theta)

            ## Calculate the residuals

            ## Each element (m,j) of B is the sum of products of m elements exp.theta.xi
            ## such that the m elements are taken from exp.theta.xi[1:j]
            ## Recursion formula: B[m,j]n=B[m,j-1] + exp.theta.xi[j]*B[m-1,j-1]

            ## B[m,j] = 0 when m > j
            B <- matrix( rep(0,si * ni), nrow=si )

            ## Create the first row of B
            B[1,1] <- exp.theta.xi[1]

            for (j in 2:ni) {
                B[1,j] <- exp.theta.xi[j] + B[1,j-1, drop = F]
            }

            ## Generate columns 2 to ni for B
            for (k in 2:ni) {
                B[2:si,k] <- B[2:si,k-1, drop = F] + exp.theta.xi[k] * B[1:(si-1),k-1, drop = F]
            }

            ## Create a matrix of powers of -exp.theta.xi
            ## D[j,k]=(-exp.theta.xi[j])^k for k=1,...,si-1 and j=1,...,ni
            D <- matrix(rep(0, ni * si), nrow=ni)
            D[, 1] <- -exp.theta.xi
            for (j in 2:si) {
                D[, j] <- -exp.theta.xi * D[, j-1, drop = F]
            }

            R <- - D[, si] - D[, (si-1):1, drop = F] %*% B[1:(si-1), ni, drop = F]

            yi.hat <- R / B[si, ni]

            ## Calculate the derivative of the kernel

            ## We want a matrix of the derivatives of yi.hat
            ## with rows representing observations
            ## and columns the partial derivatives
            dyi.hat <- matrix(rep(0, ni * nbetas), nrow=ni)

            ## We calculate derivatives separately for each beta.q
            ## with respect to beta.q

            for (q in 1:nbetas) {

                ## Calculate derivatives of exp.theta.xi with respect to beta.q
                dq.exp.theta.xi <- xi[, q, drop = F] * exp.theta.xi

                dqB <- matrix( rep(0, si * ni), nrow = si)
                ## dqB[m,j] = 0 when m > j for all q

                ## Create the first row of dqB
                dqB[1, 1] <- dq.exp.theta.xi[1]
                for (k in 2:ni) {
                    dqB[1,k] <- dqB[1, k-1, drop = F] + dq.exp.theta.xi[k]
                }

                ## Generate columns 2 to ni dqB
                for (k in 2:ni) {
                    dqB[2:si, k] <- dqB[2:si, k-1, drop = F] +
                        dq.exp.theta.xi[k] * B[1:(si-1), k-1, drop = F] +
                            exp.theta.xi[k] * dqB[1:(si-1), k-1, drop = F]
                }

                dqD <- matrix(rep(0, ni * si), nrow = ni)
                dqD[, 1] <- - xi[, q] * exp.theta.xi
                for (p in 2:si) {
                    dqD[, p] <- p * xi[, q] * D[, p]
                }

                dqR <- - dqD[, p, drop = F]

                dqR <- dqR - dqD[, (si-1):1, drop = F] %*% B[1:(si-1), ni, drop = F]

                dqR <- dqR - D[, (si-1):1, drop = F] %*% dqB[1:(si-1), ni, drop = F]

                dyi.hat[, q] <- ( dqR * B[si, ni] - R * dqB[si, ni] ) /
                    B[si, ni]^2

            }

            if (p.case) {
                eq.res <- as.vector(yi - yi.hat)
                d.eq.res <- -t(dyi.hat)
            } else {
                eq.res <- as.vector(yi - (1 - yi.hat))
                d.eq.res <- t(dyi.hat)
            }

            return( list(eq.res = eq.res,
                         d.eq.res = d.eq.res
                         ) )

        }

    }

