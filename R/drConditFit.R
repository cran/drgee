expit <- function(x){ (exp(-x) + 1)^(-1) }

drConditFit <- function(object, rootFinder = findRoots) {

    ## Find the total number of observations
    n.obs <- length(object$id)
    
    ## Find the number of main parameters
    n.beta.params <- ncol(object$ax)

    ## Find the number of outcome nuisance parameters
    n.gamma.params <- ncol(object$v)
    if( is.null(n.gamma.params) ) n.gamma.params <- 0
    n.oparams <- n.beta.params + n.gamma.params
    
    ## Find the number of outcome nuisance parameters
    n.alpha.params <- ncol(object$z)
    if( is.null(n.alpha.params) ) n.alpha.params <- 0
    n.eparams <- n.beta.params + n.alpha.params

    ## Clumsy workaround to get through the CRAN check
    idx <- NULL
    id <- NULL
    size <- NULL
    y <- NULL
    a <- NULL
    y.sum <- NULL
    a.sum <- NULL
    y.disc <- NULL
    a.disc <- NULL
    ya.disc <- NULL
    
    ya.dt <- data.table(idx = seq_len(n.obs),
                        id = object$id,
                        y = as.integer(object$y),
                        a = as.integer(object$a) )
    
    names(ya.dt) <- c("idx", "id", "y", "a")

    setkey(ya.dt, "id")
    
    ## Sum by cluster id
    ya.clust.dt <- ya.dt[, list( y.sum = sum(y), a.sum = sum(a), size = .N), by = id][, c("y.sum", "a.sum", "size"), with=F]

    ## Identify outcome discordant clusters
    ya.clust.dt[, y.disc := ifelse( y.sum < size & y.sum > 0, 1, 0)]
    ## Identify exposure discordant clusters
    ya.clust.dt[, a.disc := ifelse( a.sum < size & a.sum > 0, 1, 0)]
    ## Identify doubly discordant clusters
    ya.clust.dt[, ya.disc := y.disc * a.disc]

    ## Extract info about the number of outcome discordant,
    ## exposure discordant and doubly discordant clusters
    disc.info <- colSums( ya.clust.dt[, c("y.disc", "a.disc", "ya.disc"), with=F] )

    ## Find doubly discordant pairs
    ## all possible combinations within pairs
    ya.pairs.dt <- with( ya.dt,
                        merge(ya.dt,
                              ya.dt,
                              allow.cartesian = TRUE,
                              suffixes = c(".1", ".2")
                              )[(idx.1 != idx.2) & (y.1 != y.2) & (a.1 != a.2)]
                        )
    
    ## Find the ids
    id.dd.1 <- object$id[ya.pairs.dt$idx.1]
    ## For the first member of each pair, 
    y.dd.1 <- object$y[ya.pairs.dt$idx.1,, drop = FALSE]
    colnames(y.dd.1) <- object$y.names
    ## find the exposure
    a.dd.1 <- object$a[ya.pairs.dt$idx.1,, drop = FALSE]
    colnames(a.dd.1) <- object$a.names
    ## find the interaction terms
    x.dd.1 <- object$x[ya.pairs.dt$idx.1,, drop = FALSE]
    ## find the covariates in the outcome model
    v.dd.1 <- object$v[ya.pairs.dt$idx.1,, drop = FALSE]
    ## find the covariates in the exposure model
    z.dd.1 <- object$z[ya.pairs.dt$idx.1,, drop = FALSE]

    ## For the second member of each pair, 
    ## find the interaction terms
    x.dd.2 <- object$x[ya.pairs.dt$idx.2,, drop=FALSE]
    ## find the covariates in the outcome model
    v.dd.2 <- object$v[ya.pairs.dt$idx.2,, drop=FALSE]
    ## find the covariates in the exposure model
    z.dd.2 <- object$z[ya.pairs.dt$idx.2,, drop=FALSE]

    ## For each pair,
    ## find the sum of interaction terms
    x.dd.sum <- x.dd.1 + x.dd.2
    
    ## Find the independent variables in the prospective model
    a.1.x.sum.dd <- x.dd.sum * as.vector(a.dd.1)
    colnames(a.1.x.sum.dd) <- object$ax.names
    v.pseudo <- cbind( x.dd.2, (v.dd.1 - v.dd.2) )
    v.pseudo.names <- paste("on", 1:ncol(v.pseudo), sep="")
    colnames(v.pseudo) <- v.pseudo.names

    ## Find the independent variables in the retrospective model
    y.1.x.sum.dd <- x.dd.sum * as.vector(y.dd.1)
    colnames(y.1.x.sum.dd) <- object$yx.names
    z.pseudo <- cbind( x.dd.2, (z.dd.1 - z.dd.2) )
    z.pseudo.names <- paste("en", 1:ncol(z.pseudo), sep="")
    colnames(z.pseudo) <- z.pseudo.names

    ## find the original id for each pair
    id.dd <- object$id[ya.pairs.dt$idx.1]
    
    ## The number of pseudopairs
    disc.info <- c(disc.info,  pseudo.disc = nrow(ya.pairs.dt) )

    pseudopairs <- list(y = y.dd.1,
                        a = a.dd.1,
                        x = x.dd.sum,
                        ax = a.1.x.sum.dd, 
                        v = v.pseudo, 
                        z = z.pseudo, 
                        yx = y.1.x.sum.dd, 
                        id = id.dd,
                        clustname = object$clustname,
                        y.names = object$y.names, 
                        a.names = object$a.names, 
                        x.names = object$x.names, 
                        ax.names = object$ax.names, 
                        v.names = v.pseudo.names, 
                        z.names = z.pseudo.names, 
                        yx.names = object$yx.names, 
                        olink = "logit",
                        elink = "logit",
                        cond = TRUE,
                        oterms = NULL,
                        eterms = NULL)
    
    class(pseudopairs) <- "drgeeData"
                        
    pseudo.fit <- dreFit( object = pseudopairs )
    pseudo.fit$disc.info <- disc.info

    return( pseudo.fit )
    
}
