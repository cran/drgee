findRoots <-
    function(beta.init, eq.func, d.eq.func = NULL, arg.list, ...){

        optim.object <- do.call(nleqslv, c(list(x = beta.init,
                                                fn = eq.func,
                                                jac = d.eq.func,
                                                arg.list = arg.list),
                                           list(...)))

        beta.hat <- optim.object$x
        
        if (optim.object$termcd > 2) {
            warning(paste("\nnleqslv: ", optim.object$message))
        }

        return( list(roots = beta.hat, optim.object = optim.object))

    }
