geeFit <-
function(y, x, link=c("identity","log","logit"), ...){

	link <- match.arg(link)

	if(link=="identity"){
		fit <- try( glm.fit(x=x, y=y, family=gaussian(), ...) )
	}else if(link=="log"){
		fit <- try( glm.fit(x=x, y=y, family=quasipoisson(), ...) )
	}else if(link=="logit"){
		fit <- try( glm.fit(x=x, y=y, family=binomial(), ...) )
	}

	if(class(fit)=='try-error'){

            d.res=matrix(rep(NA,nrow(x)*ncol(x)), nrow=nrow(x))
            colnames(d.res) <- colnames(x)

            return( list(coefficents=rep(NA, ncol(x)),
                         fitted.values=matrix(rep(NA,nrow(y)), ncol=ncol(y)),
                         res=rep(NA,nrow(y)), d.res=d.res) )

        }else{

            if(link=="identity"){
		d.res <- -x
            }else if(link=="log"){
		d.res <- -x * fit$fitted.values
            }else if(link=="logit"){
		d.res <- -x * fit$fitted.values / (1+exp(fit$linear.predictors))
            }

            colnames(d.res) <- colnames(x)

            return(list(coefficients = fit$coefficients, fitted.values =
                        fit$fitted.values, res = as.vector(y-fit$fitted.values),
                        d.res = d.res))
    }

}
