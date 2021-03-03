# Distribution for pricing params from B. Bolker
## want distribution: lower limit = 1
## q0.05 ~ 1.05
## q0.5 ~ 1.25
## q0.95 ~ 2

## Johnson distribution?

## three-parameter distribution
library(rmutil)
## generalized gamma

get_ggamma <- function(q=c(0.05,0.5,0.95),val=c(1.05,1.25,2),shift=1,...) {
  obj <- function(p) {
    ## generalized gamma parameterized by shape (s), scale (m), family (f)
    sum((qggamma(q,s=p[1],m=p[2],f=p[3])-(val-shift))^2)
  }
  opt <- optim(fn=obj,par=c(s=1,m=1,f=1),...)
  if (opt$convergence!=0) warning("didn't converge")
  return(opt$par)
}

pars <- get_ggamma(control=list(maxit=2000))
# with(as.list(pars),curve(dggamma(x-1,s,m,f),from=1.0001,to=3,
#                          main=expression(paste("Distribution of the markup rate, ", xi))))
#with(as.list(pars),abline(v=qggamma(c(0.05,0.5,0.95),s,m,f)+1,lty=2))

# gamma
# pars <- get_ggamma(q=c(0.05,0.5,0.95),val=c(0.1,0.2,1),control=list(maxit=2000),shift=0)
# with(as.list(pars),curve(dggamma(1.1-x,s,m,f),from=0,to=1.0999))

pars <- get_ggamma(q=c(0.05,0.4,0.95),val=c(0.1,0.2,1),control=list(maxit=2000),shift=0)
# with(as.list(pars),curve(dggamma(1.05-x,s,m,f),from=0,to=1.04,
#                          main=expression(paste("Distribution of ", gamma))))
