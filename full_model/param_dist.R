### PARAMETER DISTRIBUTIONS
library(rmutil)

# uncertainty parameters:
unc_parms <- c(
  eta_mean  = 0.4,
  eta_sd    = 0.12,
  # Labour productivity 
  alpha_mean  = 0.0206,
  alpha_sd    = 0.0112,
  # Equilibrium climate sensitivity
  S_mean      = 1.107,
  S_sd        = 0.264,
  # Pre-industrialisation CO2 concentration in the UP layer
  CO2_UP_preind_mean = 5.8855763,
  CO2_UP_preind_sd   = 0.2512867
)


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



# Parameter distributions
draw_eta <- function() {
  mu     <- parms[['eta_mean']] 
  sigma  <- parms[['eta_sd']] 
  draw <- rnorm(1, mu, sigma)
  while (draw < 0 | draw > 4 ){
    draw <- rnorm(1, mu, sigma)
  }
  return(draw)
}


draw_markup <- function() {
  pars <- get_ggamma(q=c(0.05,0.4,0.95),val=c(0.2,0.55,1.5),
                     control=list(maxit=2000),shift=0)
  draw <- rggamma(1,pars['s'],pars['m'],pars['f'])+1
  while (draw < 1 | draw > 3){
    draw <- rggamma(1,pars['s'],pars['m'],pars['f'])+1
  }
  return(draw)
}



draw_gamma <- function() {
  pars <- get_ggamma(q=c(0.05,0.4,0.95),val=c(0.1,0.09,1),
                     control=list(maxit=2000),shift=0)
  draw <- 1 - rggamma(1,pars['s'],pars['m'],pars['f'])
  while (draw < 0 | draw > 1){
    draw <- 1.05 - rggamma(1,pars['s'],pars['m'],pars['f'])
  }
  return(draw)
}


## Labour productivity, alpha, follows normal distribution
draw_alpha<- function() {
  mu     <- parms[['alpha_mean']] 
  sigma  <- parms[['alpha_sd']] 
  draw <- rnorm(1, mu, sigma)
  while (draw < 0 ){
    draw <- rnorm(1, mu, sigma)
  }
  return(draw)
}

## Equilibrium climate sensitivity, S, follows log-normal dist
draw_ECS <- function() {
  mu     <- parms[['S_mean']] 
  sigma  <- parms[['S_sd']] 
  return(rlnorm(1, mu, sigma))
}

## Pre-industriation CO2 concentration in the UP layer also follows
## log normal dist
draw_C_UP_preind <- function() {
  mu     <- parms[['CO2_UP_preind_mean']]
  sigma  <- parms[['CO2_UP_preind_sd']]
  return(rlnorm(1, mu, sigma))
}

# update parameters
Parms <- c(Parms,unc_parms)
