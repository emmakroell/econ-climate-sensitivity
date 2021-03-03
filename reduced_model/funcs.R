### Functions
# Phillips curve
# phill <- function(lambda, pars){
#   Phill <- pars[['phi1']] / ((1 - lambda)^2) - pars[['phi0']]
#   return(Phill)
# }
phill <- function(lambda, pars){
  Phill <- pars[['phi_0']] + pars[['phi_1']] * lambda
  return(Phill)
}

# # Linear investment function
inv <- function(x, pars) {
     Inv <- pmax(pars[['kappa_min']], pmin(pars[['kappa_const']] +
              pars[['kappa_slope']] * x, pars[['kappa_max']]))
     return(Inv)
   }


# inv <- function(profit_share, pars) {
#   Inv <- pars[['kappa0']] + exp(pars[['kappa1']] +
#                                   pars[['kappa2']] * profit_share)
#   return(Inv)
# }

# dividend function
div <- function(x, pars){
  Div <- pmax(pars[['div_min']], pmin(pars[['div_const']] +
                pars[['div_slope']] * x, pars[['div_max']]))
  return(Div)
}


# Inflation
infl <- function(omega, pars){
  Infl <-  pars[['eta_p']] * (pars[['markup']] * omega - 1)
  return(Infl)
}

## Workforce dynamics
pop_dyn <- function(N, pars) {
  d_N <- pars[['delta_N']]*(1 - N/pars[['N_max']])
  return(d_N)
}


transform_lambda <- function(lambda){
  lambda_trnfm <- tan(pi*(lambda - 0.5))
  return(lambda_trnfm)
}

transform_omega <- function(omega){
  omega_trnfm <- log(omega)
  return(omega_trnfm)
}

return_lambda <- function(lambda){
  lambda <- atan(lambda) / pi  + 0.5
  return(lambda)
}

return_omega <- function(omega){
  omega <- exp(omega)
  return(omega)
}

# =============================================================================
# need a list of functions for parallel computing:
functions <- c('phill',
               'inv',
               'div',
               'infl',
               'pop_dyn',
               'transform_lambda',
               'transform_omega',
               'return_lambda',
               'return_omega')

