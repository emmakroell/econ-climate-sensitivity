### Functions
#=============================================================================
# Phillips curve
# phill <- function(lambda, pars){
#   Phill <- pars[['phi1']] / ((1 - lambda)^2) - pars[['phi0']]
#   return(Phill)
# }

phill <- function(lambda, pars){
  Phill <- pars[['phi_0']] + pars[['phi_1']] * lambda
  return(Phill)
}

# Investment function
# # Linear investment function
inv <- function(x, pars) {
  Inv <- pmax(pars[['kappa_min']], pmin(pars[['kappa_const']] +
                                          pars[['kappa_slope']] * x, pars[['kappa_max']]))
  return(Inv)
}

## Workforce dynamics
pop_dyn <- function(N, pars) {
  d_N <- pars[['delta_N']]*(1 - N/pars[['N_max']])
  return(d_N)
}

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


#=============================================================================
## CLIMATE MODULE
# Exogenous forcing function
F_exo <- function(t, pars) {
  current_time <- t + pars[['F_exo_start_year']]
  future <- pars[['F_exo_end_year']] - pars[['F_exo_start_year']]
  if(as.numeric(current_time) <= pars[['F_exo_end_year']]) {
    F_exo_slope <- (pars[['F_exo_end']] - pars[['F_exo_start']])/future
    F_exo <- pars[['F_exo_start']] + t*F_exo_slope
  } else {
    F_exo <- pars[['F_exo_end']]
  }
  return(F_exo)
}

# Exogenous forcing long function
F_exo_vec <- function(t, pars) {
  n_steps <- length(t)
  F_exo <- rep(0, n_steps)
  for (k in 1:n_steps) {
    F_exo[k] <- F_exo(t = t[k], pars = pars)
  }
  return(F_exo)
}


dam_curve <- function(Temp, pars, options) {
  switch(options$damage_scenario,
         'None' = 0 * Temp,
         'Nordhaus' = 1 - 1/(1 + pars[['dam_2']] * Temp^2),
         '10at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_ten']] * Temp^pars[['dam_exp']]),
         '20at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_twenty']] * Temp^pars[['dam_exp']]),
         '30at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_thirty']] * Temp^pars[['dam_exp']]),
         '40at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_forty']] * Temp^pars[['dam_exp']]),
         '50at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_fifty']] * Temp^pars[['dam_exp']])
  )
}

# Carbon pricing (three steps)
carbon_price <- function(p_Car,pars,options,time,state){
  # the price of carbon cannot be greater than the price of the backstop technology
  if(p_Car >= state[['p_BS']]) {
    p_Car = state[['p_BS']]
  }
  
  # Stern-Stiglitz carbon pricing scheme:
  if (options$p_carb_scheme == 'Stern-Stiglitz'){
    if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start']])) {
      d_p_Car <- (pars[['p_Car_val_1']] - pars[['p_Car_start']]) /
        (pars[['p_Car_step_year_1']] - time[['start']])
    } else {
      d_p_Car  <- (pars[['p_Car_val_2']] - pars[['p_Car_val_1']]) /
        (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
    }
    
    # No carbon pricing
  }  else if (options$p_carb_scheme == 'None'){
    d_p_Car   <- pars[['g_p_BS']] * p_Car
  }
  return(d_p_Car)
}

# need a list of functions for parallel computing:
functions <- c('phill',
               'inv',
               'div',
               'infl',
               'pop_dyn',
               'transform_lambda',
               'transform_omega',
               'return_lambda',
               'return_omega',
               'F_exo',
               'F_exo_vec',
               'dam_curve',
               'carbon_price')
