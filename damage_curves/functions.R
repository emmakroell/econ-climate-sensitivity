### Functions
# Old functions file from 2019 to work with damage curve Rmd file (damage_curves.Rmd)
#=============================================================================
## MACROECONOMIC MODULE
# Production function
f_prod <- function(K, pars) {
  Y <- K/pars[['nu']]
  return(Y)
}

# Phillips curve
f_phill_N <- function(lambda, pars, active){
  switch(active$Phillips,
         'linear' = pars[['Phill_slope_N']] * lambda + pars[['Phill_lin_const_N']],
         'asymptotic' = pars[['Phill_denom_N']] / ((1 - lambda)^2) + pars[['Phill_asymp_const_N']]
  )
}

f_phill_S <- function(lambda, pars, active){
  switch(active$Phillips,
         'linear' = pars[['Phill_slope_S']] * lambda + pars[['Phill_lin_const_S']],
         'asymptotic' = pars[['Phill_denom_S']] / ((1 - lambda)^2) + pars[['Phill_asymp_const_S']]
  )
}

# Linear dividend function
f_div <- function(pi, pars) {
  div <- pmax(pars[['div_min']], pmin(pars[['div_max']], pars[['div_constant']] + pars[['div_slope']]*pi))
  return(div)
}

# Linear investment function
f_inv <- function(pi, pars, active) {
  switch(active$Invest,
         'linear' = pmax(pars[['I_lin_min']], pmin(pars[['I_lin_max']], pars[['I_lin_constant']]+pars[['I_lin_slope']]*pi)),
         'arctan' = pars[['kappa_0']] + pars[['kappa_1']] * atan(pars[['kappa_2']] * pi + pars[['kappa_3']])
  )
}

# Taylor rule function
f_taylor <- function(infl, pars) {
  taylor <- pmax(0, pars[['r_longtarget']] - pars[['r_taylor']] * 
                   pars[['i_target']] + (1 + pars[['r_taylor']]) * infl)
  return(taylor)
}

# Workforce dynamics function
f_workforce <- function(N, delta, max) {
  d_N <- delta * N * (1 - N/max)
  return(d_N)
}

## Labor productivity dynamics
f_labour_prod_N <- function(a, alpha, Temp, pars, active) {
  switch(active$Lab_Prod_N,
         "constant" = a * pars[['alpha_N']],
         "Burke"    = a * (pars[['alpha_1']] * (pars[['T_preind']] + Temp) + 
                             pars[['alpha_2']] * (pars[['T_preind']] + Temp)^2)
  )
}

f_labour_prod_S <- function(a, Temp, pars, active) {
  switch(active$Lab_Prod_S,
         "constant" = a * pars[['alpha_S']],
         "Burke"    = a * (pars[['alpha_1']] * (pars[['T_preind']] + Temp) + 
                             pars[['alpha_2']] * (pars[['T_preind']] + Temp)^2)
  )
}

#=============================================================================
## CLIMATE MODULE
# Exogenous forcing function
f_F_exo <- function(t, pars) {
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
f_F_exo_vec <- function(t, pars) {
  n_steps <- length(t)
  F_exo <- rep(0, n_steps)
  for (k in 1:n_steps) {
    F_exo[k] <- f_F_exo(t = t[k], pars = pars)
  }
  return(F_exo)
}

# Industrial radiative forcing function
f_F_ind <- function(CO2_AT, pars) {
  F_ind <- (pars[['F2xCO2']]/log(2))*log(CO2_AT/pars[['CO2_AT_preind']])
  return(F_ind)
}

#==============================================================================
# DAMAGES
# Climate Damage Functions
f_Dam_N <- function(Temp, pars, active) {
  switch(active$Damage_N,
         'None' = 0 * Temp,
         'Nordhaus' = 1 - 1/(1 + pars[['dam_1']] * Temp + pars[['dam_2']] * Temp^2),
         'Weitzman' = 1 - 1/(1 + pars[['dam_1']] * Temp + pars[['dam_2']] * Temp^2 +
                               pars[['dam_weitzman']] * Temp^pars[['dam_exp']]),
         'Stern' = 1 - 1/(1 + pars[['dam_1']] * Temp + pars[['dam_2']] * Temp^2 +
                            pars[['dam_stern']]*Temp^pars[['dam_exp']]),
         '10at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_ten']] * Temp^7),
         '20at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_twenty']] * Temp^7),
         '30at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_thirty']] * Temp^7),
         '40at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_forty']] * Temp^7),
         '50at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_fifty']] * Temp^7)
  )
}

f_Dam_S <- function(Temp, pars, active) {
  switch(active$Damage_S,
         'None' = 0 * Temp,
         'Nordhaus' = 1 - 1/(1 + pars[['dam_1']] * Temp + pars[['dam_2']] * Temp^2),
         'Weitzman' = 1 - 1/(1 + pars[['dam_1']] * Temp + pars[['dam_2']] * Temp^2 +
                               pars[['dam_weitzman']] * Temp^pars[['dam_exp']]),
         'Stern' = 1 - 1/(1 + pars[['dam_1']] * Temp + pars[['dam_2']] * Temp^2 +
                            pars[['dam_stern']]*Temp^pars[['dam_exp']]),
         '10at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_ten']] * Temp^7),
         '20at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_twenty']] * Temp^7),
         '30at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_thirty']] * Temp^7),
         '40at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_forty']] * Temp^7),
         '50at4' = 1 - 1/(1 + pars[['dam_2']] * Temp^2 + pars[['dam_fifty']] * Temp^7)
  )
}

# Damages allocated to output
f_Dam_Y <- function(Dam, pars, active) {
  D_Y <- Dam * (1 - active$Damage_K)
  return(D_Y)
}

# Damages allocated to capital
f_Dam_K <- function(Dam, pars, active) {
  D_K <- active$Damage_K * Dam
  return(D_K)
}

#==============================================================================
# Carbon pricing (three steps)
f_p_Car_N <- function(p_Car,pars,active,time,state){
  # the price of carbon cannot be greater than the price of the backstop technology
  if(p_Car >= state[['p_BS']]) {
    p_Car = state[['p_BS']]
  }
  # first, check that carbon pricing has started
  if(state[['counter']] >= (time[['start_p_Car_N']] - time[['start']])){
    # set emissions reduction rate
    n_N         <- min((p_Car  / ((1 - active$s_A_N) * state[['p_BS']]))^
                         (1 / (pars[['theta_2']] - 1)), 1)
    
    # Low carbon pricing scheme:
    if (active$p_Car_N == 'Low'){
      if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start_p_Car_N']])) {
        d_p_Car_N    <- (pars[['p_Car_L_val_1']] - pars[['p_Car_N_start']]) /
          (pars[['p_Car_step_year_1']] - time[['start_p_Car_N']])
      } else if (state[['counter']] <= (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])){
        d_p_Car_N  <- (pars[['p_Car_L_val_2']] - pars[['p_Car_L_val_1']]) /
          (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
      } else {
        d_p_Car_N <- (pars[['p_Car_L_val_3']] - pars[['p_Car_L_val_2']]) /
          (time[['end']] - pars[['p_Car_step_year_2']])
      }
      
      # Medium carbon pricing scheme:
    } else if (active$p_Car_N == 'Medium'){
      if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start_p_Car_N']])) {
        d_p_Car_N    <- (pars[['p_Car_M_val_1']] - pars[['p_Car_N_start']]) /
          (pars[['p_Car_step_year_1']] - time[['start_p_Car_N']])
      } else if (state[['counter']] <= (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])){
        d_p_Car_N  <- (pars[['p_Car_M_val_2']] - pars[['p_Car_M_val_1']]) /
          (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
      } else {
        d_p_Car_N <- (pars[['p_Car_M_val_3']] - pars[['p_Car_M_val_2']]) /
          (time[['end']] - pars[['p_Car_step_year_2']])
      }
      
      # High carbon pricing scheme:
    }  else if (active$p_Car_N == 'High'){
      if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start_p_Car_N']])) {
        d_p_Car_N    <- (pars[['p_Car_H_val_1']] - pars[['p_Car_N_start']]) /
          (pars[['p_Car_step_year_1']] - time[['start_p_Car_N']])
      } else if (state[['counter']] <= (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])){
        d_p_Car_N  <- (pars[['p_Car_H_val_2']] - pars[['p_Car_H_val_1']]) /
          (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
      } else {
        d_p_Car_N <- (pars[['p_Car_H_val_3']] - pars[['p_Car_H_val_2']]) /
          (time[['end']] - pars[['p_Car_step_year_2']])
      } 
      
      # Canadian federal backstop carbon pricing scheme
    } else if (active$p_Car_N == 'Canada'){
      if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start_p_Car_N']])) {
        d_p_Car_N    <- (pars[['p_Car_C_val_1']] - pars[['p_Car_N_start']]) /
          (pars[['p_Car_step_year_1']] - time['start_p_Car_N'])
      } else if (state[['counter']] <= (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])){
        d_p_Car_N  <- (pars[['p_Car_C_val_2']] - pars[['p_Car_C_val_1']]) /
          (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
      } else {
        d_p_Car_N <- (pars[['p_Car_C_val_3']] - pars[['p_Car_C_val_2']]) /
          (time[['end']] - pars[['p_Car_step_year_2']])
      }
      
      # No carbon pricing
    }  else if (active$p_Car_N == 'None'){
      d_p_Car_N   <- pars[['g_p_BS']] * p_Car
      
      # if carbon pricing hasn't started yet
    }} else {
      n_N         <- pars[['n_N_init']]
      d_p_Car_N   <- pars[['g_p_BS']] * p_Car
    }
  # return emissions reduction rate and carbon pricing gradient
  res <- data.frame(
    n_N = n_N,
    d_p_Car_N = d_p_Car_N
  )
  return(res)
}


f_p_Car_S <- function(p_Car,pars,active,time,state){
  # the price of carbon cannot be greater than the price of the backstop technology
  if(p_Car >= state[['p_BS']]) {
    p_Car = state[['p_BS']]
  }
  # first, check that carbon pricing has started
  if(state[['counter']] >= (time[['start_p_Car_S']] - time[['start']])){
    # set emissions reduction rate
    n_S         <- min((p_Car  / ((1 - active$s_A_S) * state[['p_BS']]))^
                         (1 / (pars[['theta_2']] - 1)), 1)
    
    # Low carbon pricing scheme
    if (active$p_Car_S == 'Low'){
      if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start_p_Car_S']])) {
        d_p_Car_S    <- (pars[['p_Car_L_val_1']] - pars[['p_Car_S_start']]) /
          (pars[['p_Car_step_year_1']] - time[['start_p_Car_S']])
      } else if (state[['counter']] <= (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])){
        d_p_Car_S  <- (pars[['p_Car_L_val_2']] - pars[['p_Car_L_val_1']]) /
          (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
      } else {
        d_p_Car_S <- (pars[['p_Car_L_val_3']] - pars[['p_Car_L_val_2']]) /
          (time[['end']] - pars[['p_Car_step_year_2']])
      }
      
      # Medium carbon pricing scheme
    } else if (active$p_Car_S == 'Medium'){
      if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start_p_Car_S']])) {
        d_p_Car_S    <- (pars[['p_Car_M_val_1']] - pars[['p_Car_S_start']]) /
          (pars[['p_Car_step_year_1']] - time[['start_p_Car_S']])
      } else if (state[['counter']] <= (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])){
        d_p_Car_S  <- (pars[['p_Car_M_val_2']] - pars[['p_Car_M_val_1']]) /
          (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
      } else {
        d_p_Car_S <- (pars[['p_Car_M_val_3']] - pars[['p_Car_M_val_2']]) /
          (time[['end']] - pars[['p_Car_step_year_2']])
      }
      
      # High carbon pricing scheme
    }  else if (active$p_Car_S == 'High'){
      if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start_p_Car_S']])) {
        d_p_Car_S    <- (pars[['p_Car_H_val_1']] - pars[['p_Car_S_start']]) /
          (pars[['p_Car_step_year_1']] - time[['start_p_Car_S']])
      } else if (state[['counter']] <= (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])){
        d_p_Car_S  <- (pars[['p_Car_H_val_2']] - pars[['p_Car_H_val_1']]) /
          (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
      } else {
        d_p_Car_S <- (pars[['p_Car_H_val_3']] - pars[['p_Car_H_val_2']]) /
          (time[['end']] - pars[['p_Car_step_year_2']])
      } 
      
      # Canadian federal backstop carbon pricing scheme
    } else if (active$p_Car_S == 'Canada'){
      if (state[['counter']] <= (pars[['p_Car_step_year_1']] - time[['start_p_Car_S']])) {
        d_p_Car_S    <- (pars[['p_Car_C_val_1']] - pars[['p_Car_S_start']]) /
          (pars[['p_Car_step_year_1']] - time[['start_p_Car_S']])
      } else if (state[['counter']] <= (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])){
        d_p_Car_S  <- (pars[['p_Car_C_val_2']] - pars[['p_Car_C_val_1']]) /
          (pars[['p_Car_step_year_2']] - pars[['p_Car_step_year_1']])
      } else {
        d_p_Car_S <- (pars[['p_Car_C_val_3']] - pars[['p_Car_C_val_2']]) /
          (time[['end']] - pars[['p_Car_step_year_2']])
      }
      
      # no carbon pricing
    }  else if (active$p_Car_S == 'None'){
      d_p_Car_S   <- pars[['g_p_BS']] * p_Car
      
      # if the carbon price is not yet in effect
    }} else {
      n_S         <- pars[['n_S_init']]
      d_p_Car_S   <- pars[['g_p_BS']] * p_Car
    }
  
  # return emissions reduction rate and carbon price gradient
  res <- data.frame(
    n_S = n_S,
    d_p_Car_S = d_p_Car_S
  )
  return(res)
}
