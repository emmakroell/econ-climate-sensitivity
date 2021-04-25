simulation <- function(time, init_state, parms, options = options,
                       method = 'lsoda', montecarlo = F) {
  # timer
  start_time <- proc.time()
  # time sequence for ode
  time_seq <- seq(0, time[['end']] - time[['start']], by = time[['step']])
  
  ### INITIAL STATE ####
  # convert some initial conditions in standard state variables to initial
  # conditions for the more intuitive differential equations
  state_0 <- function(IS = init_state, params = parms){
    with(as.list(c(IS, params)), {
      # initial carbon price (eq. 26)
      p_Car  <- n^(theta - 1) * (1 - options$subsidy) * p_BS
      # initial sigma (eq. 15)
      sigma  <- E_ind / ((1 - n) * K / nu)
      
      # create vector of initial conditions for numerical integrator
      state       <- c(
        K       = K,
        w       = w,
        debt    = debt,
        a       = a,
        p       = p,
        pop     = pop,
        sigma   = sigma,
        g_sigma = g_sigma,
        E_land  = E_land,
        CO2_AT  = CO2_AT,
        CO2_UP  = CO2_UP,
        CO2_LO  = CO2_LO,
        Temp    = Temp,
        Temp_LO = Temp_LO,
        p_BS    = p_BS,
        p_Car   = p_Car, 
        counter = 0
      )
      return(state)
    })
  }
  
  init_vars <- state_0()
  
  # update parameters
  C     <- 5/(1/parms['C_init'] + parms['beta_C']*(parms['S'] - parms['S_ref']))
  
  parms <- c(parms,
             p_Car_start = init_vars[['p_Car']], 
             F_exo_start_year = time[['start']],
             F_exo_end_year = time[['end']],
             n_init = init_state[['n']],
             C = C) 
  
  
  #============================================================================
  ## DIFFERENTIAL SYSTEM
  # Function for deSolve:
  Keen_w_prices <- function(t, state, parms) {
    with(as.list(c(parms,state)), {
      
      ### Auxiliary equations
      # Total output 
      Y_0 <- K / nu
      # Labour 
      L <- Y_0 / a
      # Damages from climate change 
      Dam <- dam_curve(Temp, parms, options)
      
      # emisssions reduction rate 
      n <- min((p_Car/((1 - options$subsidy) * p_BS))^(1/(theta - 1)), 1)
      # Abatement costs 
      A  <-  (sigma * p_BS / (1000 * theta)) * n^theta  
      # Productive output 
      Y <- (1-Dam)*(1-A)*Y_0
      # Industrial emissions 
      E_ind  <- sigma * (1 - n) * Y_0 
      # Total emissions
      E  <- E_ind + E_land 
      # Carbon tax
      Tax <- (p_Car * conv10to15 * E_ind)
      # Subsidy 
      Subsidy <- options$subsidy * A * Y_0   
      # Net transfer 
      Net_transfer <- Subsidy - Tax
      # Profits 
      Profits <- p * Y - w * L - r * debt + p * Net_transfer
      # profit share
      profit_share <- Profits / (p * Y)
      # wage share 
      omega <- w * L / (p * Y)
      # employment rate
      lambda <- L / pop
      # inflation
      i  <- infl(omega,parms)
      # investment 
      investmt <- inv(profit_share,parms) * Y
      # dividends
      dividends <- div(profit_share,parms) * p * Y
      # Phillips curve function
      phill <- phill(lambda,parms)
      
      ### Differential system 
      # capital 
      d_K <- investmt - delta * K
      # labour productivity
      d_a <- alpha * a
      # population
      d_pop <- pop * pop_dyn(pop,parms) #beta * pop
      # debt 
      d_debt <- p * investmt - Profits + dividends
      # wages 
      d_w <- w * (phill + gamma * i)
      # prices 
      d_p <- i * p
      # Emissions reduction 
      d_sigma     <- g_sigma * sigma
      d_g_sigma   <- delta_g_sigma * g_sigma
      # land based emissions gradient 
      d_E_land <- delta_E_land * E_land 
      # Backstop technology price
      d_p_BS      <- g_p_BS * p_BS
      # Carbon price 
      d_p_Car <- carbon_price(p_Car,parms,options,time,state)
      ### Climate Module
      # CO2 accumulation
      C02         <- matrix(c(CO2_AT,CO2_UP,CO2_LO),3,1)
      # Emissions
      Emissions   <- matrix(c(E/3.666,0,0),3,1)
      
      # Carbon cycle
      Diffusion <- matrix(c(-phi_12, phi_12 * (CO2_AT_preind / CO2_UP_preind),
                            0, phi_12, -phi_12 *(CO2_AT_preind / CO2_UP_preind) - 
                              phi_23, phi_23 * (CO2_UP_preind / CO2_LO_preind),
                            0, phi_23, -phi_23 * (CO2_UP_preind /CO2_LO_preind)),
                          3,3,byrow = TRUE)
      
      # compute C02 gradients
      d_C02    <- Emissions + Diffusion %*% C02
      d_CO2_AT <- d_C02[1,1]
      d_CO2_UP <- d_C02[2,1]
      d_CO2_LO <- d_C02[3,1]
      
      # Radiative forcing
      #browser()
      F_exo   <- F_exo(t = counter, pars = parms)
      # Industrial forcing 
      F_ind   <- (F2xCO2 / log(2)) * log(CO2_AT / CO2_AT_preind)
      # Total forcing 
      Forcing <- F_ind + F_exo 
      
      # Temperature change
      d_Temp <- Forcing/C - F2xCO2/(C * S) * Temp -  (gamma_star / C) *
        (Temp - Temp_LO)
      d_Temp_LO <- (gamma_star / C0) * (Temp - Temp_LO)
      
      # counter
      d_counter <- 1
      
      grad <- c(
        d_K       =  d_K ,
        d_w       =  d_w, 
        d_debt    =  d_debt,
        d_a       =  d_a,
        d_p       =  d_p,
        d_pop     =  d_pop,
        d_sigma   = d_sigma,
        d_g_sigma = d_g_sigma,
        d_E_land  = d_E_land,
        d_CO2_AT  = d_CO2_AT,
        d_CO2_UP  = d_CO2_UP,
        d_CO2_LO  = d_CO2_LO,
        d_Temp    = d_Temp,
        d_Temp_LO = d_Temp_LO,
        d_p_BS    = d_p_BS,
        d_p_Car   =  d_p_Car,
        d_counter =  d_counter
      )
      
      # ===========================================================================
      ## Return the list of gradients 
      return(list(grad))
    })
  }
  
  #==============================================================================
  ## SIMULATION
  simulation <- as.data.frame(ode(func = Keen_w_prices, 
                                  method = method, 
                                  y = init_vars, 
                                  parms = c(parms,time), 
                                  times = time_seq))
  
  # ============================================================================
  # REMAINING VARIABLES COMPUTATION
  
  compute.remaining <- function(sim = simulation, params = c(parms,time)){
    with(as.list(c(sim,params)), {
      ### Auxiliary equations
      # Total output (eq. 1)
      Y_0 <- K / nu
      # Damages from climate change
      Dam <- dam_curve(Temp, parms, options)
      # emissions reduction rate
      n <- (pmin((p_Car/((1 - options$subsidy) * p_BS)) ^
                   (1 / (theta - 1)), 1) - n_init) + n_init
      # Abatement costs
      A  <-  (sigma * p_BS / (1000 * theta)) * n^theta  
      # Remaining output
      Y <- (1-Dam)*(1-A)*Y_0
      # Labour
      L <- Y / a
      # Profits 
      Profits <- p * Y - w * L - r * debt
      # profit share
      profit_share <- Profits / (p * Y)
      # wage share 
      omega <- w * L / (p * Y)
      # debt share
      debt_share <- debt / (p * Y)
      # employment rate
      lambda <- L / pop
      # inflation
      i  <- infl(omega,parms)
      # investment 
      investmt <- inv(profit_share,parms) * Y
      # dividends
      dividends <- div(profit_share,parms) * p * Y
      # Phillips curve function
      phill <- phill(lambda,parms)
      
      # Emissions 
      E_ind  <- sigma * (1 - n) * Y / ((1 - Dam) * (1 - A))   # industrial emissions North
      E  <- E_ind + E_land  # total emissions
      
      
      # Radiative forcing
      F_exo   <- F_exo_vec(t = counter, pars = parms)
      F_ind   <- (F2xCO2 / log(2)) * log(CO2_AT / CO2_AT_preind)
      Forcing <- F_ind + F_exo 
      
      if(montecarlo == T){
        indices <- which((start + simulation$time) %in% seq(2016,2100,1))
        # 2016 start is hard coded (prob bad)
        sim_output <- data.frame(
          Y        = Y[indices],
          omega    = omega[indices],
          lambda   = lambda[indices],
          profit_share = profit_share[indices],
          i        = i[indices],
          p        = simulation$p[indices],
          debt_share = debt_share[indices], 
          K        = simulation$K[indices],
          E        = E[indices],
          Temp     = simulation$Temp[indices],
          p_Car      = simulation$p_Car[indices]
        )
        
      }else{
        
        # return results as a data frame
        sim_output <- data.frame(
          year     = start + time,
          counter = counter,
          K = K,
          w = w,
          p = p,
          debt = debt,
          a = a,
          pop = pop,
          lambda   = lambda,
          omega    = omega,
          debt_share = debt_share,
          profit_share = profit_share,
          dividends = dividends,
          i        = i,
          p_Car = p_Car,
          p_BS = p_BS,
          Temp = Temp,
          E = E)}
      
      return(sim_output)
    })
  }
  
  # Return data frame
  simulation_output <- compute.remaining()
  return(simulation_output)
}
functions <- c(functions, 'simulation')