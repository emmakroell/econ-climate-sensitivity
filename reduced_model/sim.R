simulation_red <- function(time, init_state, options = Options, parms,
                           method = 'lsoda', montecarlo=F) {
  
  # timer
  start_time <- proc.time()
  # time sequence for ode
  time_seq <- seq(0, time[['end']] - time[['start']], by = time[['step']])
  
  # if using log(omega) and tan(pi(lambda-1/2))
  if (options$transform_vars == TRUE){
    init_vars <- c(
      lambda  = transform_lambda(init_state[[1]]),
      omega   = transform_omega(init_state[[2]]),
      debt    = init_state[[3]],
      pop     = init_state[[4]],
      counter = 0
    )
  }
  
  # otherwise, use regular variables
  else {
    init_vars <- c(init_state,
                   counter =0)
  }
  
  #============================================================================
  ## DIFFERENTIAL SYSTEM
  # Function for deSolve:
  Keen_w_prices <- function(t, state, parms) {
    with(as.list(c(parms,state)), {
      
      ### If using transformed variables
      if (options$transform_vars == TRUE){
        # Record transformed variables:
        lambda_trsnfd <- lambda 
        omega_trsnfd <- omega   
        # Obtain normal variables
        lambda <- return_lambda(lambda_trsnfd)
        omega <- return_omega(omega_trsnfd)
      }
      
      ### Auxiliary equations
      # Profit share
      profit_share <- 1 - omega - r * debt
      # inflation
      i  <- infl(omega, parms)
      # investment function
      kappa <- inv(profit_share, parms)
      # growth rate of the economy
      g  <-  kappa / nu - delta
      # Phillips curve function
      phill <- phill(lambda, parms)
      # dividends function
      div <- div(profit_share, parms)
      # population growth
      beta <- pop_dyn(pop,parms)
      
      ### Differential system 
      # If using transformed variables
      if (options$transform_vars == TRUE){
        # employment rate (d_tan_lambda)
        d_lambda <- (1 + lambda_trsnfd^2) * pi * lambda * (g - alpha - beta)
        # wage share (d_log_omega)
        d_omega  <- phill - alpha - (1 - gamma) * i
      }
      # Otherwise, use d_lambda and d_omega
      else {
        # employment rate
        d_lambda <- lambda * (g - alpha - beta)
        # wage share
        d_omega  <- omega * (phill - alpha - (1 - gamma) * i)
      }
      # debt share
      d_debt <- kappa - profit_share + div - debt * (i + g)
      # population (stock value)
      d_pop <- pop * beta
      # Counter 
      d_counter <- 1
      
      # ===========================================================================
      ## Return the list of gradients 
      return(list(c(
        d_lambda   =  d_lambda ,
        d_omega    =  d_omega, 
        d_debt     =  d_debt,
        d_pop      =  d_pop,
        d_counter  =  d_counter)))
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
  # REMAINING VARIABLES COMPUTATION'
  
  compute.remaining <- function(sim = simulation, params = c(parms,time)){
    with(as.list(c(sim,params)), {
      #### If using transformed variables
      if (options$transform_vars == TRUE){
        # record transformed variables:
        lambda_trsnfd <- lambda
        omega_trsnfd <- omega
        # recover original variables
        lambda <- return_lambda(lambda_trsnfd)
        omega <- return_omega(omega_trsnfd)
      }
      
      profit_share <- 1 - omega - r * debt   # profit share
      i  <- infl(omega,params)     # inflation
      div <- div(profit_share, parms)   # dividends
      kappa <- inv(profit_share, parms)  # investment
      
      if(montecarlo == T){
        to.extract <- which((start + simulation$time) %in% seq(2020,2120,1))
        sim_output <- data.frame(
          omega    = omega[to.extract],
          lambda   = lambda[to.extract],
          pi       = profit_share[to.extract],
          i        = i[to.extract],
          debt_share = debt[to.extract]
        )
        
      }else{
        # return results as a data frame
        sim_output <- data.frame(
          year     = start + time,
          counter  = counter,
          lambda   = lambda,
          omega    = omega,
          debt_share = debt,
          pop = pop,
          profit_share = profit_share,
          kappa    = kappa,
          i        = i,
          div      = div)}
      
      return(sim_output)
    })
  }
  # Return data frame
  simulation_output <- compute.remaining()
  return(simulation_output)
}

functions <- c(functions, 'simulation_red')