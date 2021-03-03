#MONTE CARLO SIMULATION FUNCTION

## Draw a simulation
draw_simulation <- function(time, init_state, parms, options, montecarlo, iter_num) {
  
  # Drawing of the parameters
  parms['eta_p']  <-  draw_eta()
  parms['markup'] <-  draw_markup()
  parms['gamma']  <-  draw_gamma()
  parms['alpha']          <-  draw_alpha()
  parms['S']              <-  draw_ECS()
  parms['CO2_UP_preind']  <-  draw_C_UP_preind()
  
  # Run of a simulation
  Draw_simu <- try(simulation(time       = time,
                              init_state = init_state, 
                              parms       = parms, 
                              options    = options,
                              montecarlo = T))
  
  
  # Return
  return(list(draw = c(parms['eta_p'], parms['markup'], parms['gamma']),
              parms['alpha'], parms['S'], parms['CO2_UP_preind'],
              sim = Draw_simu))
}

## Process several iterations of drawings
data_simulations <- function(iter             = 10,
                             time             = Time,
                             init_state       = IC, 
                             parms            = Parms, 
                             options          = Options,
                             montecarlo       = T) {
  
  # Begin of calculation
  start_time  <- proc.time()
  
  # Start parallel computing
  mc.cores <- getOption("mc.cores", max(1, detectCores()-1))
  cl <- makeCluster(mc.cores)
  registerDoParallel(cl)
  
  # Computation
  data <- foreach(k=1:iter, .packages=c('deSolve','rmutil'),
                  .export = functions, .combine='cbind') %dopar% {
                    drawing <- draw_simulation(time        = time,
                                               init_state  = init_state, 
                                               parms       = parms, 
                                               options    = options,
                                               montecarlo = montecarlo,
                                               iter_num = k)
                    
                    return(list(c(draw = drawing$draw, sim = drawing$sim)))
                  }
  
  # Start parallel computing
  stopCluster(cl)
  
  # End of calculation
  end_time <- proc.time() - start_time
  print(end_time)

  # Returning
  return(data)
}


functions <- c(functions,
               'draw_eta',
               'draw_markup',
               'draw_gamma',
               'draw_alpha',
               'draw_ECS',
               'draw_C_UP_preind',
               'draw_simulation',
               'get_ggamma')


