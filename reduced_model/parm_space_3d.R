library('here')
library('deSolve')
library('RColorBrewer')
library('qrng')
library('scatterplot3d')
library('tidyverse')
library('plotly')
library('rgl')
library('gMOIP')

source('reduced_model/pars.R')   # load parameters
source('reduced_model/sim.R')    # load simulation file
source('reduced_model/funcs.R')  # load functions
# =============================================================================
# Set up simulation
Time <- c(
  start         =     2020, 
  end           =     2100,
  step          =     0.05
)

Options <- list(
  transform_vars  = FALSE 
)
#================================================================================
create.parm.grid <- function(n_pts=20, type = 'sobol'){
  if (type == 'grid'){
    eta.grid    <- seq(0, 4, length.out = n_pts)   # eta
    markup.grid <- seq(1, 3, length.out = n_pts)   # markup
    gamma.grid  <- seq(0, 1, length.out = n_pts)   # gamma
    res <- expand.grid(x = eta.grid, y = markup.grid, z = gamma.grid)
    
  } else if (type == 'uniform'){
    #res <- replicate(3, runif(n_pts^3))
    eta.grid    <- 4 * runif(n_pts)   # eta
    markup.grid <- 1 + 2 * runif(n_pts)   # markup
    gamma.grid  <- runif(n_pts)   # gamma
    res <- expand.grid(x = eta.grid, y = markup.grid, z = gamma.grid)
    
  } else if (type == 'sobol'){
    var_ranges <- matrix(c(0,1,1,2.5,0,1),ncol=2, byrow=TRUE)
    colnames(var_ranges) <- c("min","max")
    res <- qrng::sobol(n_pts^3,d=3) 
    for (i in 1:ncol(res)) {
      minval <- var_ranges[i,1]
      maxval <- var_ranges[i,2]
      res[,i] <- minval + (maxval-minval)*res[,i]
    }
  }
  return(res)
}

#================================================================================
explore.parm.space <- function(n_pts, type, trnsfm = FALSE, 
                               plot.res = FALSE, lambda_init=0.9, 
                               omega_init=0.9, debt_init = 0.3){
  # name file to save data:
  savefile <- sprintf("reduced_model/parms_res/npts_%g_type_%s_trnsfm_%s_lambs_%g_omg_%g_d_%g.Rdata",
                      n_pts, type, trnsfm, lambda_init, omega_init, debt_init)
  if (!file.exists(savefile)) {
    # Set initial conditions 
    IC <- c(
      lambda  =  lambda_init, 
      omega   =  omega_init,
      debt    =  debt_init,
      pop     = 4.83
    )
  
    # set-up grid
    grid <- create.parm.grid (n_pts = n_pts, type = type)
    outcome <- rep(NA, nrow(grid))
    cputime.start <- proc.time() # start timer
    
    # for loop
    for (i in seq(1,nrow(grid))){
      cat(i, 'out of', nrow(grid), '\n')
      Parms[['eta_p']] = grid[i,1]
      Parms[['markup']] = grid[i,2]
      Parms[['gamma']] = grid[i,3]
      
      # Transform variables:
      Options <- list(
        transform_vars  = trnsfm
      )
      
      # run simulation
      Sim <- try(simulation_red(time  = Time,
                            init_state  = IC,
                            parms   = Parms,
                            method = 'lsoda'))
      
      # if the simulation encounters an error, record this and keep going:
      if (inherits(Sim,"try-error")){
        outcome[i] = 'error'
      } else if (is.na(tail(Sim$lambda, 1)) == TRUE){
        outcome[i] = 'error'
      }
      else{
        # categorize as good / bad / error
        if ((0.4 <= tail(Sim$omega,1)) & (1 > tail(Sim$omega,1)) &
            (0.4 <= tail(Sim$lambda,1)) & (1 > tail(Sim$lambda,1)) &
            (2.7 >= tail(Sim$debt_share,1))){
          outcome[i] = 'good'
        } else if ((tail(Sim$omega,1) > 1) | (tail(Sim$lambda,1) > 1)){
          outcome[i] = 'outside_bounds'
        } else{
          outcome[i] = 'bad'
        }}}
    cputime.end <- proc.time() # end timer
    cputime <- cputime.end - cputime.start
    # Result saved in data frame with parameters
    result <- data.frame(grid,outcome)
    colnames(result) <- c("eta", "markup", "gamma", "outcome")
    # Save results
    save(result, cputime, file=savefile)
  } else { 
    # read in the data 
    load(savefile)
  }
  
  print(result)
  if (plot.res == TRUE){
    par(mfrow = c(1,1), las=1, xpd=T)
    colors <- c("#D95F02", "#1B9E77","#7570B3",'red')
    colors <- colors[as.numeric(result$outcome)]
    scatterplot3d(result[,1:3], pch = 16, color=colors,
                  main="Parameter grid search", xlab = '',
                  ylab = '', zlab = '')
    legend(x=0, y=-1, legend = levels(result$outcome),
           col = c("#D95F02", "#1B9E77","#7570B3",'red'), pch = 16,
           inset = -0.25, xpd = TRUE, horiz = TRUE)
    text(x = 7, y = 0.2, expression(xi), srt = 45)
    text(x = -0.8, y = 2.7, expression(gamma), srt = 0)
    text(x = 2.5, y = -1, expression(eta), srt = 0)
  }
  
  return(result)
}

#------------------------------------------------------------------------------
# Interactive plot
# function generates interactive 3d scatter plot based on a result array
interactive.scatter <- function(result){
  df <- result %>%
    rownames_to_column() %>%
    as_tibble() %>%
    mutate(outcome = as.factor(outcome))
  
  # Create the plot
  p <- plot_ly(
    df, x = ~eta, y = ~markup, z = ~gamma, size = 5,
    color = ~outcome, colors = c("#D95F02","#7570B3","#1B9E77")
  ) %>%
    add_markers() %>%
    layout(
      scene = list(xaxis = list(title = 'eta'),
                   yaxis = list(title = 'markup'),
                   zaxis = list(title = 'gamma'))
    )
  return(p)
}

#------------------------------------------------------------------------------
# RESULTS
result1 <- explore.parm.space(n_pts=20, lambda_init = 0.9, omega_init=0.9,
                              debt_init=0.3, type="sobol")
result1$outcome[result1$outcome == 'outside_bounds']='error'
interactive.scatter(result1)

result2 <- explore.parm.space(n_pts=20, lambda_init = 0.5, omega_init=0.6,
                              debt_init=1, type="sobol")
result2$outcome[result2$outcome == 'outside_bounds']='error'
interactive.scatter(result2)

result3 <- explore.parm.space(n_pts=20, lambda_init = 0.675, omega_init=0.578,
                              debt_init=1.53, type="sobol")
result3$outcome[result3$outcome == 'outside_bounds']='error'
interactive.scatter(result3)
