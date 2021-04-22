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
explore.parm.space.3d <- function(n_pts, type, end_time = 2300,
                               stopping_points = c(2100,2200,2300),
                               lambda_init=0.9, omega_init=0.9,
                               debt_init = 0.3){
  # name file to save data:
  savefile <- sprintf("reduced_model/parms_res/npts_%g_type_%s_lambs_%g_omg_%g_d_%g_end%g.Rdata",
                      n_pts, type, lambda_init, omega_init, debt_init, end_time)
  if (!file.exists(savefile)) {
    # Set initial conditions 
    IC <- c(
      lambda  =  lambda_init, 
      omega   =  omega_init,
      debt    =  debt_init,
      pop     = 4.83
    )
    
    # Set up simulation
    Time <- c(
      start         =     2016, 
      end           =     end_time,
      step          =     0.05
    )
    
    Options <- list(
      transform_vars  = FALSE 
    )
  
    # set-up grid
    grid <- create.parm.grid (n_pts = n_pts, type = type)
    # list to save the results
    result <- vector("list",nrow(grid))
    
    # for loop
    for (i in seq(1,nrow(grid))){
      cat(i, 'out of', nrow(grid), '\n')
      Parms[['eta_p']] = grid[i,1]
      Parms[['markup']] = grid[i,2]
      Parms[['gamma']] = grid[i,3]
      
      # run simulation
      Sim <- try(simulation_red(time  = Time,
                            init_state  = IC,
                            parms   = Parms,
                            options = Options,
                            method = 'lsoda'))
      
      # record result for all stopping points and price params
      result[[i]] <- as_tibble(Sim) %>% 
        filter(year %in% stopping_points) %>% 
        mutate(eta = Parms[['eta_p']],
               markup = Parms[['markup']],
               gamma = Parms[['gamma']])
      }
    # Save results
    save(result, file=savefile)
    } else { 
    # read in the data 
    load(savefile)
      }
  return(result)
  }

#-----------------------------------------------------------------------------
# function to unpack list into a tibble
# there is probably a way to find this with purr but I didn't find it
##' @param res output of compute.basin.reduced
##' @return res reshaped into a single tibble
flatten.result <- function(res){
  out <- res[[1]]
  for (i in 2:length(res)){
    out <- rbind(out,res[[i]])
  }
  return(out)
}

# function to categorize result
##' @param res flattened output of compute.basin.reduced (single tibble with all years)
##' @return res with an additional variable 'outcome' containing categorized result
##' NOTE: if outcome is 'not assigned' this should be investigated
categorize.result <- function(res){
  result_categorized <- res %>%
    mutate(outcome = 'not assigned') %>% 
    mutate(outcome = ifelse((0.4 <= omega) & (1 > omega) & (0.4 <= lambda) &
                              (1 > lambda) & (2.7 >= debt_share), "good",outcome)) %>% 
    mutate(outcome = ifelse((omega > 1) | (lambda > 1), 'outside_bounds',outcome)) %>% 
    mutate(outcome = ifelse((omega < 0.4) | (lambda < 0.4) | (debt_share > 2.7), 
                            'bad', outcome)) 
  return(result_categorized)
}

# function generates interactive 3d scatter plot based on a result array
##' @param result tibble from categorize.result
interactive.scatter <- function(result){
  df <- result %>% mutate(outcome = as.factor(outcome))
  
  # Create the plot
  p <- plot_ly(
    df, x = ~eta, y = ~markup, z = ~gamma, size = 5,
    color = ~outcome, colors = c("#D95F02","#7570B3","#1B9E77")
  ) %>%
    add_markers() %>%
    layout(
      scene = list(xaxis = list(title = 'lambda'),
                   yaxis = list(title = 'omega'),
                   zaxis = list(title = 'd'))
    )
  return(p)
}


#------------------------------------------------------------------------------
# RESULTS
result1 <- explore.parm.space.3d(n_pts=20, lambda_init = 0.9, omega_init=0.9,
                                 debt_init=0.3, type="sobol", end_time = 2300,
                                 stopping_points = c(2100,2200,2300))

result1 <- result1 %>% flatten.result() %>% categorize.result()

# looking at good outcomes out of curiosity
result1 %>% filter(outcome == "good", year == 2300) %>%
  ggplot(aes(omega,lambda,colour=debt_share)) +
  geom_point() + theme_bw()

#interactive.scatter(result1)

result2 <- explore.parm.space.3d(n_pts=20, lambda_init = 0.5, omega_init=0.6,
                              debt_init=1, type="sobol", end_time = 2300,
                              stopping_points = c(2100,2200,2300))

result2 <- result2 %>% flatten.result() %>% categorize.result()

# looking at good outcomes out of curiosity
result2 %>% filter(outcome == "good", year == 2300) %>%
  ggplot(aes(omega,lambda,colour=debt_share)) +
  geom_point() + theme_bw()
