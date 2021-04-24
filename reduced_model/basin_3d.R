# FILE FOR FIGURE 2
# Packages
library('here')
library('deSolve')
# library('RColorBrewer')
library('qrng')
library('tidyverse')
# library('plotly')
# library('rgl')
# library('gMOIP')

## SET 'ncores' here for parallel computation
## (left equal to 1 by default, unless edited, for safety)
ncores <- getOption("mc.cores", 1)
## ncores <- 10

# Source model code
source('reduced_model/pars.R')    # load parameters
source('reduced_model/funcs.R')   # load functions
source('reduced_model/sim.R')     # load simulation file 

#------------------------------------------------------------------------------
# function creates grid of initial conditions
create.ic.grid <- function(n_pts=20, type = 'grid'){
  if (type == 'grid'){
    lambda.grid  <- seq(0.2, 0.99, length.out = n_pts)   
    omega.grid   <- seq(0.2, 0.99, length.out = n_pts)  
    debt.grid    <- seq(0.1, 2.7, length.out = n_pts)   
    res <- expand.grid(x = lambda.grid, y = omega.grid, z = debt.grid)
    
  } else if (type == 'sobol'){
    var_ranges <- matrix(c(0.2,0.99,0.2,0.99,0.1,2.7),ncol=2, byrow=TRUE)
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

#------------------------------------------------------------------------------
# main function for basin of attraction
##' @param n_pts number of points to run per dimension. Total runs is n_pts^3
##' @param end_time when to end the simulation
##' @param stopping_points vector of times to report results (should include end_time)
##' @param eta value of eta (fixed in all runs)
##' @param markup value of markup (fixed in all runs)
##' @param gamma value of gamma (fixed in all runs)
compute.basin.reduced <- function(n_pts, eta=0.192, markup = 1.18,
                             gamma=0.9, end_time = 2500,
                             stopping_points = c(2100,2300,2500),
                             type='sobol'){
  # name file to save data:
  savefile <- sprintf("reduced_model/basin_res/npts_%g_type_%s_eta_%g_mark_%g_gam_%g_end_%g.Rdata",
                      n_pts, type, eta, markup, gamma, end_time)
  if (!file.exists(savefile)) {
    # set-up grid
    grid <- create.ic.grid(n_pts = n_pts, type = type)

    # list to save the results
    result <- vector("list",nrow(grid))

    # set parameters
    Parms[['eta_p']] = eta
    Parms[['markup']] = markup
    Parms[['gamma']] = gamma

   loop_fun <- function(i) {
                # Set initial conditions 
      IC <- c(
        lambda = grid[i,1], 
        omega  = grid[i,2],
        debt   = grid[i,3],
        pop    = 4.83
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
      
      # Run simulation
      Sim <- try(simulation_red(time = Time,
                               init_state = IC,
                               parms      = Parms,
                               options    = Options,
                               method     = 'lsoda'))
      
      # record result for all stopping points and original IC
      res <- as_tibble(Sim) %>% 
        filter(year %in% stopping_points) %>% 
        mutate(lambda.ic = IC[['lambda']],
               omega.ic = IC[['omega']],
               debt.ic = IC[['debt']])
             return(res)
      } ## loop_fun

    # for loop
    if (ncores==1) {
        for (i in seq(1,nrow(grid))){
            #cat(i, 'out of', nrow(grid), '\n')
            results[[i]] <- loop_fun(i)
        }
    } else {
        papply <- if (!require("pbmcapply")) parallel::mclapply else pbmclapply

        ## pbmclapply for progress bar (may only work in interactive mode?)
        results <- papply(seq(nrow(grid)), loop_fun,
                                      mc.cores = ncores)
    }   # Save results
    save(result, file=savefile)
  } else { 
    # read in the data 
    load(savefile)
  }
  return(result)
}

#------------------------------------------------------------------------------
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
    df, x = ~lambda.ic, y = ~omega.ic, z = ~debt.ic, size = 5,
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
result_3d1 <- compute.basin.reduced(n_pts=20, markup=1.18, end_time = 2500,
                                    stopping_points = c(2100,2300,2500))

# result1 <- result_3d1 %>% flatten.result() %>% categorize.result()
# 
# # Find equilibrium
# result1 %>% filter(outcome == "good", year == 2500) %>%
#   ggplot(aes(omega,lambda,colour=debt_share)) +
#   geom_point() + theme_bw()
# # appears to be approaching a single equilibrium
# 
# # compute equilibrium using median
# equilibrium1 <- result1 %>% filter(outcome == "good", year == 2500) %>% 
#   dplyr::summarise(lambda.eq = median(lambda),
#                     omega.eq = median(omega),
#                     debt.eq = median(debt_share))

# plot output
#interactive.scatter(result1)


#----
result_3d2 <- compute.basin.reduced(n_pts=20, markup=1.3, end_time = 2500,
                                    stopping_points = c(2100,2300,2500))

# result2 <- result_3d2 %>% flatten.result() %>% categorize.result()
# 
# # Find equilibrium
# result2 %>% filter(outcome == "good", year == 2500) %>%
#   ggplot(aes(omega,lambda,colour=debt_share)) +
#   geom_point() + theme_bw()
# # appears to be approaching a single equilibrium
# 
# # compute equilibrium using median
# equilibrium2 <- result2 %>% filter(outcome == "good", year == 2500) %>% 
#   dplyr::summarise(lambda.eq = median(lambda),
#                     omega.eq = median(omega),
#                     debt.eq = median(debt_share))

#----
result_3d3 <- compute.basin.reduced(n_pts=20, markup=1.875, end_time = 2500,
                                    stopping_points = c(2100,2300,2500))

# result3 <- result_3d3 %>% flatten.result() %>% categorize.result()
# 
# # Find equilibrium
# result3 %>% filter(outcome == "good", year == 2500) %>%
#   ggplot(aes(omega,lambda,colour=debt_share)) +
#   geom_point() + theme_bw()
# # appears to be approaching a single equilibrium
# 
# # compute equilibrium using median
# equilibrium3 <- result3 %>% filter(outcome == "good", year == 2500) %>% 
#   dplyr::summarise(lambda.eq = median(lambda),
#                    omega.eq = median(omega),
#                    debt.eq = median(debt_share))


