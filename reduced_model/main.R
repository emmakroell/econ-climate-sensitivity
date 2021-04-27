# MAIN FILE
# runs the reduced form model (no climate)

# Packages
library(deSolve)
library(tidyverse)

# Set-up 
source('reduced_model/pars.R')   # load parameters
source('reduced_model/funcs.R')  # load functions
source('reduced_model/sim.R')    # load simulation file
#==============================================================================
### Set initial conditions 
IC <- c(
  lambda  =  0.675, 
  omega   =  0.578,
  debt    =  1.53,
  pop     =  4.83
)

# Set up simulation
Time <- c(
  start         =     2016, 
  end           =     2300,
  step          =     0.05
)

Options <- list(
  transform_vars  = FALSE  # T / F
)

#================================================================================
## Run a scenario
Sim <- simulation_red(time       = Time,
                      init_state = IC,
                      parms      = Parms,
                      method     = 'lsoda')

#==============================================================================
# Make a plot:
as_tibble(Sim) %>% 
  select(year, lambda, omega, debt_share, profit_share, i, pop) %>% 
  pivot_longer(cols = lambda:pop, names_to = 'variable', values_to = 'value') %>% 
  mutate(variable = factor(variable,levels = c("lambda", "omega", "profit_share",
                                               "debt_share", "i", "pop"))) %>% 
  ggplot(aes(year,value,colour = variable)) + geom_line(size=1.1) +
  facet_wrap(~variable, scale = "free") + theme_bw() +
  theme(legend.position = "none")

