# MAIN FILE
# Packages
library("deSolve")
library("RColorBrewer")

# Colour palettes for graphs
greens <- brewer.pal(n = 6, name = "Greens")
colourful <- brewer.pal(n = 3, name = "Set1")
colour2 <- brewer.pal(n = 3, name = "Set2")

# Set-up
source('full_model/pars.R')     # load parameters
source('full_model/funcs.R')    # load functions
source('full_model/sim.R')      # load simulation file
#==============================================================================
# ### Set initial conditions
lambda_init <- 0.675
omega_init <-  0.578
debt_share_init <-1.53
Y <- 59.74
nu <- 2.7


IC <- c(
  K       =  Y * nu,
  w       =  omega_init * Y / ( lambda_init * 4.83),
  p       =  1,
  debt    =  debt_share_init * Y,
  a       =  Y / ( lambda_init * 4.83),
  pop     =  4.83,
  g_sigma =  -0.0105,  # Growth rate of the emissions intensity of the economy
  n       =  0.03,     # Emissions reduction rate
  E_ind   =  35.85,    # Industrial CO2-e emissions, in Gt CO2-e
  CO2_AT  =  851,      # CO2-e concentration in the atmosphere layer, in Gt C
  CO2_UP  =  460,      # CO2-e concentration in the biosphere and upper ocean layer, in Gt C
  CO2_LO  =  1740,     # CO2-e concentration in the deep ocean layer, in Gt C
  E_land  =  2.6,      # Exogenous land use change CO2-e emissions, in Gt CO2-e
  p_BS    =  547.22,   # Backstop technology price level
  Temp    =  0.85,     # Temperature in the atmosphere, biosphere, and upper ocean layer, in degrees Celsius
  Temp_LO =  0.0068    # Temperature in the deep ocean layer, in degrees Celsius
)


#=================================================================================
# Set up simulation
Time <- c(
  start         =     2016, 
  end           =     2100,
  step          =     0.05
)

Options <- list(
  invest = 'lin',  # exp / arctan / lin
  damage_scenario = 'Nordhaus',
  p_carb_scheme = "Stern-Stiglitz",
  subsidy = 0.5   # fraction of abatement costs subsidized by government
)


#================================================================================
Sim <- simulation(time       = Time,
                  init_state = IC,
                  parms      = Parms,
                  options    = Options,
                  method     = 'lsoda')

as_tibble(Sim) %>% mutate(Y = K / 2.7) %>% 
  select(year, Y, debt_share, lambda, p_Car, E, Temp) %>% 
  pivot_longer(cols = Y:Temp, names_to = 'variable', values_to = 'value') %>% 
  mutate(variable = factor(variable, levels = c("Y", "debt_share", "lambda", "p_Car", "E", "Temp"))) %>% 
  ggplot(aes(year,value,colour = variable)) + geom_line(size=1.1) +
  facet_wrap(~variable, scale = "free") + theme_bw()
ggsave("plot.png",height = 5, width =7)

