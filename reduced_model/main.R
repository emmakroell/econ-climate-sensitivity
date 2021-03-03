# MAIN FILE
# runs the reduced form model (no climate)

# Packages
library("deSolve")
library("RColorBrewer")
# Colour palettes for graphs
greens <- brewer.pal(n = 6, name = "Greens")
colourful <- brewer.pal(n = 3, name = "Set1")
colour2 <- brewer.pal(n = 3, name = "Set2")

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
  pop     = 4.83
)

# Set up simulation
Time <- c(
  start         =     2020, 
  end           =     2100,
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
# PLOTS
par(mfrow = c(1,1), las=1, xpd=T)
plot(x = Sim$year, y = Sim$lambda, type = 'l',col = colourful[1],
     xlab = '', ylab = '', ylim = c(-0.1,1),
     main=expression(paste("Key economic variables")))
lines(x = Sim$year, y = Sim$omega, type = 'l',col = colourful[2])
lines(x = Sim$year, y = Sim$profit_share, type = 'l',col = colourful[3])
legend(x=2088, y=0.15, legend = c(expression(paste(lambda)),
                                  expression(paste(omega)),
                                  expression(paste(pi))),
       col=colourful, lty=1:1,box.lty=0, cex=0.75)

plot(x = Sim$year, y = Sim$debt, type = 'l',col = 'purple', xlab = '', 
     ylab = '', main=expression(paste("Debt share (d)")))

# plot(x = Sim$year, y = Sim$i, type = 'l',
plot(x = Sim$year, y = Sim$i, type = 'l',col = 'blue', xlab = '',
     ylab = '', main=expression(paste("Inflation (i)")))
