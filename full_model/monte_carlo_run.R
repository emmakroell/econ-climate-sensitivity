# MONTE CARLO CODE
# Packages 
library(doParallel)
library(tidyverse)
library(gridExtra)


set.seed(720)
source('full_model/pars.R')     # load parameters
source('full_model/funcs.R')    # load functions
source('full_model/sim.R')      # load simulation file
source("full_model/param_dist.R") # load parameter distributions
source("full_model/monte_carlo_sim.R") # load monte carlo functions

# Number of iterations
iter <- 100

# Set up simulations
lambda_init <-0.675
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

Time <- c(
  start         =     2016, 
  end           =     2100,
  step          =     0.05
)

#==============================================================================
# No feedback scenario
dam <- 'None'

Options <- list(
  invest = 'lin',  # exp / arctan / lin
  damage_scenario = dam,
  p_carb_scheme = "None",
  subsidy = 0   # fraction of abatement costs subsidized by government
)

## MONTE CARLO SIMULATION
# name file to save data:
savefile <- sprintf("full_model/MC_res/niter_%g_dam_%s.Rdata",
                    iter, dam)
if (!file.exists(savefile)){
  mc_nofeedback <- data_simulations(iter)
  bad_runs <- vapply(mc_nofeedback,is.character,FUN.VALUE=logical(1))
  mc_nofeedback <- mc_nofeedback[!bad_runs]
  save(mc_nofeedback, file=savefile)
} else{
  load(savefile)
}

# Reformulate data using tidyverse
mc_nofeedback_tibble <- map(mc_nofeedback, ~tibble(
  year=2016:2100, 
  lambda=.$sim.lambda,
  debt_share = .$sim.debt_share,
  growth_Y =(.$sim.Y-lag(.$sim.Y))/lag(.$sim.Y),
  infl =.$sim.i,
  emissions = .$sim.E,
  temp = .$sim.Temp)) 

## combine elements into a single data frame
mc_nofeedback_long <- (bind_rows(mc_nofeedback_tibble,.id="rep")
        ## collapse statevars to long format
        %>% gather(statevar,value,-c(rep,year))
)

ggplot(mc_nofeedback_long, aes(year,value,group=rep))+
  geom_line(alpha=0.5)+
  facet_wrap(~statevar,scale="free")

## summarize and draw lower/upper quantiles
mc_nofeedback_summarized <- mc_nofeedback_long %>%
  group_by(year,statevar) %>%
  summarise(median=median(value,na.rm=T),
            lwr=quantile(value,0.05,na.rm=T),
            upr=quantile(value,0.95,na.rm=T)) %>% 
  mutate(statevar = factor(statevar,levels = c('lambda', 'growth_Y', 'debt_share',
                                               'infl', 'emissions', 'temp'))) 


ggplot(mc_nofeedback_summarized,aes(year,y=median,ymin=lwr,ymax=upr))+
  geom_line(colour='blue') + geom_ribbon(colour=NA,alpha=0.3,fill='blue') +
  facet_wrap(~statevar,labeller = as_labeller(c(lambda = "Employment rate",
                                                debt_share = "Debt share",
                                                growth_Y = "Output growth",
                                                infl = "Inflation",
                                                emissions = "Emissions",
                                                temp = "Temperature")),scale="free")

#==============================================================================
# Model with damages and policy
dam <- 'Nordhaus'
Options <- list(
  invest = 'lin',  # exp / arctan / lin
  damage_scenario = dam,
  p_carb_scheme = "Stern-Stiglitz",
  subsidy = 0.5   # fraction of abatement costs subsidized by government
)

## MONTE CARLO SIMULATION
# name file to save data:
savefile <- sprintf("full_model/MC_res/niter_%g_dam_%s.Rdata",
                    iter, dam)
if (!file.exists(savefile)){
  mc_dam <- data_simulations(iter)
  bad_runs <- vapply(mc_dam,is.character,FUN.VALUE=logical(1))
  mc_dam <- mc_dam[!bad_runs]
  save(mc_dam, file=savefile)
} else{
  load(savefile)
}

# Reformulate data using tidyverse
mc_dam_tibble <- map(mc_dam, ~tibble(
  year = 2016:2100, 
  lambda = .$sim.lambda,
  debt_share = .$sim.debt_share,
  growth_Y =(.$sim.Y-lag(.$sim.Y))/lag(.$sim.Y),
  infl =.$sim.i,
  emissions = .$sim.E,
  temp = .$sim.Temp)) 

## combine elements into a single long tibble
mc_dam_long <- (bind_rows(mc_dam_tibble,.id="rep")
        ## collapse statevars to long format
        %>% gather(statevar,value,-c(rep,year))
)

ggplot(mc_dam_long, aes(year,value,group=rep))+
  geom_line(alpha=0.5)+
  facet_wrap(~statevar,scale="free")

## summarize and draw lower/upper quantiles
mc_dam_summarized <- mc_dam_long %>% 
  group_by(year,statevar) %>%
  summarise(median=median(value,na.rm=T),
            lwr=quantile(value,0.05,na.rm=T),
            upr=quantile(value,0.95,na.rm=T)) %>% 
  mutate(statevar = factor(statevar,levels = c('lambda', 'growth_Y', 'debt_share',
                                               'infl', 'emissions', 'temp'))) 

ggplot(mc_dam_summarized,aes(year,y=median,ymin=lwr,ymax=upr))+
  geom_line(colour='red') + geom_ribbon(colour=NA,alpha=0.3,fill='red') +
  facet_wrap(~statevar,labeller = as_labeller(c(lambda = "Employment rate",
                                                debt_share = "Debt share",
                                                growth_Y = "Output growth",
                                                infl = "Inflation",
                                                emissions = "Emissions",
                                                temp = "Temperature")),scale="free")

#==============================================================================
# COMBINED PLOT
ggplot(data=mc_dam_summarized)+geom_line(aes(year,y=median),colour='red') +
  geom_ribbon(colour=NA,aes(year,y=median,ymin=lwr,ymax=upr),alpha=0.3,fill='red') +
  geom_line(data=mc_nofeedback_summarized,aes(year,y=median),colour='blue') +
  geom_ribbon(data=mc_nofeedback_summarized,colour=NA,aes(year,y=median,ymin=lwr,ymax=upr),alpha=0.3,fill='blue')+
  facet_wrap(~statevar,labeller = as_labeller(c(lambda = "Employment rate",
                                                debt_share = "Debt share",
                                                growth_Y = "Output growth",
                                                infl = "Inflation",
                                                emissions = "Emissions",
                                                temp = "Temperature")),scale="free") +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  

#==============================================================================
# PROBABILITIES
# compute probability temp < 2, debt < 2.7, and both
# temp under 2
(temp_nf <- mc_nofeedback_long %>% 
  filter(statevar == 'temp', year == 2100) %>% 
  mutate(less_2 = ifelse(value < 2, TRUE, FALSE)) %>% 
  group_by(less_2) %>% 
  summarise(count = n()))

(temp_dam <- mc_dam_long %>% 
  filter(statevar == 'temp', year == 2100) %>% 
  mutate(less_2 = ifelse(value < 2, TRUE, FALSE)) %>% 
  group_by(less_2) %>% 
  summarise(count = n()))

col1 <- c(temp_nf$count[2]/sum(temp_nf$count),
          temp_dam$count[2]/sum(temp_dam$count))

# temp under 3 - just to see
temp3_nf <- mc_nofeedback_long %>% 
  filter(statevar == 'temp', year == 2100) %>% 
  mutate(less_3 = ifelse(value < 3, TRUE, FALSE)) %>% 
  group_by(less_3) %>% 
  summarise(count = n())

temp3_dam <- mc_dam_long %>% 
  filter(statevar == 'temp', year == 2100) %>% 
  mutate(less_3 = ifelse(value < 3, TRUE, FALSE)) %>% 
  group_by(less_3) %>% 
  summarise(count = n())

(probs_3d <- c(temp3_nf$count[2]/sum(temp3_nf$count),
          temp3_dam$count[2]/sum(temp3_dam$count)))

#-------------------------------------------------------------------
# debt share under 2.7
(debt_share_nf <- mc_nofeedback_long %>% 
   filter(statevar == 'debt_share', year == 2100) %>% 
   mutate(less_2.7 = ifelse(value < 2.7, TRUE, FALSE)) %>% 
   group_by(less_2.7) %>% 
   summarise(count = n()))

(debt_share_dam <- mc_dam_long %>% 
    filter(statevar == 'debt_share', year == 2100) %>% 
    mutate(less_2.7 = ifelse(value < 2.7, TRUE, FALSE)) %>% 
    group_by(less_2.7) %>% 
    summarise(count = n()))

col2 <- c(debt_share_nf$count[2]/sum(debt_share_nf$count),
          debt_share_dam$count[2]/sum(debt_share_dam$count))

#----------------------------------------------------------------
# temp under 2 and debt share under 2.7
(both_nf <- mc_nofeedback_long %>% 
  filter(statevar %in% c('temp','debt_share'), year == 2100) %>% 
  pivot_wider(names_from = statevar, values_from=value) %>% 
  mutate(less_both = ifelse(temp < 2  & debt_share <= 2.7, TRUE, FALSE)) %>% 
  group_by(less_both) %>% 
  summarise(count = n()))

(both_dam <- mc_dam_long %>% 
  filter(statevar %in% c('temp','debt_share'), year == 2100) %>% 
  pivot_wider(names_from = statevar, values_from=value) %>% 
  mutate(less_both = ifelse(temp < 2  & debt_share <= 2.7, TRUE, FALSE)) %>% 
  group_by(less_both) %>% 
  summarise(count = n()))

col3 <- c(both_nf$count[2]/sum(both_nf$count),
          both_dam$count[2]/sum(both_dam$count))

# create table 2
table <- cbind(col1,col2, col3)
colnames(table) <- c("Temp < 2", "Debt share < 2.7", "Both")
table

#----------------------------------------------------------------
# MEDIANS
mc_nofeedback_long %>% 
  filter(statevar == 'temp', year == 2100) %>%
  summarize(temp = median(value))

mc_dam_long %>% 
  filter(statevar == 'temp', year == 2100) %>%
  summarize(temp = median(value))

#-------------------------------------------------------------------------
# Make a custom facet plot with a ylim on debt_share only
# this is probably a sub-optimal way to do this
for (var in c('lambda', 'growth_Y', 'debt_share', 'infl', 'emissions', 'temp')){
  this_plot <- ggplot(data=(mc_dam_summarized %>% filter(statevar==var)))+
    geom_line(aes(year,y=median),colour='red') +
    geom_ribbon(colour=NA,aes(year,y=median,ymin=lwr,ymax=upr),alpha=0.3,fill='red') +
    geom_line(data=(mc_nofeedback_summarized  %>% filter(statevar==var)),
              aes(year,y=median),colour='blue') +
    geom_ribbon(data=(mc_nofeedback_summarized %>% filter(statevar==var)),
                colour=NA,aes(year,y=median,ymin=lwr,ymax=upr),alpha=0.3,fill='blue') +
    ggtitle(var) + xlab('') + ylab('')
  assign(paste0('plot_',var), this_plot)
}

# add titles and the limits on debt_share
plot_lambda <- plot_lambda + ggtitle("Employment rate")
plot_growth_Y <- plot_growth_Y + ggtitle("Output growth")
plot_debt_share <- plot_debt_share + 
  coord_cartesian(xlim = c(2016, 2100), ylim = c(0,5)) +
  ggtitle("Debt share")
plot_infl <- plot_infl + ggtitle("Inflation")
plot_emissions <- plot_emissions + ggtitle("Emissions")
plot_temp <- plot_temp + ggtitle("Temperature")

# make into facet like plot using gridExtra package
grid.arrange(plot_lambda,plot_growth_Y,plot_debt_share,
             plot_infl, plot_emissions, plot_temp, nrow=2)

#ggsave("figure4.pdf",height=7,width=5)
