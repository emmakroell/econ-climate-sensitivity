library(tidyverse)
# compare pricing params to data from M. Grasselli
params <- read.csv('parameters/Grasselli_parameters.csv')
params <- as.data.frame(params)
# generalized gamma dist code from B. Bolker
source("full_model/gen_gamm_dist.R") 

unc_parms <- c(
  eta_mean  = 0.4,
  eta_sd    = 0.12
)

# eta
x <- seq(0, 1, length=1000)
y <- dnorm(x, mean=unc_parms[['eta_mean']], sd=unc_parms[['eta_sd']])
xy <- as.data.frame(cbind(x,y))
ggplot() + geom_histogram(data=params, mapping=aes(x=eta), binwidth = 0.05) +
  geom_line(data = xy, mapping=aes(x=x,y=y),colour='red')+ 
  labs(title=expression('Distribution of'~ eta))


# gamma
pars <- get_ggamma(q=c(0.05,0.4,0.95),val=c(0.1,0.09,1),control=list(maxit=2000),
                   shift=0)
x <- seq(0, 0.9999, length=1000)
y <- with(as.list(pars),dggamma(1-x,s,m,f))
# print s, m, f for gamma distribution:
with(as.list(pars),cat('s=',s,'m=',m,'f=',f))

1 - with(as.list(pars),qggamma(0.5,s,m,f))   # mean
xy <- as.data.frame(cbind(x,y))
ggplot() + geom_histogram(data=params, mapping=aes(x=gamma), binwidth = 0.05) +
  geom_line(data = xy, mapping=aes(x=x,y=y),colour='red')+ 
  labs(title=expression('Distribution of'~gamma))


#xi
x <- seq(1.0001, 3, length=1000)
pars <- get_ggamma(q=c(0.05,0.4,0.95),val=c(0.2,0.55,1.5),control=list(maxit=2000),
                   shift=0)
y <- with(as.list(pars),dggamma(x-1,s,m,f))
# print s, m, f for markup distribution
with(as.list(pars),cat('s=',s,'m=',m,'f=',f))

1 + with(as.list(pars),qggamma(0.5,s,m,f))  ### mean
xy <- as.data.frame(cbind(x,y))
ggplot() +
  geom_histogram(data=params, mapping=aes(x=xi), binwidth = 0.05) +
  geom_line(data = xy, mapping=aes(x=x,y=y),colour='red') + 
  labs(title=expression('Distribution of'~xi))

 


