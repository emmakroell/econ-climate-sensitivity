# Logistic regression + PRCC
library('sensitivity')
library('here')

# Set-up 
source('reduced_model/main.R')   # load main
source('reduced_model/pars.R')   # load parameters
source('reduced_model/funcs.R')  # load functions
source('reduced_model/sim.R')    # load simulation file
source('reduced_model/param_dist_noclimate.R') 
parms <- Parms
set.seed(720)

# Set up simulation
Time <- c(
  start         =     2016, 
  end           =     2300,
  step          =     0.05
)

Options <- list(
  transform_vars  = FALSE 
)

#------------------------------------------------------------------------------
# function to compute PRCC
categorize.result <- function(n,lambda_init,omega_init,debt_share_init){
  # n: number of model runs
  # lambda_init: initial employment rate
  # omega_init: initial wage share
  # debt_share_init: initial debt share
  filename <- sprintf("reduced_model/prcc_data/npts_%g_lambs_%g_omg_%g_d_%g.Rdata",
                      n, lambda_init, omega_init, debt_share_init)
  
  if (!file.exists(filename)){
    # draw pricing parameters from distributions
    eta    <- replicate(n,draw_eta())
    markup <- replicate(n,draw_markup())
    gamma  <- replicate(n,draw_gamma())
    alpha  <- replicate(n,draw_alpha())
    
    # Set initial conditions 
    IC <- c(
      lambda  =  lambda_init, 
      omega   =  omega_init,
      debt    =  debt_share_init,
      pop     = 4.83
    )

    # arrays to store outcomes
    outcome <- rep(NA,n)
    good.outcome <- rep(NA,n)
    
    good.eta <- rep(NA,n)
    good.markup<- rep(NA,n)
    good.gamma <- rep(NA,n)
    good.alpha <- rep(NA,n)
    
    for (i in 1:n){
      Parms[['eta_p']] = eta[i]
      Parms[['markup']] = markup[i]
      Parms[['gamma']] = gamma[i]
      Parms[['alpha']] = alpha[i]

      ## Run a scenario
      Sim <- try(simulation_red(time       = Time,
                                init_state = IC,
                                parms      = Parms,
                                method     = 'lsoda'))
      if (class(Sim) == "try-error"){
        eta[i]     <- NA
        markup[i]  <- NA
        gamma[i]   <- NA
        alpha[i]   <- NA
      } else{
        if (is.na(tail(Sim$lambda, 1))){
          eta[i]     <- NA
          markup[i]  <- NA
          gamma[i]   <- NA
          alpha[i]   <- NA
        }else if ((tail(Sim$lambda, 1) >= 0.4) & (tail(Sim$omega, 1) >= 0.4) &
                  (tail(Sim$debt_share, 1) <= 2.7)){
          outcome[i] <- 1
          good.eta[i] <- eta[i]
          good.markup[i] <- markup[i]
          good.gamma[i] <- gamma[i]
          good.alpha[i] <- alpha[i]
          good.outcome[i] <- tail(Sim$lambda, 1)
        } else {
          outcome[i] <- 0
        } 
      }
    }
    save(eta,markup,gamma,alpha,outcome,good.eta,good.markup,
         good.gamma,good.alpha,good.outcome,file=filename)
  } else{
    # read in the data
    load(filename)
  }
  return(data.frame(eta,markup,gamma,alpha,outcome,good.eta,
                    good.markup,good.gamma,good.alpha,good.outcome))
}

#------------------------------------------------------------------------------
result <- categorize.result(n=1000,lambda_init=0.675,omega_init=0.578,
                            debt_share_init=1.53)

# Logistic regression
# Fit the model
model <- glm(outcome ~., data = na.omit(result[1:5]), family = binomial)
# Summarize the model
summary(model)

# scale vars
scaled_inputs <- scale(result[1:4])
data_scaled <- data.frame(eta = scaled_inputs[,1],
                          markup = scaled_inputs[,2],
                          gamma = scaled_inputs[,3],
                          alpha = scaled_inputs[,4],
                          outcome = result[5])

# Fit the model
model <- glm(outcome ~., data = na.omit(data_scaled), family = binomial)
# Summarize the model
logistic_res1 <- model
(logistic_res1S <- summary(model))

col1_1 <- logistic_res1S$coefficients[2:5,1]

#PRCC
inputs <- na.omit(result[6:9])
output <- na.omit(result[10])
pcc(inputs, output,rank=TRUE)
# again we see little effect of climate params
# scaled
# inputs_scaled <- as.data.frame(scale(inputs))
# pcc(inputs_scaled, output,rank=TRUE)
# # scaling doesn't affect PRCC

# for the table:
inputs <- na.omit(result[6:9])
output <- na.omit(result[10])
(prcc_res1 <- pcc(inputs, output,rank=TRUE,nboot=1000))
col2_1 <- prcc_res1$PRCC

# create table 1
table_1 <- cbind(col1_1,col2_1)
colnames(table_1) <- c("Logistic reg.", "PRCC")
table_1


