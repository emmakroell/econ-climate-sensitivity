# Logistic regression + PRCC
library('sensitivity')

# Set-up 
source('full_model/main.R')   # load parameters
source('full_model/pars.R')   # load parameters
source('full_model/funcs.R')  # load functions
source('full_model/sim.R')    # load simulation file
source('full_model/param_dist.R') 
parms <- Parms
set.seed(720)

#------------------------------------------------------------------------------
# function to compute PRCC
categorize.result <- function(n,lambda_init,omega_init,debt_share_init,climate){
  # n: number of model runs
  # lambda_init: initial employment rate
  # omega_init: initial wage share
  # debt_share_init: initial debt share
  filename <- sprintf("full_model/prcc_data/npts_%g_lambs_%g_omg_%g_d_%g_climate_%s.Rdata",
                      n, lambda_init, omega_init, debt_share_init,climate)
  
  if (!file.exists(filename)){
    # draw pricing parameters from distributions
    eta    <- replicate(n,draw_eta())
    markup <- replicate(n,draw_markup())
    gamma  <- replicate(n,draw_gamma())
    alpha  <- replicate(n,draw_alpha())
    ECS    <- replicate(n,draw_ECS())
    C_UP_preind <- replicate(n,draw_C_UP_preind())
    X <- data.frame(eta,markup,gamma,alpha,ECS,C_UP_preind)

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
    
   
    # Set up simulation
    Time <- c(
      start         =     2020, 
      end           =     2100,
      step          =     0.05
    )
    
    if (climate == TRUE){
      Options <- list(
        invest = 'lin',  # exp / arctan / lin
        damage_scenario = "10at4",
        subsidy = 0.5   # fraction of abatement costs subsidized by government
      )
      Parms['carbon_slope'] = 1
    } else{
      Options <- list(
        invest = 'lin',  # exp / arctan / lin
        damage_scenario = "None",
        subsidy = 0   # fraction of abatement costs subsidized by government
      )
      Parms['carbon_slope'] = 0
    }
    
    # arrays to store outcomes
    outcome <- rep(NA,n)
    good.outcome <- rep(NA,n)
    
    good.eta <- rep(NA,n)
    good.markup<- rep(NA,n)
    good.gamma <- rep(NA,n)
    good.alpha <- rep(NA,n)
    good.ECS <- rep(NA,n)
    good.C_UP <- rep(NA,n)

    for (i in 1:n){
      Parms[['eta_p']] = eta[i]
      Parms[['markup']] = markup[i]
      Parms[['gamma']] = gamma[i]
      Parms[['alpha']] = alpha[i]
      Parms[['S']] = ECS[i]
      Parms[['CO2_UP_preind']] = C_UP_preind[i]
      ## Run a scenario
      Sim <- try(simulation(time       = Time,
                            init_state = IC,
                            parms      = Parms,
                            options    = Options,
                            method     = 'lsoda'))
      if (class(Sim) == "try-error"){
        eta[i]     <- NA
        markup[i]  <- NA
        gamma[i]   <- NA
        alpha[i]   <- NA
        ECS[i]     <- NA
        C_UP_preind[i] <- NA
      } else{
        if (is.na(tail(Sim$lambda, 1))){
        eta[i]     <- NA
        markup[i]  <- NA
        gamma[i]   <- NA
        alpha[i]   <- NA
        ECS[i]     <- NA
        C_UP_preind[i] <- NA
      }else if ((tail(Sim$lambda, 1) >= 0.4) & (tail(Sim$omega, 1) >= 0.4) &
                (tail(Sim$debt_share, 1) <= 2.7)){
          outcome[i] <- 1
          good.eta[i] <- eta[i]
          good.markup[i] <- markup[i]
          good.gamma[i] <- gamma[i]
          good.alpha[i] <- alpha[i]
          good.ECS[i] <- ECS[i]
          good.C_UP[i] <- C_UP_preind[i]
          good.outcome[i] <- tail(Sim$lambda, 1)
        } else {
          outcome[i] <- 0
        } 
      }
    }
    save(eta,markup,gamma,alpha,ECS,C_UP_preind,outcome,good.eta,good.markup,
         good.gamma,good.alpha,good.ECS,good.C_UP,good.outcome,file=filename)
  } else{
    # read in the data
    load(filename)
  }
  return(data.frame(eta,markup,gamma,alpha,ECS,C_UP_preind,outcome,
                    good.eta,good.markup,good.gamma,good.alpha,good.ECS,
                    good.C_UP,good.outcome))
}

#------------------------------------------------------------------------------
# No feedback
result1 <- categorize.result(n=1000,lambda_init=0.675,omega_init=0.578,
                             debt_share_init=1.53,climate=FALSE)

# Logistic regression
# Fit the model
model <- glm(outcome ~., data = na.omit(result1[1:7]), family = binomial)
# Summarize the model
summary(model)

# scale vars
scaled_inputs <- scale(result1[1:6])
data_scaled <- data.frame(eta = scaled_inputs[,1],
                          markup = scaled_inputs[,2],
                          gamma = scaled_inputs[,3],
                          alpha = scaled_inputs[,4],
                          ECS = scaled_inputs[,5],
                          C_UP = scaled_inputs[,6],
                          outcome = result1[7])

# Fit the model
model <- glm(outcome ~., data = na.omit(data_scaled), family = binomial)
# Summarize the model

logistic_res2 <- model
(logistic_res2S <- summary(model))
# climate coeffs aren't statistically significant - good since they shouldn't
# affect the model with no feedback!
# alpha also not sign

col1_2 <- logistic_res2S$coefficients[2:7,1]

#PRCC
inputs <- na.omit(result1[8:13])
output <- na.omit(result1[14])
pcc(inputs, output,rank=TRUE)
# again we see little effect of climate params
# scaled
# inputs_scaled <- as.data.frame(scale(inputs))
# pcc(inputs_scaled, output,rank=TRUE)
# # scaling doesn't affect PRCC

# for the table:
inputs <- na.omit(result1[8:13])
output <- na.omit(result1[14])
(prcc_res2 <- pcc(inputs, output,rank=TRUE,nboot=1000))
col2_2 <- prcc_res2$PRCC

# create table 2
table_2 <- cbind(col1_2,col2_2)
colnames(table_2) <- c("Logistic reg.", "PRCC")
table_2


#------------------------------------------------------------------------------
# With climate and damages
result2 <- categorize.result(n=1000,lambda_init=0.675,omega_init=0.578,
                             debt_share_init=1.53,climate=TRUE)

# Logistic regression
# Fit the model
model <- glm(outcome ~., data = na.omit(result2[1:7]), family = binomial)
# Summarize the model
summary(model)

# scale vars
scaled_inputs <- scale(result2[1:6])
data_scaled <- data.frame(eta = scaled_inputs[,1],
                          markup = scaled_inputs[,2],
                          gamma = scaled_inputs[,3],
                          alpha = scaled_inputs[,4],
                          ECS = scaled_inputs[,5],
                          C_UP = scaled_inputs[,6],
                          outcome = result2[7])

# Fit the model
model <- glm(outcome ~., data = na.omit(data_scaled), family = binomial)
# Summarize the model

logistic_res3 <- model
(logistic_res3S <- summary(model))
col1_3 <- logistic_res3S$coefficients[2:7,1]
# gamma not statistically significant

#PRCC
inputs <- na.omit(result2[8:13])
output <- na.omit(result2[14])
(prcc_res3 <- pcc(inputs, output,rank=TRUE,nboot=1000))
col2_3 <- prcc_res3$PRCC

# create table 2
table_3 <- cbind(col1_3,col2_3)
colnames(table_3) <- c("Logistic reg.", "PRCC")
table_3

