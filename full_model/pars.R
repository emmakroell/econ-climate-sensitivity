### Parameters
Parms <- c(
  # ECONOMIC PARAMTERS
  alpha = 0.02,      # productivity growth rate (Bovari et al. 2018a,b code)
  delta_N = 0.031,  # population growth rate (Bovari et al. 2018a,b)
  N_max = 7.056,     # population max (Bovari et al. 2018a,b)
  delta = 0.04,      # depreciation rate (Bovari et al. 2018a,b)
  nu = 2.7,          # capital to output ratio (Bovari et al. 2018a,b)
  r = 0.02,          # interest rate on loans (Bovari et al. 2018b)
  
  eta_p = 0.192,     # adjustment speed for prices (Bovari et al. 2018b)
  markup = 1.875,    # mark-up value (Bovari et al. 2018b)
  gamma = 0.9,       # inflation sensitivity in the bargaining equation (Grasselli estimate)
  
  # linear Phillips curve (Bovari et al. 2018a,b)
  phi_0 = -0.292,   
  phi_1 = 0.469,
  
  # linear investment curve (Bovari et al. 2018a,b)
  kappa_const = 0.0318,
  kappa_slope = 0.575,
  kappa_min = 0,
  kappa_max = 0.3,
  
  # linear dividend function (Bovari et al. 2018b)
  div_const = -0.078,
  div_slope = 0.553,
  div_min = 0,
  div_max = 0.3,
  
  #==========================================================================
  # CLIMATE PARAMETERS  (Bovari et al. 2018a,b)
  # C02 levels
  CO2_AT_preind = 588, # preindustrial concentration of CO2 in the 
  # atmosphere layer, in Gt C
  CO2_UP_preind = 360, # preindustrial concentration of CO2 in the 
  # bioshpere and upper ocean layer, in Gt C
  CO2_LO_preind = 1720, # preindustrial concentration of CO2 in the 
  # lower ocean layer, in Gt C
  phi_12 = 0.0239069, #0.024,  # Transfer coefficient for carbon from the atmosphere
  # to the upper ocean and biosphere
  phi_23 = 0.0013409, #0.001,  # Transfer coefficient for carbon from the upper ocean
  # and biosphere to the lower ocean
  
  # CO2 emissions
  delta_E_land = -0.022, # Growth rate of land use change CO2-e emissions
  delta_g_sigma = -0.001, # Variation rate of the growth of emission intensity
  
  # Radiative Forcing
  F2xCO2 = 3.681,  # Change in radiative forcing resulting from a doubling
  # of CO2 concentration from pind levels, in W/m^2
  F_exo_start = 0.5, # Value of exogenous radiative forcing in 2016
  F_exo_end = 1, # Value of exogenous radiative forcing in 2100
  F_exo_start_year = 2016,
  F_exo_end_year = 2100,
  
  # Temperature
  T_preind = 13.74,  # Preindustrial temperature, in degrees Celsius
  C_init = 1/0.098,  # Heat capacity of the atmosphere, biosphere, 
  # and upper ocean, in SI
  C0 = 3.52,  # Heat capacity of the deeper ocean
  gamma_star = 0.0176, # Heat exchange coefficient between temperature layers, in SI
  S = 3.1,  # Equilibrium climate sensitivity, in degrees Celsius
  #C = 49.76115, # Heat capacity
  
  #==========================================================================
  # DAMAGE PARAMETERS
  # from the literature
  dam_exp       = 6.754,       # Parameter for damage function - Weitzman (2012) 
  dam_1         = 0,           # Damage function parameter, in  per degrees Celsius - DICE model
  dam_2         = 0.00236,     # Damage function parameter, in  per degrees Celsius squared - DICE model
  dam_weitzman  = 0.00000507,  # Damage function parameter, in  per degrees Celsius 
  # to the power of damage_exp - Weitzman (2012)
  dam_stern     = 0.0000819, # Damage function parameter, in  per degrees Celsius 
  # to the power of damage_exp - Deitz and Stern (2015)
  
  # calculated damage curves parameters - inverse polynomial of degree 7
  dam_ten     = 4.476997e-06,   # calculated to produce 10% damages at 4 degrees C
  dam_twenty  = 1.29541e-05,    # calculated to produce 20% damages at 4 degrees C
  dam_thirty  = 2.385324e-05,   # calculated to produce 30% damages at 4 degrees C
  dam_forty   = 3.838542e-05,   # calculated to produce 40% damages at 4 degrees C
  dam_fifty   = 5.873047e-05,   # calculated to produce 50% damages at 4 degrees C

  #==========================================================================
  # CARBON PRICE PARAMETERS 
  theta = 2.6,        # Abatement cost function parameter (Bovari et al. 2018a,b)
  g_p_BS  = -0.0051,  # Exogenous growth rate of the price of backstop technology  (Bovari et al. 2018a,b)
  beta_C  = 0.0124,   #(Bovari et al. 2018a,b)
  S_ref   = 2.9,      #(Bovari et al. 2018a,b)
  carbon_parm = 0.02,  # Nordhaus path
  carbon_slope = 2,
  conv10to15     = 1.160723971/1000,   # currency conversion #(Bovari et al. 2018a,b code/private correspondence)
  p_Car_step_year_1 = 2020,
  p_Car_step_year_2 = 2030,
  p_Car_val_1 = 1.05543*1.05405*1.05494*1.07863*1.02471*80,  # conversion from AFD code
  p_Car_val_2 = 1.05543*1.05405*1.05494*1.07863*1.02471*100
)