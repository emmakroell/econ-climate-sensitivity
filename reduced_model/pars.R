Parms <- c(
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
  div_max = 0.3
)