#libraries
library('here')
library('deSolve')
library('RColorBrewer')
library('qrng')
library('scatterplot3d')
library('tidyverse')
library('plotly')
library('rgl')
library('gMOIP')

# Colour palettes for graphs
greens <- brewer.pal(n = 6, name = "Greens")
colourful <- brewer.pal(n = 3, name = "Set1")
colour2 <- brewer.pal(n = 3, name = "Set2")

source('full_model/pars.R')     # load parameters
source('full_model/funcs.R')    # load functions
source('full_model/sim.R')      # load simulation file

# Set up simulation
Time <- c(
  start         =     2020, 
  end           =     2100,
  step          =     0.05
)

Options <- list(
  invest = 'lin',  # exp / arctan / lin
  damage_scenario = '10at4',
  subsidy = 0.5   # fraction of abatement costs subsidized by government
)
#================================================================================
create.parm.grid <- function(n_pts=20, type = 'sobol'){
  if (type == 'grid'){
    eta.grid    <- seq(0, 4, length.out = n_pts)   # eta
    markup.grid <- seq(1, 3, length.out = n_pts)   # markup
    gamma.grid  <- seq(0, 1, length.out = n_pts)   # gamma
    res <- expand.grid(x = eta.grid, y = markup.grid, z = gamma.grid)
    
  } else if (type == 'uniform'){
    eta.grid    <- 4 * runif(n_pts)   # eta
    markup.grid <- 1 + 2 * runif(n_pts)   # markup
    gamma.grid  <- runif(n_pts)   # gamma
    res <- expand.grid(x = eta.grid, y = markup.grid, z = gamma.grid)
    
  } else if (type == 'sobol'){
    var_ranges <- matrix(c(0,1,1,2.5,0,1),ncol=2, byrow=TRUE)
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

#================================================================================
explore.parm.space <- function(n_pts, type, plot.res = FALSE, 
                               lambda_init=0.9,omega_init=0.9,
                               debt_init = 0.3){
  # name file to save data:
  savefile <- sprintf("full_model/parms_res/npts_%g_type_%s_lambs_%g_omg_%g_d_%g.Rdata",
                      n_pts, type, lambda_init, omega_init, debt_init)
  if (!file.exists(savefile)) {
    # Set initial conditions 
    lambda_init <-lambda_init
    omega_init <-   omega_init
    debt_share_init <-debt_init
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
    # set-up grid
    grid <- create.parm.grid(n_pts = n_pts, type = type)
    outcome <- rep(NA, nrow(grid))
    cputime.start <- proc.time() # start timer
    
    # for loop
    for (i in seq(1,nrow(grid))){
      cat(i, 'out of', nrow(grid), '\n')
      
      Parms[['eta_p']] = grid[i,1]
      Parms[['markup']] = grid[i,2]
      Parms[['gamma']] = grid[i,3]
      
      Sim <- simulation(time       = Time,
                        init_state = IC,
                        parms      = Parms,
                        options    = Options,
                        method     = 'lsoda')
      
      # if the simulation encounters an error, record this and keep going:
      if (inherits(Sim,"try-error")){
        outcome[i] = 'error'
      } else if (is.na(tail(Sim$lambda, 1)) == TRUE){
        outcome[i] = 'error'
      }
      else{
        # categorize as good / bad / error
        if ((0.4 <= tail(Sim$omega,1)) & (1 > tail(Sim$omega,1)) &
            (0.4 <= tail(Sim$lambda,1)) & (1 > tail(Sim$lambda,1)) &
            (2.7 >= tail(Sim$debt_share,1))){
          outcome[i] = 'good'
        } else if ((tail(Sim$omega,1) > 1) | (tail(Sim$lambda,1) > 1)){
          outcome[i] = 'outside_bounds'
        } else{
          outcome[i] = 'bad'
        }}}
    
    cputime.end <- proc.time() # end timer
    cputime <- cputime.end - cputime.start
    # Result saved in data frame with parameters
    result <- data.frame(grid,outcome)
    colnames(result) <- c("eta", "markup", "gamma", "outcome")
    # Save results
    save(result, cputime, file=savefile)
  } else { 
    # read in the data 
    load(savefile)
  }
  if (plot.res == TRUE){
    par(mfrow = c(1,1), las=1, xpd=T)
    colors <- c("#D95F02", "#1B9E77","#7570B3",'red')
    colors <- colors[as.numeric(result$outcome)]
    scatterplot3d(result[,1:3], pch = 16, color=colors,
                  main="Parameter grid search", xlab = '',
                  ylab = '', zlab = '')
    legend(x=0, y=-1, legend = levels(result$outcome),
           col = c("#D95F02", "#1B9E77","#7570B3",'red'), pch = 16,
           inset = -0.25, xpd = TRUE, horiz = TRUE)
    text(x = 7, y = 0.2, expression(xi), srt = 45)
    text(x = -0.8, y = 2.7, expression(gamma), srt = 0)
    text(x = 2.5, y = -1, expression(eta), srt = 0)
  }
  return(result)
}
#------------------------------------------------------------------------------
# Interactive plot
# function generates interactive 3d scatter plot based on a result array
interactive.scatter <- function(result){
  df <- result %>%
    rownames_to_column() %>%
    as_tibble() %>%
    mutate(outcome = as.factor(outcome))
  
  # Create the plot
  p <- plot_ly(
    df, x = ~eta, y = ~markup, z = ~gamma, size = 5,
    color = ~outcome, colors = c("#D95F02","#7570B3","#1B9E77")
  ) %>%
    add_markers() %>%
    layout(
      scene = list(xaxis = list(title = 'eta'),
                   yaxis = list(title = 'markup'),
                   zaxis = list(title = 'gamma'))
    )
  return(p)
}

#------------------------------------------------------------------------------
# PLOTS
result1 <- explore.parm.space(n_pts=20, lambda_init = 0.9, omega_init=0.9,
                              debt_init=0.3, type="sobol")
result1$outcome[result1$outcome == 'outside_bounds']='error'
interactive.scatter(result1)

result2 <- explore.parm.space(n_pts=20, lambda_init = 0.5, omega_init=0.6,
                              debt_init=1, type="sobol")
result2$outcome[result2$outcome == 'outside_bounds']='error'
interactive.scatter(result2)

result3 <- explore.parm.space(n_pts=20, lambda_init = 0.675, omega_init=0.578,
                              debt_init=1.53, type="sobol")
result3$outcome[result3$outcome == 'outside_bounds']='error'
interactive.scatter(result3)

#------------------------------------------------------------------------------
###' @param x a 3-column coordinate matrix
###' @param color colour of 3D hull/points
###' @param alpha transparency (0-1)
###' @param do_points plot separate points rather than hull, e.g. for small/disconnected sets
pfun <- function(x, color, alpha, do_points=FALSE, pointsize=NA) {
  ## cat(color,alpha,"\n")
  if (do_points) {
    rgl::points3d(x[,1],x[,2],x[,3],color=color, alpha=alpha,
                  size=pointsize)
  } else {
    gMOIP::plotHull3D(x, drawLines=FALSE,
                      argsPolygon3d=list(color=color,
                                         alpha=alpha))
  }
}

###' @param x data frame/tibble with lambda, omega, d, outcome
###' @param colvec
###' @param alphavec
###' @param do_points
hull3d <- function(x, colvec=c("orange","purple","green"),
                   alphavec=c(0.1,0.3,0.3),
                   do_points=rep(FALSE, 3),
                   ...) {
  m1 <- x %>% split(.$outcome) %>% map(~as.matrix(.[,1:3]))
  pmap(c(list(m1,colvec,alphavec, do_points), list(...)),  pfun)
  rgl::bbox3d()
  rgl::axes3d()
  invisible(NULL)
}


#----------------------------------------------------------------------
# from inspection:
p3dsave <- list(zoom = 1, userProjection = structure(c(1, 0, 0, 0, 0, 1,
                                                       0, 0, 0, 0, 1, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)),
                userMatrix = structure(c(-0.763124227523804,-0.212693452835083,
                                         0.61024808883667, 0, 0.646242499351501, -0.256196916103363,
                                         0.718841850757599, 0, 0.00345077319070697, 0.942934095859528, 
                                         0.332961738109589, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)),
                scale = c(1.19023811817169,0.793492019176483, 1.19023811817169),
                FOV = 30)


###' @param x data frame/tibble with xi, eta, gamma, outcome
###' @param colvec 
###' @param alphavec
###' @param do_points
###' @param save_image T/F save image as png in full_model/parms_res
convex_hull_plot <- function(x, name, colvec=c("orange","purple","green"),
                             alphavec=c(0.1,0.3,0.3),
                             do_points=rep(FALSE, 3),
                             save_image = TRUE,
                             ...) {
  ## parallel color/alpha adjustment
  set_alpha <- function(colvec, alphavec) {
    purrr::map2_chr(colvec, alphavec, ~adjustcolor(col=.x, alpha.f=.y))
  }
  open3d()
  ## set up big window (actual pixel sizes matter for rendering!)
  par3d(windowRect = c(0, 0, 1000, 1000))
  ## make 'reflective' component black (could try gray?)
  material3d(specular="black")  ## tried lit=FALSE but don't like it
  # plot the convex hull(s) (and/or points)
  hull3d(x,colvec = colvec, alphavec = alphavec,
         do_points = do_points, pointsize=7)
  aspect3d(1,1,1)  ## reset to equal aspect ratio (seems like a good idea)
  ## restore saved projection info
  if (exists("p3dsave")) par3d(p3dsave)
  
  # axes labels
  mtext3d(edge="x+-",text=expression(eta),line=2,cex=2)
  mtext3d(edge="y+-",text=expression(xi),line=2,cex=2)
  mtext3d(edge="z+-",text=expression(gamma),line=2,cex=2)
  # # if expression() doesn't work:
  # mtext3d(edge="x+-",text=intToUtf8(951),line=2,cex=2)
  # mtext3d(edge="y+-",text=intToUtf8(958),line=2,cex=2)
  # mtext3d(edge="z+-",text=intToUtf8(947),line=2,cex=2)
  
  # narrow window
  par3d(windowRect = c(100, 150, 800, 850))
  
  if (save_image){
    filename <- paste0(name,".png")
    rgl.snapshot(here("full_model/parms_res",filename))
  }
  close3d()
}

convex_hull_plot(result1, save_image=FALSE, name = "good", 
                 do_points=list(FALSE,TRUE,FALSE))
convex_hull_plot(result2, save_image=FALSE, name = "bad", 
                 do_points=list(FALSE,TRUE,FALSE))


