# FILE FOR FIGURE 2
# Packages
library('here')
library('deSolve')
library('RColorBrewer')
library('qrng')
library('scatterplot3d')
library('tidyverse')
library('plotly')
library('rgl')
library('gMOIP')

# Source model code
source('reduced_model/pars.R')    # load parameters
source('reduced_model/funcs.R')   # load functions
source('reduced_model/sim.R')     # load simulation file 

# Set up simulation
Time <- c(
  start         =     2020, 
  end           =     2100,
  step          =     0.05
)

Options <- list(
  transform_vars  = FALSE  
)
#------------------------------------------------------------------------------
# function creates grid of initial conditions
create.ic.grid <- function(n_pts=20, type = 'grid'){
  if (type == 'grid'){
    lambda.grid  <- seq(0.2, 0.99, length.out = n_pts)   
    omega.grid   <- seq(0.2, 0.99, length.out = n_pts)  
    debt.grid    <- seq(0.1, 2.7, length.out = n_pts)   
    res <- expand.grid(x = lambda.grid, y = omega.grid, z = debt.grid)
    
  } else if (type == 'sobol'){
    var_ranges <- matrix(c(0.2,0.99,0.2,0.99,0.1,2.7),ncol=2, byrow=TRUE)
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

#------------------------------------------------------------------------------
# function computes basin of attraction for 3d model for a given number of points. 
# params eta, markup, and gamma can be specified
compute.basin.reduced <- function(n_pts, eta=0.192, markup = 1.18,
                             gamma=0.9, type='sobol', plot.res = TRUE){
  # name file to save data:
  savefile <- sprintf("reduced_model/basin_res/npts_%g_type_%s_eta_%g_mark_%g_gam_%g.Rdata",
                      n_pts, type, eta, markup, gamma)
  if (!file.exists(savefile)) {
    # set-up grid
    grid <- create.ic.grid(n_pts = n_pts, type = type)
    outcome <- rep(NA, nrow(grid))
    cputime.start <- proc.time() # start timer
    
    # set parameters
    Parms[['eta_p']] = eta
    Parms[['markup']] = markup
    Parms[['gamma']] = gamma
    
    # for loop
    for (i in seq(1,nrow(grid))){
      
      cat(i, 'out of', nrow(grid), '\n')
      
      # Set initial conditions 
      IC <- c(
        lambda = grid[i,1], 
        omega  = grid[i,2],
        debt   = grid[i,3],
        pop    = 4.83
      )
      
      
      # Run simulation
      Sim <- try(simulation_red(time = Time,
                               init_state = IC,
                               parms      = Parms,
                               #options    = Options,
                               method     = 'lsoda'))
      tsim <- tail(Sim,1)
      # if the simulation encounters an error, record this and keep going:
      if (inherits(Sim,"try-error") || is.na(tail(Sim$lambda, 1))) {
        outcome[i] = 'error'
      } else{
        # categorize as good / bad / error
        if ((0.4 <= tail(Sim$omega,1)) && (1 > tail(Sim$omega,1)) &&
            (0.4 <= tail(Sim$lambda,1)) && (1 > tail(Sim$lambda,1)) &&
            #(0.1 <= tail(Sim$debt_share,1)) & 
            (2.7 >= tail(Sim$debt_share,1))){
          outcome[i] = 'good'
        } else if ((tail(Sim$omega,1) > 1) | (tail(Sim$lambda,1) > 1)){
          outcome[i] = 'outside_bounds'
        } else{
          outcome[i] = 'bad'
        }}}
    cputime.end <- proc.time() # end timer
    cputime <- cputime.end - cputime.start
    cat('Time in min:',cputime/60)
    # Result saved in data frame with parameters
    result <- data.frame(grid,outcome)
    colnames(result) <- c("lambda", "omega", "debt", "outcome")
    # Save results
    save(result, cputime, file=savefile)
  } else { 
    # read in the data 
    load(savefile)
  }
  if (plot.res){
    par(mfrow = c(1,1), las=1, xpd=TRUE)
    colors <- c("#D95F02", "#1B9E77","#7570B3",'red')
    colors <- colors[as.numeric(result$outcome)]
    scatterplot3d(result[,1:3], pch = 16, color=colors,
                  main="Basin of attraction", xlab = 'lambda',
                  ylab = 'omega', zlab = 'debt')
    legend(x=0, y=-1, legend = levels(result$outcome),
           col = c("#D95F02", "#1B9E77","#7570B3",'red'), pch = 16,
           inset = -0.25, xpd = TRUE, horiz = TRUE)
  }
  return(result)
}

#------------------------------------------------------------------------------
# Interactive plot
# function generates interactive 3d scatter plot based on a result array
interactive.scatter <- function(result){
  df <- result %>%
    rownames_to_column() %>%
    as_tibble() %>%  ## as_data_frame() deprecated tibble 2.0
    mutate(outcome = as.factor(outcome))
  
  # Create the plot
  p <- plot_ly(
    df, x = ~lambda, y = ~omega, z = ~debt, size = 5,
    color = ~outcome, colors = c("#D95F02","#7570B3","#1B9E77")
  ) %>%
    add_markers() %>%
    layout(
      scene = list(xaxis = list(title = 'lambda'),
                   yaxis = list(title = 'omega'),
                   zaxis = list(title = 'd'))
    )
  return(p)
}

#------------------------------------------------------------------------------
# PLOTS
result_3d1 <- compute.basin.reduced(n_pts=20, markup=1.18, plot.res = FALSE)
result_3d1$outcome[result_3d1$outcome == 'outside_bounds']='error'
interactive.scatter(result_3d1)

result_3d2 <- compute.basin.reduced(n_pts=20, markup=1.3, plot.res = FALSE)
result_3d2$outcome[result_3d2$outcome == 'outside_bounds']='error'
interactive.scatter(result_3d2)

result_3d3 <- compute.basin.reduced(n_pts=20, markup=1.875, plot.res = FALSE)
result_3d3$outcome[result_3d3$outcome == 'outside_bounds']='error'
interactive.scatter(result_3d3)

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

## define colours/alpha so we can match legend to hull3d() defaults
c0 <- c("orange","purple","green")
a0 <- c(0.1,0.3,0.3)

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

## parallel color/alpha adjustment
set_alpha <- function(colvec, alphavec) {
    purrr::map2_chr(colvec, alphavec, ~adjustcolor(col=.x, alpha.f=.y))
}
set_alpha(c0,a0)

interactive.scatter(result_3d1)

## could also save this to an external (text or RDS?) file
p3dsave <- list(zoom = 0.907029569149017, userProjection = structure(c(1, 
0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), .Dim = c(4L, 4L
)), userMatrix = structure(c(0.811060905456543, 0.227074816823006, 
-0.539089322090149, 0, -0.584931075572968, 0.305393934249878, 
-0.751392066478729, 0, -0.00598765164613724, 0.924754917621613, 
0.38051626086235, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), scale = c(2.06813907623291, 
2.06813907623291, 0.628396093845367), FOV = 30)

## FIXME: wrap all of this into an appropriate function
## with arguments that let us do plots more easily
open3d()
## set up big window (actual pixel sizes matter for rendering!)
par3d(windowRect = c(0, 0, 1000, 1000))
## make 'reflective' component black (could try gray?)
material3d(specular="black")  ## tried lit=FALSE but don't like it
# plot the convex hull(s) (and/or points)
hull3d(result_3d1, do_points=list(FALSE,TRUE,FALSE),
       pointsize=7)
aspect3d(1,1,1)  ## reset to equal aspect ratio (seems like a good idea)
## restore saved projection info
if (exists("p3dsave")) par3d(p3dsave)

# change axes, add legend
## axis3d('x', pos = c(NA, 0, 0))
## https://stackoverflow.com/questions/29988230/using-expression-in-r-rgl-axis-labels
## Greek letter reference http://www.alanwood.net/demos/symbol.html#s0370
## ?? expression() seemed to be working until I messed around with
##  some other graphics settings ...
if (FALSE) {
    title3d(xlab=intToUtf8(955), ## expression(paste(lambda)),
            ylab=intToUtf8(969), ## expression(omega),
            zlab= 'd',
            line=1,
            cex=2)
}
## title3d isn't flexible enough.
## use it interactively to figure out which axes are which,
## then assign edges manually
## may want to make "which edges" arguments to a plotting fn
mtext3d(edge="x++",text=intToUtf8(955),line=2,cex=2)
mtext3d(edge="y--",text=intToUtf8(969),line=2,cex=2)
mtext3d(edge="z-+",text="d",line=2,cex=2)
legend3d("topright", legend = c('bad', 'error', 'good'), pch = 16,
         col = set_alpha(c0,a0),
         cex=2, inset=c(0.02))

## rotate/zoom until you're happy, then capture perspective info;
## then save/restore it above
if (FALSE) {
    p3dsave <- par3d()[c("zoom","userProjection","userMatrix","scale","FOV")]
    dput(p3dsave)
}

close3d()

#----------------------------------------------------------------------

###' @param x data frame/tibble with lambda, omega, d, outcome
###' @param colvec 
###' @param alphavec
###' @param do_points
###' @param save_image T/F save image as png in reduced_model/basin_res
convex_hull_plot <- function(x, markup, colvec=c("orange","purple","green"),
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
  mtext3d(edge="x++",text=expression(lambda),line=2,cex=2)
  mtext3d(edge="y--",text=expression(omega),line=2,cex=2)
  mtext3d(edge="z-+",text="d",line=2,cex=2)

  par3d(windowRect = c(100, 150, 800, 850))
  if (save_image){
    filename <- sprintf("markup_%g",markup) %>% 
      str_replace("[.]","_") %>% paste0(".png")
    rgl.snapshot(here("reduced_model/basin_res",filename))
  }
  close3d()
}

convex_hull_plot(result_3d1, markup = 1.18,do_points=list(FALSE,TRUE,FALSE),
                 save_image = FALSE)
convex_hull_plot(result_3d2, markup = 1.3, do_points=list(FALSE,TRUE,FALSE),
                 save_image = FALSE)
convex_hull_plot(result_3d3, markup = 1.875, do_points=list(FALSE,TRUE,FALSE),
                 save_image = FALSE)
