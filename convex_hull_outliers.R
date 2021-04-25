# FILE TO GENERATE FIGURE 1
# sources "reduced_model/parm_space_3d.R" and "full_model/parm_space_full.R"
# libraries:
library('tidyverse')
library('rgl')
library('gMOIP')

#------------------------------------------------------------------------------
# FUNCTIONS:
# convex hull with some errors separated as outliers
hull3d_alt <- function(x, colvec=c("orange", "purple", "purple", "green"),
                       alphavec=c(0.1,0.3,0.3,0.3),
                       do_points=c(FALSE, FALSE, TRUE, FALSE),
                       markup_bound, eta_bound, gamma_bound,
                       ...) {
  # select errors for solid:
  error1 <- x %>% filter(outcome == "error", markup > markup_bound,
                         eta > eta_bound, gamma < gamma_bound)
  # remaining errors for outliers:
  error2 <- x %>% filter(outcome == "error",
                         markup <= markup_bound | eta <= eta_bound | 
                           gamma >= gamma_bound) %>% 
    mutate(outcome = "error_outlier")
  # combine
  m <- rbind(filter(x, outcome=='good'), filter(x, outcome=='bad'),
             error1, error2)
  m1 <- m %>% split(.$outcome) %>% map(~as.matrix(.[,1:3]))
  pmap(c(list(m1,colvec,alphavec, do_points), list(...)),  pfun)
  rgl::bbox3d()
  rgl::axes3d()
  invisible(NULL)
}

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
###' @param markup_bound lower bound for error solid in markup
###' @param eta_bound lower bound for error solid in eta
###' @param gamma_bound upper bound for error solid in gamma
###' @param colvec 
###' @param alphavec
###' @param do_points
###' @param save_image T/F save image as png in FINAL/reduced_model/parms_res
###' @param name name for image file
convex_hull_plot <- function(x, markup_bound=1.6, eta_bound=0.2, gamma_bound = 0.8,
                             colvec=c("orange","purple","purple","green"),
                             alphavec=c(0.3,0.3,0.3,0.3),
                             do_points=c(FALSE, FALSE, TRUE, FALSE),
                             save_image = FALSE, name = NULL,
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
  #browser()
  hull3d_alt(x,colvec = colvec, alphavec = alphavec,
         do_points = do_points, pointsize=7, markup_bound=markup_bound,
         eta_bound = eta_bound, gamma_bound = gamma_bound)
  aspect3d(1,1,1)  ## reset to equal aspect ratio (seems like a good idea)
  ## restore saved projection info
  if (exists("p3dsave")) par3d(p3dsave)
  
  # axes labels
  mtext3d(edge="x+-",text=expression(bar(eta)),line=2,cex=2)
  mtext3d(edge="y+-",text=expression(bar(xi)),line=2,cex=2)
  mtext3d(edge="z+-",text=expression(bar(gamma)),line=2,cex=2)
  
  # narrow window
  par3d(windowRect = c(100, 150, 800, 850))
  
  if (save_image){
    filename <- paste0(name,".png")
    rgl.snapshot(filename)
  }
  close3d()
}

#------------------------------------------------------------------------------
# GENERATE PLOTS FOR FIGURE 1

# reduced model
source("reduced_model/parm_space_3d.R")
unique(result1$outcome)
result1_to_plot <- result1 %>% filter(year == 2300) %>% 
  select(eta,markup,gamma,outcome) %>% 
  mutate(outcome = ifelse(outcome == "outside_bounds","error",outcome))

convex_hull_plot(result1_to_plot, save_image=TRUE, name = "reduced_good_ic")

unique(result2$outcome)
result2_to_plot <- result2 %>% filter(year == 2300) %>% 
  select(eta,markup,gamma,outcome) %>% 
  mutate(outcome = ifelse(outcome == "outside_bounds","error",outcome))

convex_hull_plot(result2_to_plot, save_image=TRUE, name = "reduced_bad_ic")

# full model
source("full_model/parm_space_full.R")
unique(result3$outcome)
result3_to_plot <- result3 %>% filter(year == 2300) %>% 
  select(eta,markup,gamma,outcome) %>% 
  mutate(outcome = ifelse(outcome == "outside_bounds","error",outcome))

convex_hull_plot(result3_to_plot, save_image=TRUE, name = "full_good_ic")

unique(result4$outcome)
result4_to_plot <- result4 %>% filter(year == 2300) %>% 
  select(eta,markup,gamma,outcome) %>% 
  mutate(outcome = ifelse(outcome == "outside_bounds","error",outcome))
convex_hull_plot(result4_to_plot, save_image=TRUE, name = "full_bad_ic")




