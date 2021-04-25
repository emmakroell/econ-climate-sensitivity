# create convex hull plots
# run after basin_full.R

source('full_model/basin_full.R')  
library('rgl')
library('gMOIP')


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
c0 <- c("orange","green")
a0 <- c(0.1,0.3)

###' @param x data frame/tibble with lambda, omega, d, outcome
###' @param colvec
###' @param alphavec
###' @param do_points
hull3d <- function(x, colvec=c("orange","purple","green"),
                   alphavec=c(0.1,0.3),
                   do_points=rep(FALSE, 2),
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



interactive.scatter(result1 %>% filter(year == 2300))

## could also save this to an external (text or RDS?) file
p3dsave <- list(zoom = 0.907029569149017, userProjection = structure(c(1,
                                                                       0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), .Dim = c(4L, 4L
                                                                       )), userMatrix = structure(c(0.811060905456543, 0.227074816823006,
                                                                                                    -0.539089322090149, 0, -0.584931075572968, 0.305393934249878,
                                                                                                    -0.751392066478729, 0, -0.00598765164613724, 0.924754917621613,
                                                                                                    0.38051626086235, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)), scale = c(2.06813907623291,
                                                                                                                                                                   2.06813907623291, 0.628396093845367), FOV = 30)



#----------------------------------------------------------------------

###' @param x data frame/tibble with lambda, omega, d, outcome
###' @param colvec
###' @param alphavec
###' @param do_points
###' @param save_image T/F save image as png in reduced_model/basin_res
convex_hull_plot <- function(x, markup, colvec=c("orange","purple","green"),
                             alphavec=c(0.3,0.3,0.3),
                             do_points=rep(FALSE, 3),
                             save_image = TRUE,...) {
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
  # # if expression() doesn't work:
  # mtext3d(edge="x++",text=intToUtf8(955),line=2,cex=2)
  # mtext3d(edge="y--",text=intToUtf8(969),line=2,cex=2)
  mtext3d(edge="z-+",text="d",line=2,cex=2)
  
  par3d(windowRect = c(100, 150, 800, 850))
  if (save_image){
    filename <- sprintf("markup_%g_full",markup) %>%
      str_replace("[.]","_") %>% paste0(".png")
    rgl.snapshot(here("full_model/basin_res",filename))
  }
  close3d()
}

result1_to_plot <- result1 %>% filter(year == 2300) %>% select(lambda.ic,omega.ic,debt.ic,outcome)
result2_to_plot <- result2 %>% filter(year == 2300) %>% select(lambda.ic,omega.ic,debt.ic,outcome)
result3_to_plot <- result3 %>% filter(year == 2300) %>% select(lambda.ic,omega.ic,debt.ic,outcome)



convex_hull_plot(result1_to_plot, markup = 1.18,do_points=list(FALSE,FALSE),
                 colvec=c("orange","green"),alphavec=c(0.3,0.3),save_image = FALSE)
convex_hull_plot(result2_to_plot, markup = 1.3, do_points=list(FALSE,FALSE),
                 colvec=c("orange","green"),alphavec=c(0.3,0.3),save_image = FALSE)
convex_hull_plot(result3_to_plot, markup = 1.875, do_points=list(TRUE,FALSE),
                 colvec=c("orange","green"),alphavec=c(0.3,0.3),save_image = FALSE)
