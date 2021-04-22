# create convex hull plots
# run after basin_3d.R


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
                   alphavec=c(0.3,0.3,0.3),
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
  # # if expression() doesn't work:
  # mtext3d(edge="x++",text=intToUtf8(955),line=2,cex=2)
  # mtext3d(edge="y--",text=intToUtf8(969),line=2,cex=2)
  mtext3d(edge="z-+",text="d",line=2,cex=2)

  par3d(windowRect = c(100, 150, 800, 850))
  if (save_image){
    filename <- sprintf("markup_%g",markup) %>%
      str_replace("[.]","_") %>% paste0(".png")
    rgl.snapshot(here("reduced_model/basin_res",filename))
  }
  close3d()
}

convex_hull_plot(result1, markup = 1.18,do_points=list(FALSE,TRUE,FALSE),
                 save_image = FALSE)
convex_hull_plot(result2, markup = 1.3, do_points=list(FALSE,TRUE,FALSE),
                 save_image = FALSE)
convex_hull_plot(result3, markup = 1.875, do_points=list(FALSE,TRUE,FALSE),
                 save_image = FALSE)