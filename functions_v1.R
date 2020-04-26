#
#
#

# x = snapShots[[i]]
tMax = 1
# hcols=hyCols
plot_hyphae = function(x, tMax=100, hcols=rev(brewer.pal(11, "Spectral")), ...) {
  segCol = data.frame(colorRamp(hcols, alpha=0.5)(x[,"t"]/max(tMax,1)))
  colnames(segCol) = c("red", "green", "blue", "alpha")
  segCols = rgb(red=segCol$red, green=segCol$green, blue=segCol$blue, maxColorValue = 255)
  plot(NA, xlab="", ylab="", ...)
  for(i in 1:dim(x)[1])
    segments(x$x0[i], x$y0[i], x$x[i], x$y[i], col=segCols[i])
  #points(x$x0, x$y0, col = segCols, pch = 20, cex=0.3)
}

#
# Helper function for hyphal length
#
hyphal_length = function(x) {
  return( sqrt( (x["x"]-x["x0"])^2 + (x["y"]-x["y0"])^2 ) )
}
seg_length = function(x) {
  return( sqrt( (x[3]-x[1])^2 + (x[4]-x[2])^2 ) )
}

corners = function(x, w) {
  u = c(x[1]-w, x[2]-w)
  v = c(x[1]+w, x[2]+w)
  return(t(matrix(c(u,v),2,2)))
}

insideRAE = function(x, r) {
  xSat = r[1,1] <= x[1] & x[1] <= r[2,1]
  ySat = r[1,2] <= x[2] & x[2] <= r[2,2]
  return(xSat & ySat)
}

insideBoundingBox = function(x, bb) {
  x = as.numeric(x)
  bb = as.numeric(bb)
  xSat = bb[1] <= x[1] & x[1] <= bb[3]
  ySat = bb[2] <= x[2] & x[2] <= bb[4]
  return(xSat & ySat)
}

hitsBB = function(bb, x, full=FALSE) {
  stIn = insideBoundingBox(x=x[1:2], bb=bb)
  enIn = insideBoundingBox(x=x[3:4], bb=bb)
  s1 = c(bb[c(3,2)], bb[c(1,2)])
  s2 = c(bb[c(1,2)], bb[c(1,4)])
  s3 = c(bb[c(1,4)], bb[c(3,4)])
  s4 = c(bb[c(3,4)], bb[c(3,2)])
  botX = doSegmentsIntersect(segment1=s1, segment2=x)
  lefX = doSegmentsIntersect(segment1=s2, segment2=x)
  topX = doSegmentsIntersect(segment1=s3, segment2=x)
  rigX = doSegmentsIntersect(segment1=s4, segment2=x)
  r = c(b=botX, l=lefX, t=topX, r=rigX, s=stIn, e=enIn)
  if(!full)
    r = sum(r) > 0
  return(r)
}

##############################################
# Map the RAE hit by each hyphae
hyphae_hits = function(hl, bbs) {
  m = length(hl)
  h2b = matrix(0, m, dim(bbs)[1])
  for(j in 1:m) {
    if(j %% 1000 == 0) print(j)
    hi = hl[[j]][1:4]
    hiBB = getBoundingBox(P0=hi[1:2], P1=hi[3:4])

    xSAT = (hiBB[1] <= bbs[,3]) & (bbs[,1] <= hiBB[3])
    ySAT = (hiBB[2] <= bbs[,4]) & (bbs[,2] <= hiBB[4])
    bbInds = which(xSAT & ySAT)
    bbIndsHits = bbInds[apply(bbs[bbInds,,drop=FALSE], 1, hitsBB, x=hi)]
    h2b[j, bbIndsHits] = 1
  }
  return(h2b)
}

##############################################
#
grid_bounding_boxes = function(w=10, xrng=c(-50,50), yrng=c(-50,50)) {
  x0 = seq(xrng[1], xrng[2], w)
  y0 = seq(yrng[1], yrng[2], w)
  centers = cbind(rep(x0, length(y0)), rep(y0, each=length(x0)))
  bbs = cbind(centers-(w/2), centers+(w/2))
  return(bbs)
}


RAEintersection <- function(m, b, side, bb){
  if(length(side) == 1){
    if(side %in% c("r", "l")){
      x = bb[side]
      p = c(x, m*x+b)
    }
    else {
      y = bb[side]
      if(is.finite(m)) {
        p = c((y-b)/m, y)
      }
      else {
        p = c(b, y)
      }
    }
  }
  return(p)
}

###########################################################
# Calculate the density in each RAE

hyphal_length_by_RAE = function(hl, h2bbs, bbs, plotting=FALSE) {
  hPerBox = colSums(h2bbs)

  d = array(0, dim(bbs)[1]) # density for each RAE
  for(j in which(hPerBox>0) ) { # All RAE that have one or more hyphae
    if(j %% 10 ==0) print(j)
    hInds = which(h2bbs[,j] > 0)
    if(plotting) polygon(x=bbs[j, c(1,1,3,3)], y=bbs[j,c(2,4,4,2)])

    if(length(hInds) > 0) { ## should be
      #segments(x0=hi[1], y0=hi[2], x1=hi[3], y1=hi[4], col=rgb(0.6,0.6,0.6,0.3))
      # bbj = getBoundingBox(P0=bbs[j,1:2], P1=bbs[j,3:4]) # same as next line
      bbj = bbs[j, ]
      names(bbj) = c("l", "b", "r", "t")
      for(i in hInds) {
        hi = hl[[i]][1:4]
        fl = hitsBB(bb=bbs[j,], x=hi, full=TRUE)

        sel = sum(fl[5:6])
        #print(sel)
        if(sel == 2) { # both points in RAE
          d[j] = d[j] + hl[[i]]["l"]
          if(plotting) segments(x0=hi[1], y0=hi[2], x1=hi[3], y1=hi[4], lty=1)
        }
        else {
          # linear fit to hyphae line
          bm = c(NA, NA)
          if(abs(hi[3]-hi[1]) > 0) { ## If not a vertical segment
            bm[2] = (hi[4] - hi[2]) / (hi[3] - hi[1])
            bm[1] = hi[4] - hi[3]*bm[2]
          }
          else
            bm[1] = hi[1] # slope undefined (vertical segment). Just set b to x

          # Reset either hi(x0,y0) or hi(x,y) to the coordinate where
          # the hyphal segment and the border(s) intersect
          if(sel == 1) { # Solve one intersection
            side = names(which(fl[1:4]))[1]
            xy = RAEintersection(m=bm[2], b=bm[1], side = side, bb = bbj)
            if(fl["s"])
              hi[3:4] = xy
            else
              hi[1:2] = xy
            if(plotting) segments(x0=hi[1], y0=hi[2], x1=hi[3], y1=hi[4], lty=2)
          }
          else { # Solve 2 intersections (whole hyphae crosses this RAE)
            sides = names(which(fl[1:4]))
            # if line hits a corner/vertex, then 2 sides must be non-adjacent
            if(length(sides)==3) {
              print(paste("Warning: hyphae", i, "hits a corner"))
              if(sum(fl[c(1,3)]) == 2) sides = c("b","t")
              else sides = c("l","r")
            }
            # reset hi to be the coordinates where the hyphae line and the borders intersects
            hi[1:2] = RAEintersection(m=bm[2], b=bm[1], side = sides[1], bb = bbj)
            hi[3:4] = RAEintersection(m=bm[2], b=bm[1], side = sides[2], bb = bbj)
          }
          if(sum(!is.finite(hi))>0) { print(paste("Non-finite hi", i, hi[3]-hi[1], hi[4]-hi[2])) }
          d[j] = d[j] + seg_length(hi)
          if(plotting) segments(x0=hi[1], y0=hi[2], x1=hi[3], y1=hi[4], lty=3)
        }
      }
    }
  }
  return(d)
}

##############################################
tipExtension <- function(ktip1, ktip2, Kt, l) {
  # source: Lejeune et al 1995, Morphology of Trichoderma reesei QM 9414 in Submerged Cultures
  return(ktip1+ktip2*(l/(l+Kt)))
}

tipExtensionMonod <- function(ktip1, ktip2, Kt, l, S, Ks) {
  # source: Lejeune et al 1995, Morphology of Trichoderma reesei QM 9414 in Submerged Cultures
  return( ( ktip1+ktip2*(l/(l+Kt)) ) * S/(S+Ks) )
}

perpendicularDistance <- function(x, xc, yc, R){
  x1 = x[, "x0"]
  x2 = x[, "x"]
  y1 = x[, "y0"]
  y2 = x[, "y"]

  d = (abs((y2-y1)*xc-(x2-x1)*yc+x2*y1-y2*x1))/(sqrt((y2-y1)**2+(x2-x1)**2))

  return(unique(c(which(d < R), which(d==R))))

  # return(d <= R)
}

hyphae_hits_substrate <- function(hl, bbs){
  m = length(hl)
  h2b = matrix(0, m, dim(bbs)[1])
  for (j in 1:m) {
    hi = hl[[j]][c(3, 4)]   #x,y value
    xSAT = (hi[1] <= bbs[,3]) & (bbs[,1] <= hi[1]) 
    ySAT = (hi[2] <= bbs[,4]) & (bbs[,2] <= hi[2])
    bbInds = which(xSAT & ySAT)
    
    h2b[j, bbInds] = 1
    
  }
  return(h2b)
  
}
