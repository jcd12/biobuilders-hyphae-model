
# R functions to determine if two line segments intersect
# Note: this is very different question from whether
# two lines intersect (terminology is important!)

# These functions are based upon the excellent discussion here,
# written by Martin Thoma:
# martin-thoma.com/how-to-check-if-two-line-segments-intersect/
# The starting point was the code here:
# github.com/MartinThoma/algorithms/tree/master/crossingLineCheck/Geometry/src
# though I did not use everything in that directory.
# If found this via this link, which has other useful discussion:
# stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

# Most of these function have names identical or similar to the
# java originals.  The java was ported to R and tested using
# the graphical test included here.  It seems to work.
# Released under GPL-3 (OK'd by Martin)

# Prof. Bryan A. Hanson
# DePauw University, Greencastle IN 46135
# April 2013

getBoundingBox <- function(P0, P1) {

  # P0, P1 each have c(x,y)
  # ll = lower left
  # ur = upper right

  llx <- min(P0[1], P1[1])
  lly <- min(P0[2], P1[2])
  urx <- max(P0[1], P1[1])
  ury <- max(P0[2], P1[2])

  bb <- c(llx, lly, urx, ury)
}

doBoxesIntersect <- function(box1, box2) {

  ans <- FALSE
  chk1 <- box1[1] <= box2[3]
  chk2 <- box1[3] >= box2[1]
  chk3 <- box1[2] <= box2[4]
  chk4 <- box1[4] >= box2[2]
  if (chk1 && chk2 && chk3 && chk4) ans <- TRUE
  ans
}

isPointOnLine <- function(segment, point) {
  # segment is c(x1, y1, x2, y2)
  # translate segment to origin
  newseg <- c(0.0, 0.0, segment[3] - segment[1], segment[4] - segment[2])
  newpt <- c(point[1] - segment[1], point[2] - segment[2])
  # calc a modified cross product:
  # a.x * b.y - b.x * a.y
  # if zero, point is on segment
  # basically, you have two vectors sharing 0,0 as one end
  ans <- newseg[3]*newpt[2] - newpt[1]*newseg[4]
  return(isTRUE(all.equal(abs(ans), 0)))
}

isPointRightOfLine <- function(segment, point) {
  # see notes in isPointOnLine
  newseg <- c(0.0, 0.0, segment[3] - segment[1], segment[4] - segment[2])
  newpt <- c(point[1] - segment[1], point[2] - segment[2])
  ans <- newseg[3]*newpt[2] - newpt[1]*newseg[4]
  return(ans < 0)
}

lineSegmentTouchesOrCrossesLine <- function(segment1, segment2) {
  # segments given as c(x1, y1, x2, y2)
  ans <- 	(isPointOnLine(segment1, segment2[1:2]) ||
             isPointOnLine(segment1, segment2[3:4]) ||
             xor(isPointRightOfLine(segment1, segment2[1:2]),
                 isPointRightOfLine(segment1, segment2[3:4])))
  return(ans)
}

doSegmentsIntersect <- function(segment1, segment2) {
  # segments given as c(x1, y1, x2, y2)
  box1 <- getBoundingBox(segment1[1:2], segment1[3:4])
  box2 <- getBoundingBox(segment2[1:2], segment2[3:4])
  return(doBoxesIntersect(box1, box2) &&
           lineSegmentTouchesOrCrossesLine(segment1, segment2) &&
           lineSegmentTouchesOrCrossesLine(segment2, segment1))
}

# Graphical Test using random line segments

intersectTest <- function() {
  # Creates two random line segments and checks to see
  # if they intersect, shows result graphically
  x <- runif(4, -5, 5)
  y <- runif(4, -5, 5)
  xy <- data.frame(x = x, y = y) # points are in rows

  # Set up a plotting region, then add points
  plot(-10:10, -10:10, type = "n")
  points(xy, cex = 0.5)
  # Next label each segment
  text(xy[seq(1, nrow(xy), by = 2),], labels = 1:(nrow(xy)/2), pos = 4, cex = 0.75)

  s <- seq(1, nrow(xy), by = 2)

  for (n in s) { # draw the bounding boxes
    segments(xy[n,1], xy[n,2], xy[n+1,1], xy[n+1,2], col = "red")
    polygon(x = c(xy[n,1], xy[n,1], xy[n+1,1], xy[n+1,1]),
            y = c(xy[n,2], xy[n+1,2], xy[n+1,2], xy[n,2]))
  }

  for(n in 1:(length(s)-1)) {
    # run the bb intersection test, report to console
    cross <- doBoxesIntersect(
      getBoundingBox(P0 = c(xy[n,1], xy[n,2]), P1 = c(xy[n+1,1], xy[n+1,2])),
      getBoundingBox(P0 = c(xy[n+2,1], xy[n+2,2]), P1 = c(xy[n+3,1], xy[n+3,2])))

    cat("Bounding box", n, "intersects bounding box", n+1, ":", cross, "\n")

    # run the final segment intersection test, report to console
    inter <- doSegmentsIntersect(
      segment1 = c(xy[n,1], xy[n,2], xy[n+1,1], xy[n+1,2]),
      segment2 = c(xy[n+2,1], xy[n+2,2], xy[n+3,1], xy[n+3,2]))

    cat("Segment", n, "intersects segment", n+1, ":", inter, "\n")
  }

} # end of intersectTest

