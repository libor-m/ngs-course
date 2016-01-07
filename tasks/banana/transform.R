# create a 3d rotation matrix
# https://www.math.duke.edu/education/ccp/materials/linalg/rotation/rotm3.html
rotX <- function(t) matrix(c(cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1), nrow=3)
rotY <- function(t) matrix(c(1, 0, 0, 0, cos(t), sin(t), 0, -sin(t), cos(t)), nrow=3)
rotZ <- function(t) matrix(c(cos(t), 0, -sin(t), 0, 1, 0, sin(t), 0, cos(t)), nrow=3)
rot3d <- function(tx, ty, tz) rotX(tx) %*% rotY(ty) %*% rotZ(tz)

# rotate a data frame with points in rows
rot3d_df <- function(df, tx, ty, tz) {
  rmx <- rot3d(tx, ty, tz)
  res <- data.frame(t(apply(df, 1, function(x) rmx %*% as.numeric(x))))
  colnames(res) <- colnames(df)
  res
}

# load and display the banana;)
d<-read.table("banana.tsv", col.names=c("x", "y","z"))
ggplot(d, aes(x, z, colour=x)) + 
  geom_point() + 
  coord_equal() + 
  scale_colour_gradientn(colours = rainbow(7), guide="none")

rb <- rot3d_df(d, .9, .1, 2)

ggplot(rb, aes(x, y)) + geom_point() + coord_equal()
ggplot(rb, aes(x, z)) + geom_point() + coord_equal()
ggplot(rb, aes(y, z)) + geom_point() + coord_equal()

# write the result when satisfied with the transformation
write.csv(rb, file="rotated.csv", row.names=F)

rb <- read.csv("rotated.csv")

# play with shiny
library(shiny)
runApp("webapp")

# PCA
pc <- prcomp(as.matrix(rb))
ggplot(data.frame(pc$x), aes(PC1, PC2)) + geom_point() + coord_equal()
ggplot(data.frame(pc$x), aes(PC1, PC3)) + geom_point() + coord_equal()
ggplot(data.frame(pc$x), aes(PC2, PC3)) + geom_point() + coord_equal()

# extract the roation, to set it up in web
# http://nghiaho.com/?page_id=846

xrot <- function(r) {
  tx <- atan2(r[3,2], r[3,3])
  ty <- atan2(-r[3,1], sqrt(r[3,2]^2 + r[3,3]^2))
  tz <- atan2(r[2,1], r[1,1])
  fixa <- function(a) ifelse(a < 0, 2 * pi + a, a)
  c(tx, ty, tz, fixa(tx), fixa(ty), fixa(tz))
}
# test the code
rmx <- rot3d(1, 1.1, 1.4)
xrot(t(rmx))

# result works in order x, z, y axes in shiny
xrot(pc$rotation)

# nmds
library(MASS)
dmx <- dist(rb)
mds <- isoMDS(dmx)
ggplot(data.frame(mds$points), aes(X1, X2)) + geom_point() + coord_equal()

# tsne
library(tsne)
rbx <- data.frame(rb, origx=d$x)
rbs <- sample_n(rbx, 500)
dmx <- dist(rbs[,1:3])
sne <- tsne(dmx)
ggplot(data.frame(sne, rbs), aes(X1, X2, colour=origx)) + 
  geom_point() + 
  coord_equal() + 
  scale_colour_gradientn(colours = rainbow(7), guide="none")

ggplot(rbs, aes(x, z, colour=origx)) + 
  geom_point() + 
  coord_equal() + 
  scale_colour_gradientn(colours = rainbow(7), guide="none")
