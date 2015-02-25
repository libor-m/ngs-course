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

# PCA
pc <- prcomp(as.matrix(rb))
ggplot(data.frame(pc$x), aes(PC1, PC2)) + geom_point() + coord_equal()
ggplot(data.frame(pc$x), aes(PC1, PC3)) + geom_point() + coord_equal()
# ggplot(data.frame(pc$x), aes(PC2, PC3)) + geom_point() + coord_equal()

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
