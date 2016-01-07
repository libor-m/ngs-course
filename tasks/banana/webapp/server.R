library(shiny)
library(ggplot2)

# create a 3d rotation matrix
# https://www.math.duke.edu/education/ccp/materials/linalg/rotation/rotm3.html
rotX <- function(t) matrix(c(cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1), nrow=3)
rotY <- function(t) matrix(c(1, 0, 0, 0, cos(t), sin(t), 0, -sin(t), cos(t)), nrow=3)
rotZ <- function(t) matrix(c(cos(t), 0, -sin(t), 0, 1, 0, sin(t), 0, cos(t)), nrow=3)
rot3d <- function(tx, ty, tz) rotX(tx) %*% rotY(ty) %*% rotZ(tz)

# rotate a data frame with points in rows
rot3d_df <- function(df, tx, ty, tz) {
  # print(paste(tx, ty, tz))
  rmx <- rot3d(tx, ty, tz)
  res <- data.frame(t(apply(df, 1, function(x) rmx %*% as.numeric(x))))
  colnames(res) <- colnames(df)
  res
}

rb <- read.csv("data/rotated.csv")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  output$distPlot <- renderPlot({
    rot <- rot3d_df(rb, input$x, input$y, input$z)
    g <- ggplot(rot, aes(x, y)) + geom_point() + coord_equal()
    
    if(input$incZero)
      g <- g + ylim(-350, 350) + xlim(-350, 350)
    
    g
  })
  
  output$text1 <- renderText({
    class(input$x)
  })
})
