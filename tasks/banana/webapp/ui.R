library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("Try to put the most variability on the x axis!"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("x",
                  "rotation along x",
                  min = 0,
                  max = 2 * pi,
                  value = 0),
      sliderInput("y",
                  "rotation along y",
                  min = 0,
                  max = 2 * pi,
                  value = 0),
      sliderInput("z",
                  "rotation along z",
                  min = 0,
                  max = 2 * pi,
                  value = 0),
      checkboxInput("incZero", 
                    label = "Include zero.",
                    value = FALSE)
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot"),
      textOutput("text1")
    )
  )
))
