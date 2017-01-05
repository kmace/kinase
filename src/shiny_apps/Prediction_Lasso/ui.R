library(shiny)

# Define UI for application that plots random distributions
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Expression Estimation by Lasso"),

  # Sidebar with a slider input for number of observations
  sidebarPanel(
    sliderInput("complexity",
                "Complexity:",
                min = 1,
                max = 99,
                value = 5)
  ),

  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("estimationPlot")
  )
))
