library(shiny)
library(glmnet)
load('.RData')

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {

  # Expression that generates a plot of the distribution. The expression
  # is wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically
  #     re-executed when inputs change
  #  2) Its output type is a plot
  #
  output$estimationPlot <- renderPlot({
    # generate an rnorm distribution and plot it
    i = input$complexity
    plot(predict(fit1,newx=dictionary)[,i],
         target[,1],
         xlab='Prediciton',
         ylab="Actual",
         main=paste('Lambda: ', fit1$lambda[i], ' | DF: ', fit1$df[i], ' | % VarExp: ', fit1$dev.ratio[i]))
    
  })
})
