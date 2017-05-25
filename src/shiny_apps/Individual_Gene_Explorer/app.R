#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RGraphics)
library(ggrepel)

#setwd('../../sci/')
#source('individual_gene.R')
#setwd('../shiny_apps/Individual_Gene_Explorer/')
load('shiny.RData')
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Single Gene Analysis"),
   
   # Sidebar with a slider input for number of bins 
   fluidRow(
     column(4,
            selectInput("gene",
                        "Gene of Interest:",
                        choices = sort(unique(
                          as.character(t2g$name)
                        )))),
     plotOutput("genePlot")
      
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$genePlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      make_plot(input$gene, t2g$target_id[which(t2g$name == input$gene)], all_data, t2g)
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

