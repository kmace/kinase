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
load('residuals.RData')
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

make_plot = function(gene, pathway, all, t2g) {
  dat = all %>% dplyr::filter(name == gene)
  p_heatmap1 = ggplot(dat, aes(y=Condition, x = Strain)) +
    geom_tile(aes(fill=Expression)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p_heatmap2 = ggplot(dat, aes(x = Strain, y = Condition, fill = Strain_residual)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


  p_heatmap3 = ggplot(dat, aes(x = Strain, y = Condition, fill = Condition_residual)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


  p_heatmap4 = ggplot(dat, aes(x = Strain, y = Condition, fill = Full_residual)) + geom_tile() + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))



  grid.arrange(p_heatmap1, #+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()),
               p_heatmap2, #+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()),
               p_heatmap3, #+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()),
               p_heatmap4) #+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()))
}


# Define server logic required to draw a histogram
server <- function(input, output) {

   output$genePlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      make_plot(input$gene, t2g$target_id[which(t2g$name == input$gene)], measurements, t2g)
   })
}

# Run the application
shinyApp(ui = ui, server = server)
