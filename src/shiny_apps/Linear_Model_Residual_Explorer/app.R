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
library(tibble)
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
   titlePanel("Model Residuals"),

   # Sidebar with a slider input for number of bins
   fluidRow(
     column(4,
            selectInput("gene",
                        "Gene of Interest:",
                        choices = sort(unique(
                          as.character(genes$name)
                        )))),
     plotOutput("genePlot")

   )
)

make_plot = function(gene, all) {
  dat = all %>% dplyr::filter(name == gene) %>% unnest()

  p_heatmap1 = ggplot(dat, aes(x = Strain, y = Condition, fill = Base_model_aug_resid)) +
                geom_tile() +
                scale_fill_gradient2(low = ('cyan'), mid = 'black', high = 'yellow') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(fill = "Residual") +
                ggtitle("Expression ~ 1")

  p_heatmap2 = ggplot(dat, aes(x = Strain, y = Condition, fill = Strain_model_aug_resid)) +
                geom_tile() +
                scale_fill_gradient2(low = ('cyan'), mid = 'black', high = 'yellow') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(fill = "Residual") +
                ggtitle("Expression ~ Strain")

  p_heatmap3 = ggplot(dat, aes(x = Strain, y = Condition, fill = Condition_model_aug_resid)) +
                geom_tile() +
                scale_fill_gradient2(low = ('cyan'), mid = 'black', high = 'yellow') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(fill = "Residual") +
                ggtitle("Expression ~ Condition")

  p_heatmap4 = ggplot(dat, aes(x = Strain, y = Condition, fill = Full_model_aug_resid)) +
                geom_tile() +
                scale_fill_gradient2(low = ('cyan'), mid = 'black', high = 'yellow') +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                labs(fill = "Residual") +
                ggtitle("Expression ~ Strain + Condtion")



  grid.arrange(p_heatmap1, #+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()),
               p_heatmap2, #+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()),
               p_heatmap3, #+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()),
               p_heatmap4) #+ theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()))
}


# Define server logic required to draw a histogram
server <- function(input, output) {

   output$genePlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      make_plot(input$gene, genes)
   })
}

# Run the application
shinyApp(ui = ui, server = server)
