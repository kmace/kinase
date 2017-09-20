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
source('plotting.R')
addResourcePath('WGCNA_modules', './WGCNA_modules')
# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Cluster Analysis"),

  # Select module, show module stats, and show close modules
  fluidRow(
    column(width = 4,
           numericInput(inputId = "module",
                       label = "Select a Module:", min = 0,
                       max = max(unique(as.numeric(modules$module))),
                       value = 0
                      ),
           selectInput(inputId = "gene",
                       label = "Use a gene to find a Module:",
                       choices = sort(unique(modules$name))
           )
          ),
     column(width = 4,
            verbatimTextOutput("moduleSummary")
           ),
     column(width = 4,
            plotOutput("moduleSimilarity")
           )
  ),
  # Full Expression Matrix
  titlePanel("Expression"),
  fluidRow(
    plotOutput("moduleExpression")
  ),
  # Residual Matrix
  titlePanel("Residual"),
  fluidRow(
    plotOutput("moduleResidual")
  ),
  # Go Enrichment
  titlePanel("GO Enrichment"),
  fluidRow(
    dataTableOutput("GOEnrichmentTable")
  ),
  # Show module Genes
  titlePanel("Modual Membership"),
  fluidRow(
    dataTableOutput("membershipTable")
  ),
  # Show Meme output
  titlePanel("MEME Results"),
  fluidRow(
    uiOutput("memePage")
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, clientData, session) {

  observe({
    geneMod = modules %>% filter(name == input$gene) %>% pull(module)

  updateNumericInput(session, "module",
                    value = geneMod)
  })

  output$moduleSummary = renderText(get_moduleSummary(input$module, modules, tfs))

  output$moduleSimilarity = renderPlot(get_moduleSimilarity(input$module, MEsCor))

  output$moduleExpression = renderPlot(get_moduleExpression(input$module, module_raw_residuals))

  output$moduleResidual = renderPlot(get_moduleResidual(input$module, module_std_residuals))

  output$GOEnrichmentTable = renderDataTable(get_GOEnrichmentTable(input$module, tab),
                                             options = list(pageLength = 10))

  output$membershipTable = renderDataTable(get_membershipTable(input$module, modules),
                                            options = list(pageLength = 10))

  #output$memePage = renderImage(get_memePage(input$module))
  output$memePage = renderUI(get_memePage(input$module))
  #output$memePage = renderUI()
}

# Run the application
shinyApp(ui = ui, server = server)

