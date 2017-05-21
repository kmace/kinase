#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
load('shiny.RData')
source('kinase_helper.R')

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Kinase Analysis"),
  # Sidebar with a inputs 
  sidebarLayout(sidebarPanel(
    # Select Kinase and mode
    selectInput("kinase",
                "Kinase of Interest:",
                choices = sort(unique(
                  as.character(drug_meta$Strain)
                )[-7])),
    radioButtons(
      "coloring",
      "Coloring:",
      choices = c('Differential Expression', 'TF Binding')
    ),
    conditionalPanel(
      condition = "input.coloring == 'Differential Expression'",
      numericInput(
        'fc_cutoff',
        "abs log2 Fold Change cutoff:",
        value = 1,
        min = 0
      ),
      numericInput(
        'pvalue_cutoff',
        'Maximum P-value:',
        value = 0.05,
        min = 0,
        max = 1
      )
    ),
    conditionalPanel(
      condition = "input.coloring == 'TF Binding'",
      selectInput('tf_dataset',
                  'Select Dataset:',
                  choices = sort(
                    tidyr::unite(filter(tf_meta, expert == TRUE), all, name, expert, cname, sep = '::')$all
                  )),
      checkboxInput('na_rm',
                    'Remove Missing Values',
                    value = TRUE),
      checkboxInput('useCutoff',
                    'Use Cutoff',
                    value = FALSE),
      conditionalPanel(condition = "input.useCutoff == TRUE",
                       numericInput('tf_cutoff',
                                    "Probability Cutoff",
                                    value = 0))
    ), width = 2 ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("volcanoPlot")
  )
))

# Define server logic required to draw a histogram
server <- function(input, output) {
     output$volcanoPlot <- renderPlot({
       if(input$coloring == 'Differential Expression'){
       p = volcano_de(dds, 
                 input$kinase, 
                 t2g, 
                 input$pvalue_cutoff, 
                 input$fc_cutoff)
     } else if(input$coloring == 'TF Binding'){
       cutoff = NA
       if(input$useCutoff){
         cutoff = input$tf_cutoff
       }
       p = volcano_tf(dds, 
                 input$kinase, 
                 t2g, 
                 input$tf_dataset,
                 cutoff,
                 input$na_rm)
     }
       print(p)
       return(p)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

