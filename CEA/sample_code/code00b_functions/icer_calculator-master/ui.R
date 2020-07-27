library(shiny)
library(DT)

fluidPage(
  titlePanel("ICER Calculator"),
  sidebarLayout(
    sidebarPanel(
      fileInput('CE_data', 'Choose CSV File',
                accept=c('text/csv', 
                         'text/comma-separated-values,
                         text/plain', 
                         '.csv')),
      
      checkboxInput('show_labels', label="Show strategy names on plot", value=FALSE),
      
      tags$div(
        HTML("<small>
             <p>This function takes input in CSV format containing the strategies to be compared. The CSV file should contain three columns named Strategy, Cost and Effect. The strategy column contains the name of the strategy and the other two columns are numerical estimates of the absolute cost and effect of the strategy.
             <p>A table of results will be produced describing whether or not each strategy is dominanted (D), extendedly dominated (ED) or non-dominated (ND). ICERs will be calculated for the non-dominated strategies.
             <p>A graph of the cost-effectiveness plane is also produced showing the different strategies plotted against each other.
             <p>This site was produced by <a href='https://www.york.ac.uk/che/staff/research/miqdad-asaria/'>Miqdad Asaria</a> 
             <p>Source code can be found <a href='https://github.com/miqdadasaria/icer_calculator'>here</a>.</small>")
        )
    ),
    mainPanel(
      tabsetPanel(id="tabset",
        tabPanel("ICERs", dataTableOutput("icer_table")),
        tabPanel("Cost-effectiveness Plane", plotOutput("icer_plot"))
      )
    )
  )
)