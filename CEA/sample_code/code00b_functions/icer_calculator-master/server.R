library(shiny)
source("icer_calculator.R")

shinyServer(function(input, output) {

  data = reactive({
    if (is.null(input$CE_data)){
      filename="data/example_data.csv"
    } else {
      filename=input$CE_data$datapath
    }
    print(filename)
    read.csv(filename, stringsAsFactors=FALSE)
  })
    
  icer_table = reactive({
    calculate_icers(data())
  })

  
  output$icer_table = renderDataTable({
    withProgress(message = 'Calculating ICERs',{
      datatable(icer_table(),
                style = 'bootstrap',
                rownames = FALSE,
                colnames = colnames(icer_table()),
                options = list(pageLength = 25, autoWidth = TRUE, dom='ftrpi'))
    })
  })
  
  output$icer_plot = renderPlot({
    withProgress(message = 'Rendering Cost-effectiveness plane',{
      plot_strategies(icer_table(),input$show_labels)
    })
  })
  
})