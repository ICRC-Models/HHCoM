# ICER Calculator
Simple function to calculate incremental cost-effectiveness ratios (ICERs)

This function takes input in CSV format containing the strategies to be compared. The CSV file should contain three columns named Strategy, Cost and Effect. The strategy column contains the name of the strategy and the other two columns are numerical estimates of the absolute cost and effect of the strategy.

A table of results will be produced describing whether or not each strategy is dominanted (D), extendedly dominated (ED) or non-dominated (ND). ICERs will be calculated for the non-dominated strategies.

A graph of the cost-effectiveness plane is also produced showing the different strategies plotted against each other.

To run the shiny app from R you must have the shiny package installed and run the following command from the R console: shiny::runGitHub('icer_calculator', 'miqdadasaria')

Alternatively a running example can be found at the following address: https://miqdad.freeasinspeech.org.uk/icer_calculator/
