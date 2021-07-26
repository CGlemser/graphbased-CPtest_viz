library(shiny)
source("fct_createData.R")
source("gseg_Funs.R")

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$t_UI <- renderUI({
      sliderInput("t", "hypothesized change point t",
                  min = 2, max = input$n-1,
                  value = floor(input$n/2), step = 1)
  })
  
  output$tau_UI <- renderUI({
    sliderInput("tau",
                "relative position of true change point tau",
                min = 2, max = input$n-1,
                value = floor(input$n/2), step = 1)
  })
    
  storage <- reactiveValues(dat = data.table(), graphs = list(),
                            stats = list(), showTestTables = TRUE)
  
  observeEvent(input$sim_data, {
    storage$dat <- createDataSet(n = input$n, delta = input$delta,
                                       delta_var = input$delta_var,
                                       tau = input$tau, d = input$dim,
                                       distribution = "normal")
    MST <- createSimilarityGraph(storage$dat, "MST")
    MDP <- createSimilarityGraph(storage$dat, "MDP")
    NNG <- createSimilarityGraph(storage$dat, "NNG")
    storage$graphs <- list(MST = MST, MDP = MDP, NNG = NNG)
  })
  

  
  output$MST <- renderPlot({
    # since input$t is a renderUI it is not immediately available!
    if(is.null(input$t) | is.null(storage$graphs$MST)){
      NULL
    } else {
      MST_E <- find_Redges(extractEdges(storage$graphs$MST), input$t)
      edge_highlights <- switch(input$R_highlight,
        None = NULL,
        R0 = MST_E[, "J0"],
        R1 = MST_E[, "J1"],
        R2 = MST_E[, "J2"]
      )
      
      plotSimilarityGraph_app(graph = storage$graphs$MST, dt = storage$dat, 
                              edge_highlights = edge_highlights) +
          ggtitle("MST") + theme(plot.title = element_text(hjust = 0.5))
    }
  })
  
  # render R data table
  output$MST_stats <- renderTable({
    if(is.null(input$t) | is.null(storage$graphs$MST)){
      NULL
    } else {
      calcExpectations(storage$graphs$MST, input$n, input$t)
    }
  }, rownames = TRUE)
  
  output$MDP <- renderPlot({
    # since input$t is a renderUI it is not immediately available!
    if(is.null(input$t) | is.null(storage$graphs$MDP)){
      NULL
    } else {
      MDP_E <- find_Redges(extractEdges(storage$graphs$MDP), input$t)
    
      edge_highlights <- switch(input$R_highlight,
                                None = NULL,
                                R0 = MDP_E[, "J0"],
                                R1 = MDP_E[, "J1"],
                                R2 = MDP_E[, "J2"]
      )
      
      plotSimilarityGraph_app(storage$dat, graph = storage$graphs$MDP,
                              edge_highlights = edge_highlights) +
          ggtitle("MDP") + theme(plot.title = element_text(hjust = 0.5))
    }
  })
  
  # render R data table
  output$MDP_stats <- renderTable({
    if(is.null(input$t) | is.null(storage$graphs$MDP)){
      NULL
    } else {
      calcExpectations(storage$graphs$MDP, input$n, input$t)
    }
  }, rownames = TRUE)
  
  output$NNG <- renderPlot({
    # since input$t is a renderUI it is not immediately available!
    if(is.null(input$t) | is.null(storage$graphs$NNG)){
      NULL
    } else {
      NNG_E <- find_Redges(extractEdges(storage$graphs$NNG), input$t)

      # boolean vector to indicate which edges contribute to R0, R1, R2
      edge_highlights <- switch(input$R_highlight,
                                None = NULL,
                                R0 = NNG_E[, "J0"],
                                R1 = NNG_E[, "J1"],
                                R2 = NNG_E[, "J2"]
      )
      
      plotSimilarityGraph_app(storage$dat, graph = storage$graphs$NNG,
                              edge_highlights = edge_highlights) +
          ggtitle("NNG") + theme(plot.title = element_text(hjust = 0.5))
    }
  })
  
  # render R data table
  output$NNG_stats <- renderTable({
    if(is.null(input$t) | is.null(storage$graphs$NNG)){
      NULL
    } else {
    calcExpectations(storage$graphs$NNG, input$n, input$t)
    }
  }, rownames = TRUE)
  
  
  ##== test functions ==##
  ##== MST ==##
  testRes_MST <- eventReactive(input$test_MST | input$sim_data, {
    withProgress(message = "Running permutations...", value = 0, {
      runGsegTest(storage$dat, storage$graphs$MST)
    })
  })
  
  output$testRes_MST <- renderTable(
    testRes_MST(), rownames = TRUE
  )
  
  ##== MDP ==##
  testRes_MDP <- eventReactive(input$test_MDP | input$sim_data, {
    withProgress(message = "Running permutations...", value = 0, {
      runGsegTest(storage$dat, storage$graphs$MDP)
    })
  })
  
  output$testRes_MDP <- renderTable(
    testRes_MDP(), rownames = TRUE
  )
  
  ##== NNG ==##
  testRes_NNG <- eventReactive(input$test_NNG | input$sim_data, {
    withProgress(message = "Running permutations...", value = 0, {
      runGsegTest(storage$dat, storage$graphs$NNG)
    })
  })
  
  output$testRes_NNG <- renderTable(
    testRes_NNG(), rownames = TRUE
  )
}
