library(shiny)
library(thematic)
library(bslib)
library(shinycssloaders)
library(shinyjs)

thematic::thematic_shiny(font = "auto")

fluidPage(
  tags$head(tags$style(".shiny-output-error{color: white;}")),
  shinyjs::useShinyjs(),
  theme = bs_theme(version = 4, bootswatch = "flatly"),
  
  titlePanel(h1("Graph-based Change-Point Test", align = "center"),
    windowTitle = "Graph-based Change-Point Test"),
  
  # Sidebar for data generation 
  sidebarLayout(
    sidebarPanel(width = 3,
      sliderInput("n",
                  "number of observations",
                  min = 10, max = 500,
                  value = 30, step = 10),
      numericInput("dim", "dimensionality",
                  min = 1, max = 1000,
                  value = 2),
      sliderInput("delta",
                  "locational change",
                  min = 0, max = 5,
                  value = 1, step = .5),
      sliderInput("delta_var",
                  "scale change",
                  min = 0, max = 5,
                  value = 0, step = .5),
      uiOutput("tau_UI"),
      actionButton("sim_data", "Simulate data!"),
      
      br(),
      br(),
      
      uiOutput("t_UI"),

      br(),
      
      radioButtons("R_highlight", "which R_i do you want highlighted?",
                   choices = c("None", "R0", "R1", "R2"),
                   selected = "None")
      
      # insert display R0, R1, R2
    ),
    
    # Show a plot of the
    mainPanel(width = 9,
      fluidRow(
        column(7, plotOutput("MST")),
        column(5, align = "center",
          br(), br(), br(),
          tableOutput("MST_stats"),
          br(),
          actionButton("test_MST", "Re-run the test!"),
          br(), br(),
          tableOutput("testRes_MST")
        )
      ),
      
      br(), br(),
      
      fluidRow(
        column(7, plotOutput("MDP")),
        column(5, align = "center",
          br(), br(), br(),
          tableOutput("MDP_stats"),
          actionButton("test_MDP", "Re-run the test!"),
          br(), br(),
          tableOutput("testRes_MDP")
        )
      ),
      
      br(), br(),
      
      fluidRow(
        column(7, plotOutput("NNG")),
        column(5, align = "center",
          br(), br(), br(),
          tableOutput("NNG_stats"),
          actionButton("test_NNG", "Re-run the test!"),
          br(), br(),
          tableOutput("testRes_NNG")
        )
      )
    )
  )
)