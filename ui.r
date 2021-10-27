if (!require("DT")) install.packages('DT')
library(shiny)
library(shinyjs)
library(DT)

shinyUI(fluidPage(
  #shinythemes::themeSelector(),
  title = "Outpatient Data",
  
  h1("Outpatient Data"),
  
  fluidRow(
    column(6, verbatimTextOutput("debug")),
    actionButton("page","Select Page"),
    actionButton("all","Select All"),
    actionButton("clear","Clear Selection"),
    downloadButton("download","Download")
  ),
  
  fluidRow(
    column(6, DT::dataTableOutput("masterTable")),
    
    column(6, mainPanel(
      
      tabsetPanel(type = "tabs",
                  tabPanel("Histogram",
                           useShinyjs(),
                           plotOutput("hist"),
                           sliderInput("binwidth",
                                       "Size of bins",
                                       min = 0, max = 5, value = 1),
                           radioButtons("histType",
                                        "Histogram Type",
                                        choices = list("Single" = "single",
                                                       "Side-by-side" = "sbs")),
                           radioButtons("sbsOptions",
                                        "Side-by-side Factor",
                                        choices = list("Sex" = "sex",
                                                       "Ethnicity" = "ethnicity")),
                           radioButtons("histOptions",
                                        "Histogram Options",
                                        choices = list("Slope" = "slope",
                                                       "eGFR" = "ckd_epi_egfr",
                                                       "Number of eGFR Measurements" = "number_egfrs",
                                                       "Number of Follow Up Years" = "followup_years",
                                                       "Sex (0 = Male, 1 = Female)" = "sex",
                                                       "Ethnicity (0 = White, 1 = Black)" = "ethnicity",
                                                       "First Age" = "first_age",
                                                       "First eGFR Measurement" = "first_egfr",
                                                       "First CKD Stage" = "first_ckd_stage",
                                                       "Last Age" = "last_age",
                                                       "Last eGFR Measurement" = "last_egfr",
                                                       "Last CKD Stage" = "last_ckd_stage",
                                                       "Renal ICD" = "renal_icd",
                                                       "Diabetes Type" = "dm_type")),
                           
                           verbatimTextOutput("histSum")),
                  tabPanel("Line",
                           plotOutput("line"),
                           radioButtons("lineOptions",
                                        "Units for the x-axis",
                                        choices = list("Years" = "years", "Days" = "days")),
                           radioButtons("labOptions",
                                        "Which lab to plot",
                                        choices = list("eGFR" = "egfr", "BMI" = "bmi", "Blood Pressure" = "bp")),
                           verbatimTextOutput("lineExplain")),
                  tabPanel("Regression",
                           plotOutput("regression"),
                           radioButtons("regressionOptions",
                                        "Units for the x-axis (not working)",
                                        choices = list("Years" = "years", "Days" = "days")),
                           radioButtons("includeEvent",
                                        "ICD Events (not working)",
                                        choices = list("Show Events" = "graphEvents","Hide Events" = "dontGraphEvents")),
                           verbatimTextOutput("regressionExplain")),
                  tabPanel("Table",
                           tableOutput("table"),
                           #radioButtons("tableOptions"),
                           verbatimTextOutput("tableExplain"))
      )
    ))
  )
  
  # h1("Debug"),
  # fluidRow(
  #   column(6, verbatimTextOutput("debug"))
  # ),
  # 
  # titlePanel("Outpatient Data"),
  # sidebarLayout(
  #     sidebarPanel(
  #       radioButtons("graph","Graph Type",c("Table" = "table","Line" = "line")),
  #       width = 50),
  #     mainPanel(dataTableOutput("masterTable")
  #     )
  # )
    # sidebarLayout(
    #   
    #   mainPanel("Outpatient Data",
    #             dataTableOutput("gapminder_table"))
    # )
))
