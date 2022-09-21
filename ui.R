# UI SET-UP ---------------------------------------------------------------

## this contains the app presentation
library(shinythemes)
library(DT)

# This apparently required to make the column-wise checkboxes possible.
tweaks <-
    list(tags$head(tags$style(HTML("
                                 .multicol {
                                   height: 250px;
                                   -webkit-column-count: 1; /* Chrome, Safari, Opera */
                                   -moz-column-count: 1;    /* Firefox */
                                   column-count: 1;
                                   -moz-column-fill: auto;
                                   -column-fill: auto;
                                 }"))
    ))


# UI DESIGN ---------------------------------------------------------------

## Set-up main-page
profile_page <- tabPanel("Profiles",
                   sidebarLayout(
                       sidebarPanel(
                           uiOutput("ui_select_sample"),
                           uiOutput("ui_select_scale"),
                           uiOutput("ui_select_range"),
                           #uiOutput("ui_select_chrom"),
                           tweaks,
                           fluidRow(column(width = 5, uiOutput("ui_select_chrom")))
                       ),

                       mainPanel(
                           tabsetPanel(type = "tabs",
                                       #tabPanel("Test", plotOutput("plot")),
                                       tabPanel("Plot", plotOutput("cnvplot_rect")),
                                       tabPanel("Table", DTOutput("cnvtable")))
                       )
                   ))


## Placeholer for more functionality or information
dummy_page <- tabPanel("ToDo",
                   h1("TODO"))



## INIT the UI
ui <- navbarPage("CNVision",
                 theme = shinytheme("darkly"),
                 profile_page,
                 dummy_page
)


### DEVOPS ###



















