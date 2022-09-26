# UI SET-UP ---------------------------------------------------------------

## this contains the app presentation
library(shinythemes)
library(shinyFiles)
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
                         # titlePanel(
                         #     # app title/description
                         #     "muh"
                         # ),
                         # shinyUI(bootstrapPage(
                         #     shinyDirButton('folder', 'Select a folder', 'Please select a folder', FALSE)
                         # )),
                   sidebarLayout(
                       sidebarPanel(
                           shinyDirButton('folder', 'Select a folder', 'Please select a folder', FALSE),
                           uiOutput("ui_select_sample"),
                           uiOutput("ui_select_scale"),
                           uiOutput("ui_select_range"),
                           #uiOutput("ui_select_chrom"),
                           tweaks,
                           fluidRow(column(width = 5, uiOutput("ui_select_chrom"))),
                           width = 3
                       ),

                       mainPanel(
                           tabsetPanel(type = "tabs",
                                       #tabPanel("Test", plotOutput("plot")),
                                       tabPanel("Sample Plot",
                                                plotOutput("cnvplot__color_rect",
                                                           width = "100%",
                                                           height = "700px")),
                                       tabPanel("Segment Plot", plotOutput("cnvplot_seg",
                                                                               width = "100%",
                                                                               height = "700px")),
                                       tabPanel("Comparative Plot", plotOutput("cnvplot_rect",
                                                                   width = "100%",
                                                                   height = "700px")),
                                       tabPanel("Bin Table", DTOutput("cnvtable")))
                       ), fluid = FALSE
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



















