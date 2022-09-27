# UI SET-UP ---------------------------------------------------------------

## this contains the app presentation
library(shinythemes)
library(shinyFiles)
library(DT)

# This apparently required to make the column-wise checkboxes possible.
tweaks <-
    list(tags$head(tags$style(HTML("
                                 .multicol {
                                 font-size:12px;
                                 -webkit-column-count: 3;
                                 -moz-column-count: 3;
                                 column-count: 3;
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
                           width = 3,
                           #shinyDirButton('folder', 'Select a folder', 'Please select a folder', FALSE),
                           uiOutput("ui_select_sample"),
                           uiOutput("ui_select_scale"),
                           uiOutput("ui_select_range"),

                           tweaks,
                           fluidRow(column(width = 12, uiOutput("ui_select_chrom"))),
                       ),

                       mainPanel(
                           tabsetPanel(type = "tabs",
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
                 theme = shinytheme("paper"),
                 profile_page,
                 dummy_page
)


### DEVOPS ###



















