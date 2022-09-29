# SERVER SET-UP -----------------------------------------------------------

library(quantmod)
library(ggplot2)
library(dplyr)
library(purrr)
library(shinyWidgets)
library(shinyFiles)

source("functions.R")


# LOAD INPUT DATA ---------------------------------------------------------

## INIT datapath variable with default path
datapath <<- paste(getwd(), "/data/cnv_profiles", sep="")

## load hg19/hg38 chromosome information
load(file="data/ref_genomes/ref_data.rda")


## Prepare vector of sample-IDs for selection menu
samples <- NULL  # names(data_storage)
## Prepare vector for chromosome selection
chromosomes <<- c("all", as.character(ref_data[["hg38"]][["sizes"]]$chrom))
## Prepare vector for choice of transformations
scales <- c("ratio", "scaled.ratio", "log2.ratio", "log10.ratio", "quantile.ratio")


# SERVER LOGIC ------------------------------------------------------------

server <- function(input, output, session) {

    ### DIRECTORY SELECTION ###
    # get directory structure for executing server (local machine)
    volumes = getVolumes()

    # Set-up the UI interface logic
    shinyDirChoose(
        input,
        'directory',
        root = c(home = '~'),
        filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")

    )

    # make the directory input a reactive variable
    dir <- reactive(input$directory)

    # DEBUG: output current datapath to check results
    output$txt_file <- renderText({
        #parseDirPath(c(home = '~'), dir())
        datapath
    })

    # This sets the datapath variable and performs safety check for variable access
    observeEvent(ignoreNULL = TRUE,
                 eventExpr = {
                     input$directory
                 },
                 handlerExpr = {
                     req(is.list(input$directory))
                     home <- normalizePath("~")
                     datapath <<-
                         file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))

                     ## Just updates the debug output ##
                     output$txt_file <- renderText({
                         #parseDirPath(c(home = '~'), dir())
                         datapath
                     })
                     ## Just updates the debug output ##
                 })

    # This monitors the button that starts the file import
    observeEvent(ignoreNULL = TRUE,
                 eventExpr = {
                     input$do
                 }, handlerExpr =  {
                     data_storage <<- importData(datapath, ref_data[["hg38"]][["sizes"]], file.tag=".data.txt", ratio.column="seg.mean.LOWESS")

                     samples <<- names(data_storage)
                     # ## Prepare vector for chromosome selection
                     # chromosomes <<- c("all", as.character(chrom.sizes$chrom))
                     # ## Prepare vector for choice of transformations
                     # scales <<- c("ratio", "scaled.ratio", "log2.ratio", "log10.ratio", "quantile.ratio")
                 })




    ### DIRECTORY SELECTION ###
    loadData <- function() {
        if (exists("responses")) {
            responses
        }
    }

    ## INPUTS
    # UI for sample selection
    # This will initialize the list of samples on start and update it when the
    # import button is clicked
    observeEvent(ignoreNULL = FALSE,
                 eventExpr = {
                     input$do
                 }, handlerExpr =  {

                     output$ui_select_sample <- renderUI(
                         checkboxGroupInput("select_sample",
                                            label="Sample selection:",
                                            samples,
                                            selected=samples[1])
                     )

                 })

    # UI for chromosome selection
    output$ui_select_chrom <- renderUI(
        list(strong(p("Chromosome selection:")),
             tags$div(align = "left",
                      class = "multicol",
                      checkboxGroupInput("select_chrom",
                                         label    = NULL,
                                         choices  = chromosomes,
                                         selected = chromosomes[1],
                                         inline   = FALSE)
             ))
    )

    # UI for ratio scale selection
    output$ui_select_scale <- renderUI(
        selectInput("select_scale",
                    label="Ratio scale:",
                    scales,
                    selected=scales[1])
    )

    # UI for chromosome range selection (dynamic)
    output$ui_select_range <- renderUI({

        # Range can only be selected if no more than a single chromosome is used
        if( length(input$select_chrom) == 1 && input$select_chrom != "all" ) {

            max.range <- ref_data[["hg38"]][["sizes"]]$clength[which(ref_data[["hg38"]][["sizes"]]$chrom %in%
                                                       input$select_chrom)]
            # set max range to chrom length to skip further testing later on
            numericRangeInput("dynamic",
                              label="Chromosome range:",
                              value=c(0,max.range),
                              width = NULL,
                              min = 0,
                              max = max.range, # ToDo: use actually available chr range
                              step = 1)
        } else {
            return(NULL) # don't display the field
        }
    })

    ## OUTPUTS
    # UI for displaying CNV-table
    output$cnvtable <- renderDT({

        c.idx <- ifelse(input$select_chrom == "all", chromosomes[-1], input$select_chrom)

        data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% c.idx) %>%
            datatable(style = "bootstrap", options = list(pageLength = 200))

    })

    # UI for displaying rect-style plot
    # select more than one sample           - check
    # display in different colors           - check
    # select sample                         - check
    # select chromosome                     - check
    # select plot type
    # select scale / normalization range    - check
    # select chromosome range (single chr)  - check
    output$cnvplot_rect <- renderPlot({

        if( input$select_chrom %in% "all" ) {
            c.idx <- chromosomes[-1]
        } else {
            c.idx <- input$select_chrom
        }

        ## make sure that selected range works
        cnv_data <- data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% c.idx)

        if( length(c.idx) == 1 ) {
            cnv_data <- filterRange(cnv_data, input$dynamic[1], input$dynamic[2])
            plotCNV(cnv_data, p.type = 'rect', v.type = eval(input$select_scale),
                    c.min = input$dynamic[1], c.max = input$dynamic[2])
        } else {
            plotCNV(cnv_data, p.type = 'rect', v.type = eval(input$select_scale))
        }
    })


    output$cnvplot_seg <- renderPlot({

        if( input$select_chrom %in% "all" ) {
            c.idx <- chromosomes[-1]
        } else {
            c.idx <- input$select_chrom
        }

        ## make sure that selected range works
        cnv_data <- data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% c.idx)

        if( length(c.idx) == 1 ) {
            cnv_data <- filterRange(cnv_data, input$dynamic[1], input$dynamic[2])
            plotCNV(cnv_data, p.type = 'seg', v.type = eval(input$select_scale),
                    c.min = input$dynamic[1], c.max = input$dynamic[2])
        } else {
            plotCNV(cnv_data, p.type = 'seg', v.type = eval(input$select_scale))
        }
    })


    output$cnvplot_color_rect <- renderPlot({

        if( input$select_chrom %in% "all" ) {
            c.idx <- chromosomes[-1]
        } else {
            c.idx <- input$select_chrom
        }

        ## make sure that selected range works
        cnv_data <- data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% c.idx) %>%
            dplyr::mutate(., sID = factor(sID))

        if( length(c.idx) == 1 ) {
            cnv_data <- filterRange(cnv_data, input$dynamic[1], input$dynamic[2])
            plotCNV(cnv_data, p.type = 'rect_color', v.type = eval(input$select_scale),
                    c.min = input$dynamic[1], c.max = input$dynamic[2])
        } else {
            plotCNV(cnv_data, p.type = 'rect_color', v.type = eval(input$select_scale))
        }
    }, bg = "white")




    # The example output
    output$check <- renderText({

        paste(input$static, values$dyn)

    })

    ## list to store reactive values
    values <- reactiveValues()
    #### DEVELOPMENT CORNER ####


}




