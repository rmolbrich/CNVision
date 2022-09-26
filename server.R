# SERVER SET-UP -----------------------------------------------------------

library(quantmod)
library(ggplot2)
library(dplyr)
library(purrr)
library(shinyWidgets)

source("functions.R")


# LOAD INPUT DATA ---------------------------------------------------------

## load hg38 chromosome
# load chromosome lengths to get the panel-widths right
chrom.sizes <- data.table::fread("data/ref_genomes/hg38.chrom.sizes")
# ToDo: selection of complete chromosomes could be more dynamic
chrom.sizes <- chrom.sizes[1:24,] %>%
    dplyr::rename_with(., ~ c("chrom", "clength")) %>%
    dplyr::mutate(., chrom = factor(chrom, levels = c(paste("chr",1:22,sep=""),"chrX","chrY"))) %>%
    dplyr::arrange(., chrom)

## load all the files in cnv_profiles folder
input_files <- list.files(path = "data/cnv_profiles/", pattern = ".data.txt")
## import the data
data_storage <- importData(input_files, chrom.sizes, "ratio")

## Prepare vector of sample-IDs for selection menu
samples <- names(data_storage)
## Prepare vector for chromosome selection
chromosomes <<- c("all", as.character(chrom.sizes$chrom))
## Prepare vector for choice of transformations
scales <- c("ratio", "scaled.ratio", "log2.ratio", "log10.ratio", "quantile.ratio")


# SERVER LOGIC ------------------------------------------------------------

server <- function(input, output) {

    #### Directory selection logic ####

    shinyDirChoose(input, 'folder', roots=c(wd='.'), filetypes=c('', 'txt'))

    observe({
        print(input$folder)
    })


    #### Directory selection logic ####

    ## INPUTS
    # UI for sample selection
    output$ui_select_sample <- renderUI(
        checkboxGroupInput("select_sample",
                           label="Sample selection:",
                           samples,
                           selected=samples[1])
    )

    # UI for chromosome selection
    output$ui_select_chrom <- renderUI(
        list(h3("Chromosome selection:"),
             tags$div(align = 'left',
                      class = 'multicol',
                      checkboxGroupInput(inputId  = 'select_chrom',
                                         label    = NULL,
                                         choices  = chromosomes,
                                         selected = chromosomes[1],
                                         inline   = FALSE)))
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

            max.range <- chrom.sizes$clength[which(chrom.sizes$chrom %in%
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



    output$cnvplot__color_rect <- renderPlot({

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


    #### DEVELOPMENT CORNER ####

    ## this bit fixes the issue.. whatever the issue is?!
    # observe({
    #     if (input$static == "A") {
    #         values$dyn <- input$dynamic
    #     } else {
    #         values$dyn <- NULL
    #     }
    # })

    # The example output
    output$check <- renderText({

        paste(input$static, values$dyn)

    })

    ## list to store reactive values
    values <- reactiveValues()
    #### DEVELOPMENT CORNER ####


}




