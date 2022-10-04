# SERVER SET-UP -----------------------------------------------------------

library(quantmod)
library(ggplot2)
library(dplyr)
library(purrr)
library(shinyWidgets)
library(shinyFiles)

source("functions.R")


# LOAD INPUT DATA ---------------------------------------------------------
## load hg19/hg38 chromosome information
load(file="data/ref_genomes/ref_aesthetics.rda")

app.state <<- list(
    ref_select = "hg38",
    chromosomes = ref_aes[["hg38"]]$chrom,
    cselect = NULL,
    crange_min = 0,
    crange_max = Inf,
    datapath = paste(getwd(), "/data/cnv_profiles", sep=""),
    samples = NULL,

    scales = c("ratio", "scaled.ratio", "log2.ratio", "log10.ratio", "quantile.ratio")
)







# ref.select <- "hg38" # for the moment stick with hardcoding this
# ## INIT datapath variable with default path
# datapath <<- paste(getwd(), "/data/cnv_profiles", sep="")
# ## Prepare vector of sample-IDs for selection menu
# samples <- NULL  # names(data_storage)
# ## Prepare vector for chromosome selection
# chromosomes <<- ref_aes[[eval(ref.select)]]$chrom
# ## Prepare vector for choice of transformations
# scales <- c("ratio", "scaled.ratio", "log2.ratio", "log10.ratio", "quantile.ratio")


# SERVER LOGIC ------------------------------------------------------------

server <- function(input, output, session) {

    ## check if there is a data-file present
    if( length(list.files(path = paste(getwd(), "/data/saves", sep=""), pattern = "save")) == 0) {
        data_storage <<- NULL
    } else {
        # load existing data and init the ui-list
        load(file = paste(getwd(), "/data/saves/save_data.rda", sep=""))

        app.state$samples <<- names(data_storage)
        output$ui_select_sample <- renderUI(
            checkboxGroupInput("select_sample",
                               label="Sample selection:",
                               app.state$samples,
                               selected=app.state$samples[1])
        )
    }

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
        app.state$datapath
    })

    # This sets the datapath variable and performs safety check for variable access
    observeEvent(ignoreNULL = TRUE,
                 eventExpr = {
                     input$directory
                 },
                 handlerExpr = {
                     req(is.list(input$directory))
                     home <- normalizePath("~")
                     app.state$datapath <<-
                         file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))

                     ## Just updates the debug output ##
                     output$txt_file <- renderText({
                         #parseDirPath(c(home = '~'), dir())
                         app.state$datapath
                     })
                     ## Just updates the debug output ##
                 })

    # This monitors the button that starts the file import
    observeEvent(ignoreNULL = TRUE,
                 eventExpr = {
                     input$do
                 }, handlerExpr =  {
                     ## This part handles the file-import
                     data_storage <<- importData(app.state$datapath,
                                                 file.tag=".data.txt",
                                                 ratio.column="seg.mean.LOWESS")
                     app.state$samples <<- names(data_storage)
                     ## This part updates the check boxes with available samples
                     output$ui_select_sample <- renderUI(
                         checkboxGroupInput("select_sample",
                                            label="Sample selection:",
                                            app.state$samples,
                                            selected=app.state$samples[1])
                     )

                 })




    ### DIRECTORY SELECTION ###


    ## INPUTS
    # UI button for saving files
    # This will save the imported data to file for quicker access
    observeEvent(ignoreNULL = FALSE,
                 eventExpr = {
                     input$save
                 }, handlerExpr =  {

                     if( !is.null(data_storage) ) {
                         save(data_storage, file = paste(getwd(), "/data/saves/save_data.rda", sep=""))
                     }

                 })

    # UI for chromosome selection
    output$ui_select_chrom <- renderUI(
        list(strong(p("Chromosome selection:")),
             tags$div(align = "left",
                      class = "multicol",
                      checkboxGroupInput("select_chrom",
                                         label    = NULL,
                                         choices  = app.state$chromosomes,
                                         selected = NULL,
                                         inline   = FALSE)
             ))
    )

    # # OBSERVER for chromosome selection
    # observeEvent(ignoreNULL = FALSE,
    #              eventExpr = {
    #                  input$select_chrom
    #              }, handlerExpr =  {
    #                 print("I did some work")
    #                  if(is.null(input$select_chrom) ) {
    #                      print("if")
    #                      app.state$cselect <<- app.state$chromosomes
    #                  } else {
    #                      print("else")
    #                      app.state$cselect <<- input$select_chrom
    #                  }
    #
    #
    #              })

    # UI for ratio scale selection
    output$ui_select_scale <- renderUI(
        selectInput("select_scale",
                    label="Ratio scale:",
                    app.state$scales,
                    selected=app.state$scales[1])
    )

    # UI for chromosome range selection (dynamic)
    output$ui_select_range <- renderUI({

        # Range can only be selected if no more than a single chromosome is used
        if( length(input$select_chrom) == 1 && !is.null(input$select_chrom) ) {

            max.range <- ref_aes[[eval(app.state$ref_select)]] %>%
                dplyr::filter(., chrom %in% input$select_chrom) %>%
                dplyr::select(., clength) %>% unlist()

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

        data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% app.state$cselect) %>%
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

        if(is.null(input$select_chrom) ) {
            app.state$cselect <<- app.state$chromosomes
            } else {
                app.state$cselect <<- input$select_chrom
            }

        ## make sure that selected range works
        cnv_data <- data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% app.state$cselect)

        if( length(app.state$cselect) == 1 ) {
            cnv_data <- filterRange(cnv_data, input$dynamic[1], input$dynamic[2])
            plotCNV(cnv_data, ref_aes, p.type = 'rect', v.type = eval(input$select_scale),
                    c.min = input$dynamic[1], c.max = input$dynamic[2])
        } else {
            plotCNV(cnv_data, ref_aes, p.type = 'rect', v.type = eval(input$select_scale))
        }
    })


    output$cnvplot_seg <- renderPlot({

        if(is.null(input$select_chrom) ) {
            app.state$cselect <<- app.state$chromosomes
        } else {
            app.state$cselect <<- input$select_chrom
        }

        ## make sure that selected range works
        cnv_data <- data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% app.state$cselect)

        if( length(app.state$cselect) == 1 ) {
            cnv_data <- filterRange(cnv_data, input$dynamic[1], input$dynamic[2])
            plotCNV(cnv_data, ref_aes, p.type = 'seg', v.type = eval(input$select_scale),
                    c.min = input$dynamic[1], c.max = input$dynamic[2])
        } else {
            plotCNV(cnv_data, ref_aes, p.type = 'seg', v.type = eval(input$select_scale))
        }
    })

    output$cnvplot_point <- renderPlot({

        if(is.null(input$select_chrom) ) {
            app.state$cselect <<- app.state$chromosomes
        } else {
            app.state$cselect <<- input$select_chrom
        }

        ## make sure that selected range works
        cnv_data <- data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% app.state$cselect)

        if( length(app.state$cselect) == 1 ) {
            cnv_data <- filterRange(cnv_data, input$dynamic[1], input$dynamic[2])
            plotCNV(cnv_data, ref_aes, p.type = 'point', v.type = eval(input$select_scale),
                    c.min = input$dynamic[1], c.max = input$dynamic[2])
        } else {
            plotCNV(cnv_data, ref_aes, p.type = 'point', v.type = eval(input$select_scale))
        }

    }, bg = "white")


    output$cnvplot_color_rect <- renderPlot({

        if(is.null(input$select_chrom) ) {
            app.state$cselect <<- app.state$chromosomes
        } else {
            app.state$cselect <<- input$select_chrom
        }

        ## make sure that selected range works
        cnv_data <- data_storage[eval(input$select_sample)] %>%
            purrr::reduce(merge, all=TRUE) %>%
            filter(chrom %in% app.state$cselect)

        if( length(app.state$cselect) == 1 ) {
            cnv_data <- filterRange(cnv_data, input$dynamic[1], input$dynamic[2])
            plotCNV(cnv_data, ref_aes, p.type = 'rect_color', v.type = eval(input$select_scale),
                    c.min = input$dynamic[1], c.max = input$dynamic[2])
        } else {
            plotCNV(cnv_data, ref_aes, p.type = 'rect_color', v.type = eval(input$select_scale))
        }

    }, bg = "white")



    # sth old
    loadData <- function() {
        if (exists("responses")) {
            responses
        }
    }



}




