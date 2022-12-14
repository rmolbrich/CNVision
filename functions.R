# CNV ANALYSIS FUNCTIONS --------------------------------------------------

#' Reduce cnv-table to set range
#'
#'@param cnv.table table object of cnv-data
#'@param min.range minimal range on chromosome
#'@param max.range maximal range on chromosome
#'
#'@return cnv.table reduced to selected range
filterRange <- function(cnv.table, min.range, max.range) {
    return(filter(cnv.table, chrompos >= min.range,
                  chromendpos <= max.range))
}



#' Compute moving average of a value
#'
#'@param input.vec Measurements to average
#'@param window.size Filter size for trailing average
#'
#'@return Vector of averages over input, same size as input vector
movingAverage <- function(input.vec, window.size = 20) {

    avg.vec <- stats::filter(input.vec, rep(1 / window.size, window.size),
                             sides = 1, circular = TRUE)

    return(avg.vec)
}



#' Import cnv-profile tables
#'
#' For a list of cnv-profiles tables the function will import them into a list
#' named by the first filename segments separated by an underscore ('.').
#'
#' Will also perform the re-scaling and transformation operations for now.
#'
#' @param input_files filenames located in local data/cnv_profiles folder.
#' @param chrom.sizes table with chromosome lengths of chosen reference
#' @param file.tag Tag to filter input files by. Defaults to '.data.'
#' @param ratio.column For calculations, ... maybe use sth. else
#'
#' @return data_storage named list object containing cnv-tables
importData <- function(file.path, chrom.sizes, file.tag = ".data.", ratio.column = "ratio") {

    # - list all files
    # - split by sID, ref genome and bin size
    # - import 'data' and 'short' files separately
    # - format and calculate metrics

    input_files <- list.files(path = eval(file.path), pattern = eval(file.tag))

    data_storage <- NULL
    for( f.idx in input_files ) {

        data_storage[[unlist(strsplit(f.idx, ".", fixed=TRUE))[1]]] <- data.table::fread(paste("data/cnv_profiles/", f.idx, sep="")) %>%
            dplyr::mutate(., chrom = ifelse(chrom == 23, "chrX",
                                            ifelse(chrom == 24, "chrY",
                                                   paste("chr", chrom, sep="")))) %>%
            dplyr::mutate(., chrom = factor(chrom, levels = c(paste("chr",1:22,sep=""),
                                                              "chrX","chrY"))) %>%
            dplyr::arrange(., chrom) %>%
            # create additional column with ratios rescaled to I=[-1,1]
            #dplyr::mutate(., sID = unlist(strsplit(f.idx, "_", fixed=TRUE))[1]) %>%
            dplyr::mutate(., sID = unlist(strsplit(f.idx, ".", fixed=TRUE))[1]) %>%
            dplyr::mutate(., ref = unlist(strsplit(f.idx, "_", fixed=TRUE))[2]) %>%
            dplyr::mutate(., binSize = unlist(strsplit(f.idx, "_", fixed=TRUE))[3]) %>%
            dplyr::mutate(., scaled.ratio = setInterval(get(ratio.column),output.range = c(-1,1),input.midvalue = 1)) %>%
            dplyr::mutate(., log2.ratio = log2(get(ratio.column))) %>%
            dplyr::mutate(., log10.ratio = log10(get(ratio.column))) %>%
            dplyr::mutate(., color.ratio = factor(ifelse(get(ratio.column) < 1, "loss","gain"), levels = c("loss","gain"))) %>%
            dplyr::mutate(., movAvg = unlist(aggregate(log10.ratio ~ chrom + sID, data=., movingAverage)$log10.ratio)) %>%
            dplyr::mutate(., loss.ratio = ifelse(log10.ratio < 0, log10.ratio, 0)) %>%
            dplyr::mutate(., gain.ratio = ifelse(log10.ratio >= 0, log10.ratio, 0))
    }

    return(data_storage)

}


#' Scale to user-defined range
#'
#' Takes an input-vector an re-scales the values to fit in user-defined range.
#' Default target range is I=[-1,1].
#'
#' @keywords re-scaling, target range
#' @param input.vec A vector of measurements that should be rescaled
#' @param output.range A two element vector object that defines lower and upper
#' boundaries for the target interval. Range is I=[-1,1] by default.
#' @param input.midvalue use this to specify the middle of the input interval,
#' i.e., for an output interval of [-1,1] these are the values that are set to 0
#' @return A vector of same length as the input, re-scaled to target interval
setInterval <- function(input.vec, output.range = c(-1,1), input.midvalue = NULL) {
    i.min <- min(input.vec) # negative interval
    i.max <- max(input.vec) # positive interval

    if( is.null(input.midvalue) ) {
        output.vec <- ((input.vec - i.min) / (i.max - i.min)) *
            (output.range[2] - output.range[1]) + output.range[1]
    } else {
        # compute center of output range
        i.center <- output.range[2] - 0.5 * dist(output.range)
        # re-scale lower and upper half of input range
        output.vec <- sapply(input.vec, function(ival)
            ifelse(ival == input.midvalue, i.center,
                   ifelse(ival < input.midvalue,
                          ((ival - i.min) / (input.midvalue - i.min)) *
                              (i.center - output.range[1]) + output.range[1],
                          ((ival - input.midvalue) / (i.max - input.midvalue)) *
                              (output.range[2] - i.center) + i.center)))
    }
    return(output.vec)
}


#' Produce CNV profile plots from input tables
#'
#' Takes a formatted bin table an creates CNV-profile for given chromosomes
#' with user-defined plot-type and value range.
#'
#' @param cnv_data Output of SMURFSeq protocol but pre-processed (ToDo!)
#' @param p.type Select between Segment and Bar-plot with 'seg' or 'rect'
#' @param v.type Select value column in 'cnv_data' by name. Defaults to 'ratio'
#' @param c.min For chromosome range selection, this is the starting position.
#' @param c.max For chromosome range selection, this is the end-position.
plotCNV <- function(cnv_data, ref_aes, p.type = "rect", v.type = "log10.ratio", c.min, c.max) {
    ## color-palette
    copal <- pals::tableau20(20)[c(1,3,5,7,9,11,13,15,17,19)]

    ## filter aesthetics df
    frame.aes <- ref_aes[[eval(app.state$ref_select)]] %>%
        dplyr::filter(., chrom %in% unique(cnv_data$chrom))

    ## get min/max lengths of the plots right
    if(missing(c.min) || missing(c.max)) {
        c.min <- 0
        c.max <- frame.aes$clength
    }

    if( p.type == 'point' ) {

        p.obj <- ggplot() +
            geom_point(data = cnv_data,
                       aes(x = (chromendpos + chrompos)/2, y = get(v.type))) +
            ## Split by chromosome
            facet_wrap(facets=vars(chrom), nrow=6, ncol=4, scales='free_x') +
            ## Include center-line to get the widths right
            geom_segment(data = frame.aes, aes(x = c.min, y = 0, xend = c.max, yend = 0,
                                              group=chrom),
                         linetype = "dashed", size = 0.2,
                         alpha = 0.5, colour = "black") +
            ## draw rectangle for centromere dimensions
            geom_rect(data=frame.aes, aes(xmin = centroStart, xmax = centroEnd,
                                          ymin = min(cnv_data[,get(v.type)]),
                                          ymax = max(cnv_data[,get(v.type)])),
                      alpha = 0.2, fill = "blue") +
            ## Aesthetic shenanigance
            theme_cnv() +
            scale_color_manual(values=copal) +
            labs(x="HG38 5K-bins", y="Ratio of bin counts to genomic average",
                 name="Covariate")

    } else if( p.type == 'rect' ) {

        p.obj <- ggplot() +
            # geom_point(data = cnv_data,
            #            aes(x = (chromendpos + chrompos)/2, y = get(v.type))) +
            geom_rect(data = cnv_data,
                      aes(xmin = chrompos, xmax = chromendpos,
                          ymin = 0, ymax = get(v.type),
                          colour = sID, group = sID, alpha = 0.25), fill="white", size=0.1) +
            ## Split by chromosome
            facet_wrap(facets=vars(chrom), nrow=6, ncol=4, scales='free_x') +
            ## Include center-line to get the widths right
            geom_segment(data = frame.aes, aes(x = c.min, y = 0, xend = c.max, yend = 0,
                                               group=chrom),
                         linetype = "dashed", size = 0.2,
                         alpha = 0.5, colour = "black") +
            ## draw rectangle for centromere dimensions
            geom_rect(data=frame.aes, aes(xmin = centroStart, xmax = centroEnd,
                                          ymin = min(cnv_data[,get(v.type)]),
                                          ymax = max(cnv_data[,get(v.type)])),
                      alpha = 0.2, fill = "blue") +
            ## Aesthetic shenanigance
            theme_cnv() +
            scale_color_manual(values=copal) +
            labs(x="HG38 5K-bins", y="Ratio of bin counts to genomic average",
                 name="Covariate")


    } else if( p.type == 'seg' ) {

        p.obj <- ggplot() +

            geom_segment(data = cnv_data,
                         aes(x = chrompos, y = get(v.type), xend = chromendpos,
                             yend = get(v.type),
                             colour = sID, group = sID, alpha = 0.5),size=1) +
            ## Split by chromosome
            facet_wrap(facets=vars(chrom), nrow=6, ncol=4, scales='free_x') +
            ## Include center-line to get the widths right
            geom_segment(data = frame.aes, aes(x = c.min, y = 0, xend = c.max, yend = 0,
                                               group=chrom),
                         linetype = "dashed", size = 0.2,
                         alpha = 0.5, colour = "black") +
            ## draw rectangle for centromere dimensions
            geom_rect(data=frame.aes, aes(xmin = centroStart, xmax = centroEnd,
                                          ymin = min(cnv_data[,get(v.type)]),
                                          ymax = max(cnv_data[,get(v.type)])),
                      alpha = 0.2, fill = "blue") +
            ## Aesthetic shenanigance
            theme_cnv() +
            scale_color_manual(values=copal) +
            labs(x="HG38 5K-bins", y="Ratio of bin counts to genomic average",
                 name="Covariate")

    } else if( p.type == 'rect_color' ) {
        # Fancy way of doing jack.
        p.obj <- ggplot() +
            # geom_point(data = cnv_data,
            #            aes(x = (chromendpos + chrompos)/2, y = get(v.type))) +
            geom_rect(data = cnv_data,
                      aes(xmin = chrompos, xmax = chromendpos,
                          ymin = 0, ymax = get(v.type),
                          color=color.ratio, group = sID, alpha = 0.25), fill="white", size=0.1) +
            ## Split by chromosome
            facet_wrap(facets=vars(chrom), nrow=6, ncol=4, scales='free_x') +

            # draw moving average
            geom_path(data = cnv_data, aes(x = round((chromendpos + chrompos)/2), y = movAvg, group=1), color = "purple") +

            ## Include center-line to get the widths right
            geom_segment(data = frame.aes, aes(x = c.min, y = 0, xend = c.max, yend = 0,
                                               group=chrom),
                         linetype = "dashed", size = 0.2,
                         alpha = 0.5, colour = "black") +
            ## draw rectangle for centromere dimensions
            geom_rect(data=frame.aes, aes(xmin = centroStart, xmax = centroEnd,
                                          ymin = min(cnv_data[,get(v.type)]),
                                          ymax = max(cnv_data[,get(v.type)])),
                      alpha = 0.2, fill = "blue") +
            ## Aesthetic shenanigance
            theme_cnv() +
            scale_color_manual(values=copal) +
            labs(x="HG38 5K-bins", y="Ratio of bin counts to genomic average",
                 name="Covariate")
    }

    return(p.obj)

}


#' CNV theme for ggplot2 output.
#'
#' Make CNV plotting great again. Just a helper for 'plotCNV()' that takes care
#' of all the little adjustments required to produce pretty plots.
#'
#' @keywords theme CNV-profiles
#' @param legend_position Set legend position. Defaults to 'none'.
#' @param base_size Font size for better granularity.
#' @return A pretty theme for CNV-profiles.
#' @include mbecs_classes.R
theme_cnv <- function(legend_position = 'none', base_size = 14) {

    ggplot2::theme_bw() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       #axis.ticks.x=ggplot2::element_blank(),
                       legend.position=eval(legend_position),
                       #panel.background = element_rect(fill='transparent'), #transparent panel bg
                       #plot.background = element_rect(fill='transparent',
                        #                              color=NA),
                       #panel.grid.major = element_blank(), #remove major gridlines
                       #panel.grid.minor = element_blank(), #remove minor gridlines
                       rect = element_rect(fill = "transparent"),
                       strip.text.x = ggplot2::element_text(
                           size = eval(base_size), color = "black", face = "bold.italic"))
}












