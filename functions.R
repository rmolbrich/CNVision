# CNV ANALYSIS FUNCTIONS --------------------------------------------------

#' Reduce cnv-table to set range
#'
#'@param cnv.table table object of cnv-data
#'@param min.range minimal range on chromosome
#'@param ax.range maximal range on chromosome
#'
#'@return cnv.table reduced to selected range
filterRange <- function(cnv.table, min.range, max.range) {
    return(filter(cnv.table, chrompos >= min.range,
                  chromendpos <= max.range))
}


#' Import cnv-profile tables
#'
#' For a list of cnv-profiles tables the function will import them into a list
#' named by the first filename segments separated by an underscore ('_').
#'
#' Will also perform the re-scaling and transformation operations for now.
#'
#' @param input_files filenames located in local data/cnv_profiles folder.
#' @param chrom.sizes table with chromosome lengths of chosen reference
#'
#' @return data_storage named list object containing cnv-tables
importData <- function(input_files, chrom.sizes) {
    data_storage <- NULL
    for( f.idx in input_files ) {

        data_storage[[unlist(strsplit(f.idx, "_", fixed=TRUE))[1]]] <- data.table::fread(paste("data/cnv_profiles/", f.idx, sep="")) %>%
            dplyr::mutate(., chrom = ifelse(chrom == 23, "chrX",
                                            ifelse(chrom == 24, "chrY",
                                                   paste("chr", chrom, sep="")))) %>%
            dplyr::mutate(., chrom = factor(chrom, levels = c(paste("chr",1:22,sep=""),
                                                              "chrX","chrY"))) %>%
            dplyr::arrange(., chrom) %>%
            # create additional column with ratios rescaled to I=[-1,1]
            dplyr::mutate(., scaled.ratio = setInterval(ratio,output.range = c(-1,1),input.midvalue = 1)) %>%
            dplyr::mutate(., log2.ratio = log2(ratio)) %>%
            dplyr::mutate(., log10.ratio = log10(ratio)) %>%
            dplyr::mutate(., sID = unlist(strsplit(f.idx, "_", fixed=TRUE))[1]) %>%
            dplyr::left_join(., chrom.sizes, by = "chrom")
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
plotCNV <- function(cnv_data, p.type = "rect", v.type = "ratio", c.min, c.max) {

    if(missing(c.min) || missing(c.max)) {
        c.min <- 0
        c.max <- cnv_data$clength
    }

    if( p.type == 'rect' ) {

        p.obj <- ggplot(cnv_data, aes(ymin = 0)) +
            geom_rect(aes(xmin = chrompos, xmax = chromendpos,
                          ymax = get(v.type), colour = sID, fill = sID,
                          group = sID, alpha = 0.25)) +
            ## plot chromosome lengths to get panel widths right
            geom_segment(aes(x = c.min, y = 0, xend = c.max, yend = 0),
                         color="black", linetype = "dashed", size = 0.2,
                         alpha = 0.5) +
            #facet_grid(cols=vars(chrom), scales='free_x', space='free_x') +
            facet_wrap(facets=vars(chrom), nrow=2, scales='free_x') +
            theme_cnv() +
            labs(x="HG38 5K-bins", y="Ratio of bin counts to genomic average",
                 name="Covariate")

    } else if( p.type == 'seg' ) {

        p.obj <- ggplot(cnv_data, aes(ymin = -1)) +
            ## plot the bin values
            geom_segment(aes(x = chrompos, y = get(v.type), xend = chromendpos,
                             yend = get(v.type), colour = "segment",
                             fill = sID, group = sID, alpha = 0.25)) +
            ## plot chromosome lengths to get panel widths right
            geom_segment(aes(x = c.min, y = 0, xend = c.max, yend = 0),
                         color="black", linetype = "dashed", size = 0.2,
                         alpha = 0.5) +
            ## split by chromosome to produce panels
            #facet_grid(cols=vars(chrom), scales='free_x', space='free_x') +
            facet_wrap(facets=vars(chrom), nrow=2, scales='free_x') +
            theme_cnv() +
            labs(x="HG38 5K-bins", y="Ratio of bin counts to genomic average",
                 name="Covariate")

    } else {
        # Fancy way of doing jack.
        p.obj <- NULL
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
                       axis.ticks.x=ggplot2::element_blank(),
                       legend.position=eval(legend_position),
                       panel.background = element_rect(fill='transparent'), #transparent panel bg
                       plot.background = element_rect(fill='transparent',
                                                      color=NA),
                       #panel.grid.major = element_blank(), #remove major gridlines
                       #panel.grid.minor = element_blank(), #remove minor gridlines
                       rect = element_rect(fill = "transparent"),
                       strip.text.x = ggplot2::element_text(
                           size = eval(base_size), color = "black", face = "bold.italic"))
}












