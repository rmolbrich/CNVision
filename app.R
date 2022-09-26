# APP SET-UP --------------------------------------------------------------

library(shiny)


source("ui.R")
source("server.R")



# APP EXECUTION -----------------------------------------------------------

shinyApp(ui = ui, server = server)




# cnv_data <- data_storage$IEG101
#
# ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 1, circular = TRUE)}
#
# bla<-ma(cnv_data$log10.ratio)
#
# cnv_data <- cnv_data[which(chrom %in% "chr1"),]
#
# x <- 1:dim(cnv_data)[1]
# y <- cnv_data$log10.ratio
# lo <- loess(y~x, degree=2)
# plot(x,y)
# lines(predict(lo), col='red', lwd=2)
# lines(bla, col="blue",lwd=2)
#
# cnv_data <- cnv_data %>%
#     dplyr::mutate(., movAvg = unlist(aggregate(log10.ratio ~ chrom + sID, data=cnv_data, movingAverage)$log10.ratio))
# lossFrequency <- unlist(aggregate(log10.ratio ~ chrom + sID, data=cnv_data, movingAverage)$log10.ratio)
#
#
#
#
# cnv_data












#
