#' Handy Function to Save a ComplexHeatmap in .pdf format
#'
#' This function saves a ComplexHeatmap in .pdf format
#'
#' @param plot ComplexHeatmap
#' @param file_name <character.pdf> name of the file used to store the plot.
#' @param width <double> default 5;
#' Width of the plot
#' @param height <double> default 3.5;
#' height of the plot
#' @return a .pdf file of a ComplexHeatmap plot
#' @examples
#' #save_pdf(plot = ComplexHeatmap::Heatmap(seq(1, 10, length.out = 10), col = grDevices::blues9), file_name = '~/Documents/plot.pdf')
#' @export
save_pdf <- function(plot,
                     file_name,
                     width=5,
                     height=3.25){

    ### Description ###
    # Saves a pdf file of your plot
    # Useful for ComplexHeatMaps plots

    # example usage
    # generate data
    # data <- matrix(rnorm(100, 1, 0.5), nrow = 20)
    # heatmap.plot <- ComplexHeatmap::Heatmap(data)
    # save_pdf(plot = heatmap.plot, file_name = '~/Documents/plot.pdf')

    is_package_installed('ComplexHeatmap')

    # initializing the pdf
    grDevices::pdf(file   = file_name,
                   width  = width,
                   height = height)
    # plotting (otherwise you will not save the plot,
    # but you will get a 'damaged' file)
    ComplexHeatmap::draw(plot)
    # closing the plot (thus saving the plot in pdf)
    grDevices::dev.off()
}
