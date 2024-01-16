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
    pdf(file   = file_name,
        width  = width,
        height = height)
    # plotting (otherwise you will not save the plot,
    # but you will get a 'damaged' file)
    ComplexHeatmap::draw(plot)
    # closing the plot (thus saving the plot in pdf)
    dev.off()
}
