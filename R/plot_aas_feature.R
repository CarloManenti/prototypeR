plot_aas_feature <- function(aas.model,
                             type,
                             title,
                             point_size=0.5,
                             line_size=0.5,
                             text_axis_size=10,
                             plot_title_size=15,
                             hjust=0.5){
    ### Description ###
    # Handy function to plot the quality metrics for archetypes features
    # across different ks (number of archetypes)

    # example usage
    # plot_aas_feature(aas.model, 'total_var', 'Archetypes Variability')


    # checking packages
    is_package_installed('ParetoTI')
    is_package_installed('ggplot2')


    # defining variables
    point_size      = as.double(point_size)
    line_size       = as.double(line_size)
    text_axis_size  = as.integer(text_axis_size)
    plot_title_size = as.integer(plot_title_size)
    hjust           = as.double(hjust)



    # Archetypes plot
    plot <-  ParetoTI::plot_arc_var(aas.model,
                                    type = type,
                                    point_size = point_size,
                                    line_size  = line_size) +
             ggplot2::theme_bw() +
             ggplot2::theme(axis.title = ggplot2::element_text(size = text_axis_size),
                            plot.title = ggplot2::element_text(size = plot_title_size,
                                                               hjust = hjust)) +
             ggplot2::ggtitle(title)

    return(plot)}
