#' Handy Function for select_aas for Quality Metrics.
#'
#' This function provides quality metrics plots for select_aas model with .
#'
#' @param aas.model <model of select_aas> model obtained from select_aas()
#' @param type <character> type of quality metric to be plotted.
#' example : *total_var*.
#' @param title <character> title for the plot.
#' example : *Archetypes Variability*
#' @param point_size <double> default 0.5; dimensions of the plotted point.
#' @param line_size <double> default 0.5; width of the plotted lines.
#' @param text_axis_size <integer or double> default 10; dimension of the text
#' on the axis.
#' @param plot_title_size <integer or double> default 15; dimension of
#' the title.
#' @param hjust <double> default 0.5; Position of the title (by default in the
#' central).
#' @return a Quality Metric plot defined by the type variable.
#' @examples
#' plot_aas_feature(aas.model, 'total_var', 'Archetypes Variability')
#' @export
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

    # enforcing the types
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
