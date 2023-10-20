#' @title circle_plot function
#'
#' @param data input data, Seurat obj or data.frame can be use.
#' @param lab whether to show the label on the plot, default is TRUE.
#'
#' @return circle_p A ggplot2 object containing the plot.
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all
#'

# circle_plot function
circle_plot <- function(data, lab){
  circle_p <- list()
  for (sample in names(data)) {
    type_n <- length(unique(data[[sample]]$plot_cell_type))
    colours <- c4a('dynamic',type_n)
    sample_name <- unique(data[[sample]]$plot_sample)
    # can't change type order
    # if ('default' %in% order_type){
    #     data <- data
    # } else {
    #     data$order <- factor(data$order,levels=order_type)
    # }
    if (lab == T){
      geom_text <- geom_text(aes(label=nums,
                                 x=order , y=ratio), colour = "black",vjust = 0,hjust=1,size = 2)
    } else {
      geom_text <- NULL
    }
    circle_p[[sample]] <- ggplot(data[[sample]]) + geom_col(aes(x = order, y = ratio, fill = label)) +
      coord_polar("y") +
      theme_bw()+
      ggtitle(sample_name) +
      geom_text +
      labs(x = NULL, y = NULL, fill = 'Cell Types') +
      scale_fill_manual(values=colours) +
      theme(plot.title = element_text(hjust = 0.5, color = "black"),
            panel.grid=element_blank(),
            panel.border = element_blank(),
            axis.text.x=element_text(size = 10, hjust=1,colour = "black"),
            axis.text.y=element_blank(),axis.ticks.y=element_blank(),
            axis.ticks.x=element_line(),
            legend.key.height = unit(5,'pt'),legend.key.width = unit(8,'pt')) +
      scale_y_continuous(breaks = seq(from = 0,to = round(max(data[[sample]]$ratio), digits = 2)+0.02,by = 0.05),
                         labels = scales::percent,limits = c(0,round(max(data[[sample]]$ratio), digits = 2)+0.02))
  }
  cat('Circle Plot Finish!\n')
  return(circle_p)
}

#' @title Cell type ration visualization in circle plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_type reorder the cell type.
#' @param min_percent the minimum percent. use to combine some cell type with low ratio, default is 0.
#' @param lab whether to show the label on the plot, default is TRUE.
#'
#' @return circle A ggplot2 object containing the plot.
#'
#' @export
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all
#'
#' @examples
#' plot_circle(input, group_by = 'sample',cell_type = 'cluster.type')

plot_circle <- function(input, group_by=NULL, cell_type=NULL,sub_sample='all', sub_type='all', order_type = 'default',min_percent=0,lab=F){
  prerared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  data_new <- show_others(prerared,min_percent)
  circle <- circle_plot(data_new, lab)
  return(circle)
}





