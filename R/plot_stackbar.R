#' @title stackbar_plot function.
#'
#' @param data input data, Seurat obj or data.frame can be use.
#' @param order_sample reorder the samples.
#' @param group whether use group to use as legend, default is FALSE.
#'
#' @return stackbar_p A ggplot2 object containing the plot.
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all

# stackbar_plot function
stackbar_plot <- function(data, order_sample, group){
  plot_data <- do.call(rbind, data)
  # change sample order
  if ('default' %in% order_sample){
    plot_data <- plot_data
  } else {
    plot_data$plot_sample <- factor(plot_data$plot_sample,levels=order_sample)
  }
  if ( group == TRUE) {
    type_n <- length(unique(plot_data$plot_sample))
    colours <- c4a('dynamic',type_n)
    fill_legend <- c('Group')
    position <- c('fill')
    y <- plot_data$nums
    fill_x <- plot_data$plot_cell_type
    fill_y <- plot_data$plot_sample
    angle_x <- 50
    hjust_x <- 1.0
    vjust_x <- 1.0
  } else {
    type_n <- length(unique(plot_data$plot_cell_type))
    colours <- c4a('dynamic',type_n)
    fill_legend <- c('Cell Type')
    position <- c('stack')
    y <- plot_data$ratio
    fill_x <- plot_data$plot_sample
    fill_y <- plot_data$plot_cell_type
    angle_x <- 0
    hjust_x <- 0.5
    vjust_x <- 0.5
  }
  stackbar_p <- ggplot(plot_data, aes(x = fill_x, y = y, fill = fill_y)) +
    geom_col(position = position,width = 0.5, color = 'white', size = 0.5) +
    labs(x = NULL, y = NULL, fill = fill_legend) +
    scale_y_continuous(expand=c(0,0.01),labels = scales::percent) +
    scale_fill_manual(values=rev(colours[1:type_n])) + # 配置颜色
    theme_bw() +
    theme(axis.text.x=element_text(angle = angle_x, hjust= hjust_x, vjust = vjust_x ,size = 11,colour = "black"),
          axis.text.y=element_text(angle=0, hjust=0, vjust=0.5,size = 10,colour = "black"),
          axis.ticks.y=element_line(),
          panel.grid=element_blank()) +
    theme(legend.key.height = unit(5,'pt'),legend.key.width = unit(8,'pt'))
  cat('Stack BarPlot Finish!\n')
  return(stackbar_p)
}

#' @title Cell type ration visualization in bar plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_sample reorder the samples.
#' @param order_type reorder the cell type.
#' @param group whether use group to use as legend, default is FALSE.
#'
#' @return stackbar A ggplot2 object containing the plot.
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
#' plot_stackbar(input, group_by = 'sample',cell_type = 'cluster.type')

# plot stackbar function
plot_stackbar <- function(input, group_by=NULL, cell_type=NULL,sub_sample='all',sub_type='all',order_sample='default', order_type = 'default', group = FALSE){
  prerared <- prepare_data(input, group_by, cell_type,sub_sample, sub_type, order_type)
  stackbar <- stackbar_plot(prerared,order_sample, group)
  return(stackbar)
}
