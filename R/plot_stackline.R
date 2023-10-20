#' @title stackline_plot function
#'
#' @param data input data, Seurat obj or data.frame can be use.
#' @param order_sample reorder the samples.
#'
#' @return stackline_p A ggplot2 object containing the plot.
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all
#'

# stackline_plot function
stackline_plot <- function(data, order_sample){
  plot_data <- do.call(rbind, data)
  # change sample order
  if ('default' %in% order_sample){
    plot_data <- plot_data
  } else {
    plot_data$plot_sample <- factor(plot_data$plot_sample,levels=order_sample)
  }
  #
  if (length(as.character(unique(plot_data$plot_sample))) < 3) {
    type_n <- length(unique(plot_data$plot_sample))
    colours <- c4a('dynamic',type_n)
    fill_legend <- c('Group')
    fill_as <- c('plot_sample')
    x_as <- c('plot_cell_type')
  } else {
    type_n <- length(unique(plot_data$plot_cell_type))
    colours <- c4a('dynamic',type_n)
    fill_legend <- c('Cell Type')
    fill_as <- c('plot_cell_type')
    x_as <- c('plot_sample')
  }
  stackline_p <- ggplot(plot_data, mapping = aes(x = get(x_as), y = ratio, fill = get(fill_as), group= get(fill_as))) +
    geom_area(colour = 'white', size =0.3, alpha = 0.7) +
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size = 10, hjust=1,colour = "black"),
          axis.text.y=element_text(size = 10,colour = "black"),
          panel.grid=element_blank(),
          title=element_text(size=8),
          legend.key.height = unit(5,'pt'),legend.key.width = unit(8,'pt')) +
    labs(x = NULL, y = NULL, fill = fill_legend) +
    scale_fill_manual(values=rev(colours[1:type_n])) + # 配置颜色
    scale_y_continuous(labels = scales::percent)
  cat('Stack lineplot Finish!\n')
  return(stackline_p)
}

#' @title Cell type ration visualization in stack line plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_sample reorder the samples.
#' @param order_type reorder the cell type.
#'
#' @return stackline A ggplot2 object containing the plot.
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all
#'
#' @export
#'
#' @examples
#' plot_stackline(input, group_by = 'sample',cell_type = 'cluster.type')

# plot_stackline main function
plot_stackline <- function(input, group_by=NULL, cell_type=NULL,sub_sample='all',sub_type='all',order_sample='default', order_type = 'default'){
  prerared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  stackline <- stackline_plot(prerared, order_sample)
  return(stackline)
}
