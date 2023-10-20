#' @title sankey_plot function.
#'
#' @param data input data, Seurat obj or data.frame can be use.
#' @param order_sample reorder the samples.
#' @param coord_flip whether to flip the coord, default is FALSE.
#'
#' @return sankey_p A ggplot2 object containing the plot.
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all
#'

# sankey_plot function
sankey_plot <- function(data, order_sample, coord_flip){
  plot_data <- do.call(rbind, data)
  # change sample order
  if ('default' %in% order_sample){
    plot_data <- plot_data
  } else {
    plot_data$plot_sample <- factor(plot_data$plot_sample,levels=order_sample)
  }
  type_n <- length(unique(plot_data$plot_cell_type))
  colours <- c4a('dynamic',type_n)
  if (coord_flip == TRUE){
    flip <- coord_flip()
  } else {
    flip <- NULL
  }
  sankey_p <- ggplot(plot_data,
                   aes(x = plot_sample,
                       y = ratio,
                       fill = plot_cell_type,
                       stratum = plot_cell_type,
                       alluvium = plot_cell_type)) +
    geom_col(width = 0.6, color = 'white', size = 0.5) + # 柱状图边框
    geom_flow(width = 0.6, alpha = 0.4, knot.pos = 0.4, color = 'white', size = 0.5)+ # 柱间连接
    scale_fill_manual(values=rev(colours[1:type_n])) + # 配置颜色
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    labs(x = NULL,y = 'Cell Type Ratio',fill = 'Cell Type') +
    flip +
    theme(legend.key.height = unit(5,'pt'),
          legend.key.width = unit(8,'pt'),
          axis.text.x=element_text(colour = "black"),
          axis.text.y=element_text(colour = "black"))
  cat('Sankey Plot Finish!\n')
  return(sankey_p)
}


#' @title Cell type ration visualization in sankey plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_sample reorder the samples.
#' @param order_type reorder the cell type.
#' @param coord_flip whether to flip the coord, default is FALSE.
#'
#' @return sankey A ggplot2 object containing the plot.
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
#' plot_sankey(input, group_by = 'sample',cell_type = 'cluster.type')


# plot_sankey function
plot_sankey <- function(input, group_by=NULL, cell_type=NULL,sub_sample='all',sub_type='all',order_sample='default', order_type = 'default', coord_flip = F){
  prerared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  sankey <- sankey_plot(prerared, order_sample, coord_flip)
  return(sankey)
}

