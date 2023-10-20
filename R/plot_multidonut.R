#' @title mutidonut_plot function.
#'
#' @param data input data, Seurat obj or data.frame can be use.
#' @param order_sample reorder the samples.
#' @param lab whether to show the label on the plot, default is TRUE.
#'
#' @return mutidonut_p A ggplot2 object containing the plot.

#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all

# mutidonut_plot
mutidonut_plot <- function(data, order_sample, lab){
  plot_data <- do.call(rbind, data)
  # change sample order
  if ('default' %in% order_sample){
    plot_data <- plot_data
  } else {
    plot_data$plot_sample <- factor(plot_data$plot_sample,levels=order_sample)
  }
  type_n <- length(unique(plot_data$plot_cell_type))
  colours <- c4a('dynamic',type_n)
  if (lab == T){
    geom_text <- geom_text(size = 1.5, position = position_stack(vjust = 0.5))
  } else {
    geom_text <- NULL
  }
  mutidonut_p <- ggplot(plot_data, aes(x = plot_sample, y = ratio, fill = plot_cell_type, label = nums)) +
    geom_col(width= 0.6,color = "black")+
    geom_text +
    scale_x_discrete(expand=c(0,1.5)) +
    scale_fill_manual(values=rev(colours)) + # 配置颜色
    theme(panel.background = element_blank(), panel.border = element_blank(), # 去掉背景颜色和边框
          axis.title.x=element_blank(), axis.title.y=element_blank(), #去掉轴标题
          axis.text.x=element_blank(),axis.text.y=element_blank(), #去掉轴刻度值
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank(), #去掉轴刻度线
          #plot.margin=unit(rep(0.3,4),'cm'),
          legend.key.height = unit(5,'pt'),legend.key.width = unit(8,'pt')) +
    labs(x = NULL, y = NULL, fill = 'Cell Type') +
    coord_polar("y")
  cat('MutiDonutplot Finish!\n')
  return(mutidonut_p)
}


#' @title Cell type ration visualization in multidonut plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_sample reorder the samples.
#' @param order_type reorder the cell type.
#' @param min_percent the minimum percent. use to combine some cell type with low ratio, default is 0.
#' @param lab whether to show the label on the plot, default is TRUE.
#'
#' @return mutidonut A ggplot2 object containing the plot.
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
#' plot_mutidonut(input, group_by = 'sample',cell_type = 'cluster.type')

# plot_mutidonut main function
plot_mutidonut <- function(input, group_by=NULL, cell_type=NULL, sub_sample='all', sub_type='all', order_sample='default', order_type = 'default',min_percent=0,lab=T){
  prerared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  data_new <- show_others(prerared,min_percent)
  mutidonut <- mutidonut_plot(data_new, order_sample, lab)
  return(mutidonut)
}


