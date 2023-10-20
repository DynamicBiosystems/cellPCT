#' @title singledonut_plot function
#'
#' @param data input data, Seurat obj or data.frame can be use.
#' @param order_sample reorder the samples.
#' @param lab whether to show the label on the plot, default is TRUE.
#'
#' @return singledonut_p A ggplot2 object containing the plot.
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all

# singledonut_plot function
singledonut_plot <- function(data, order_sample, lab){
  plot_data <- do.call(rbind,data)
  # change sample order
  if ('default' %in% order_sample){
    plot_data <- plot_data
  } else {
    plot_data$plot_sample <- factor(plot_data$plot_sample,levels=order_sample)
  }
  type_n <- length(unique(plot_data$plot_cell_type))
  colours <- c4a('dynamic',type_n)
  if (lab == T){
    geom_text <- geom_text(aes(label=paste0(round(ratio,digits = 4)*100,'%')),position = position_stack(vjust = 0.5),
                           colour = "black",vjust = -1,hjust=0.5,size = 3)
  } else {
    geom_text <- NULL
  }
  singledonut_p <- ggplot(plot_data, aes(x = 3, y = ratio, fill = plot_cell_type)) +
    geom_col(width = 1.5,
             color = 'white') +
    facet_grid(.~plot_sample) +
    coord_polar(theta = "y") +
    xlim(c(0.2, 3.8)) +
    scale_fill_manual(values = colours) +
    labs(x = NULL, y = NULL, fill = 'Cell Type') +
    theme_void() +
    theme(
      strip.text.x = element_text(size = 12),
      #strip.background.x = element_rect('white'),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.border = element_blank(),
      axis.title.x=element_blank(), axis.title.y=element_blank(), #去掉轴标题
      axis.text.x=element_blank(),axis.text.y=element_blank(), #去掉轴刻度值
      axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
      panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    geom_text
  cat('Single Donutplot finish!\n')
  return(singledonut_p)
}

#' @title Cell type ration visualization in singledonut plot
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
#' @return singledonut A ggplot2 object containing the plot.
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
#' plot_singledonut(input, group_by = 'sample',cell_type = 'cluster.type')

# plot_singledonut main function
plot_singledonut <- function(input, group_by=NULL, cell_type=NULL, sub_sample='all', sub_type='all', order_sample='default', order_type = 'default', min_percent=0, lab=T){
  prerared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  data_new <- show_others(prerared, min_percent)
  singledonut <- singledonut_plot(data_new, order_sample, lab)
  return(singledonut)
}
