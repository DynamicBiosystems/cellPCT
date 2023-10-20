#' @title rose_plot function.
#'
#' @param data input data, Seurat obj or data.frame can be use.
#' @param lab whether to show the label on the plot, default is TRUE.
#'
#' @return rose_p A list contain ggplot2 objects
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all
#'

# rose_plot function
rose_plot <- function(data, lab){
  rose_p <- list()
  for (sample in names(data)) {
    type_n <- length(unique(data[[sample]]$plot_cell_type))
    colours <- c4a('dynamic',type_n)
    sample_name <- unique(data[[sample]]$plot_sample)
    if (min(data[[sample]]$log)>2.5){
      min_y <- min(data[[sample]]$log) - 2.5
    } else if (min(data[[sample]]$log)>1.5){
      min_y <- min(data[[sample]]$log) - 1
    } else if (min(data[[sample]]$log)>0.5){
      min_y <- min(data[[sample]]$log) - 0.5
    }
    if (lab == T){
      geom_text <- geom_text(aes(label=nums),size=2,nudge_y = 0.9)
    } else {
      geom_text <- NULL
    }
    rose_p[[sample]] <- ggplot(data[[sample]], aes(x = order, y = log)) +
      geom_col(aes(fill = label), width = 1, size = 0) +
      geom_col(aes(y=0),fill = "white", width = 1, alpha = 1, size = 0) +
      geom_col(aes(y=min_y),fill = "white", width =1, alpha = 1, size = 0) +
      coord_polar() +
      #scale_y_continuous(limits = c(0,type_n)) + # 需根据cell type数量判断
      theme(panel.background = element_blank(), panel.border = element_blank(), # 去掉背景颜色和边框
            axis.title.x=element_blank(), axis.title.y=element_blank(), #去掉轴标题
            axis.text.x=element_blank(),axis.text.y=element_blank(), #去掉轴刻度值
            axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
            legend.key.height = unit(5,'pt'),legend.key.width = unit(8,'pt'),
            plot.title = element_text(hjust = 0.5)) + #去掉轴刻度线
      labs(x = NULL, y = NULL, fill = 'Cell Types',title = sample_name) +
      geom_text +
      scale_fill_manual(values=colours) # 配色选取
  }
  cat('Rose Plot Finish!\n')
  return(rose_p)
}


#' @title Cell type ration visualization in rose plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_type reorder the cell type.
#' @param lab whether to show the label on the plot, default is TRUE.
#'
#' @return rose A list contain ggplot2 objects
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
#' plot_rose(input, group_by = 'sample',cell_type = 'cluster.type')

# plot_rose main function
plot_rose <- function(input, group_by=NULL, cell_type=NULL,sub_sample='all', sub_type='all', order_type = 'default',lab=T){
  prerared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  rose <- rose_plot(prerared, lab)
  return(rose)
}
