#' @title radar_plot function
#'
#' @param data input data, Seurat obj or data.frame can be use.
#'
#' @return radar_p A ggplot2 object containing the plot.
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import ggradar
#' @import cols4all
#'

# radar_plot function
radar_plot <- function(data){
  colours <- c4a('dynamic',dim(data)[1])
  radar_p <-ggradar(data,
                  grid.min = 0, #网格线最小值
                  grid.mid = 0.5, #网格线均值
                  grid.max = 1, #网格线最大值
                  base.size = 10,
                  axis.label.size = 3,
                  values.radar = c("","50%", "100%"), #轴标签显示
                  group.colours = colours,
                  group.point.size = 2,#分组点大小
                  group.line.width = 1, #线宽
                  background.circle.colour = 'gray', #背景填充色
                  background.circle.transparency = 0.1, #背景填充不透明度(这里改为0可去掉背景填充)
                  legend.position = 'right', #图例位置
                  legend.text.size = 6, #图例标签大小
                  fill = TRUE, #各分组是否填充色
                  fill.alpha = 0.3 #分组填充色不透明度
  )
  cat('Radar Plot Finish!\n')
  return(radar_p)
}


#' @title Cell type ration visualization in radar plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_type reorder the cell type.
#'
#' @return radar A ggplot2 object containing the plot.
#'
#' @export
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import ggradar
#' @import cols4all
#'
#' @examples
#'plot_radar(input, group_by = 'sample',cell_type = 'cluster.type')

plot_radar <- function(input, group_by=NULL, cell_type=NULL,sub_sample='all', sub_type='all', order_type = 'default'){
  prerared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  plot_data <- prepare_radar(prerared)
  radar <- radar_plot(plot_data)
  return(radar)
}
