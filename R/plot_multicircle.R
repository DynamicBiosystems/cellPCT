#' @title multicircle_plot function
#'
#' @param data input data, Seurat obj or data.frame can be use.
#' @param outdir output dir for multicircle_plot Eg. '/path/to/save/plot', default is './'
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all
#' @import circlize
#' @import gridBase
#' @import ComplexHeatmap
#' @import grid

# multicircle_plot function
multicircle_plot <- function(plot_data, outdir){
  # legend prepare
  lgd_1 <- Legend(
    at = c("Number of cells"),
    type = "grid",
    legend_gp = gpar(fill = '#CCCCFF'),
    title_position = "topleft",
    title = "Legend")

  lgd_2 <- Legend(
    at = c("Pct of cell type",'Pct of other types'),
    type = "grid",
    legend_gp = gpar(fill = c('#FEA1A1','#9AC5F4')))

  lgd_3 <- Legend(
    at = names(plot_data[['merge_data']]),
    type = "grid",
    legend_gp = gpar(fill = c4a('pastel',length(names(plot_data[['merge_data']])))),
    title_position = "topleft",
    title = "Groups")
  # legend merge
  merged_legend <- packLegend(lgd_1, lgd_2,lgd_3,row_gap = unit(0.2, "cm"))
  # multicircle plot
  options(warn = -1)
  pdf(paste0(outdir,'/Cell_Type_Multicircle.pdf'), width = 13, height = 13)
  copy <- dev.cur() # save plot
  png(paste0(outdir,'/Cell_Type_Multicircle.png'), res = 300, width = 3200, height = 3200)
  dev.control("enable")
  circle_size = unit(1, 'snpc')
  circos.clear()
  circos.par(gap.degree = 0.5, start.degree = 90)
  circos.initialize(sectors=plot_data[['circle_1']]$type, xlim = as.matrix(plot_data[['circle_1']][,c(2,3)]))
  circos.track(
    ylim = c(0, 1), track.height = 0.07, bg.border = F, bg.col = c4a('dynamic',dim(plot_data[['merge_data']])[1]),
    panel.fun = function(x, y) {
      ylim = get.cell.meta.data('ycenter')#ylim、xlim
      xlim = get.cell.meta.data('xcenter')
      sector.name = get.cell.meta.data('sector.index')#sector.name
      brk_lab <- c(10,100,1000,10000)
      circos.axis(h = 'top', major.at = plot_data[['brk']],labels = brk_lab,  labels.cex = 0.8, labels.niceFacing = TRUE)
      circos.text(xlim, ylim, sector.name, cex = 1, niceFacing = TRUE)
    })
  circos.genomicTrackPlotRegion(plot_data[['merge_scale_circle_2']], track.height = 0.1, bg.border = NA, stack = TRUE,#圈图的高度、颜色等设置
                                panel.fun = function(region, value, ...) {
                                  circos.genomicRect(region, value, col = c('#CCCCFF'), border = NA, ...)#区块的长度反映了富集基因的数量，颜色与 p 值有关
                                  ylim = get.cell.meta.data('ycenter')
                                  xlim = plot_data[['merge_scale_circle_2']][get.cell.meta.data('sector.numeric.index'),3]/2
                                  sector.name = plot_data[['merge_scale_circle_2']][get.cell.meta.data('sector.numeric.index'),4]
                                  circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = TRUE)
                                })
  circos.genomicTrackPlotRegion(plot_data[['circle_3']], track.height = 0.09, bg.border = NA, stack = TRUE,
                                panel.fun = function(region, value, ...) {
                                  circos.genomicRect(region, value, col = c('#FFB6C1','#87CEFA'), border = NA, ...)
                                  ylim = get.cell.meta.data('cell.bottom.radius') - 0.5
                                  xlim = plot_data[['circle_3_1']][get.cell.meta.data('sector.numeric.index'),3]/2
                                  sector.name = plot_data[['circle_3_1']][get.cell.meta.data('sector.numeric.index'),4]
                                  circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = TRUE)#将文字标签（上调基因比例）添加在图中指定位置处
                                  xlim = (plot_data[['circle_3_2']][get.cell.meta.data('sector.numeric.index'),2] + plot_data[['circle_3_2']][get.cell.meta.data('sector.numeric.index'),3])/2
                                  sector.name = plot_data[['circle_3_2']][get.cell.meta.data('sector.numeric.index'),4]
                                  circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = TRUE)
                                })
  circos.genomicTrackPlotRegion(
    plot_data[['circle_4']], ylim = c(0, max(plot_data[['circle_4']]$ratio)+0.05), track.height = 0.4, bg.col = NA, bg.border = NA,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,col = c4a('pastel',length(names(plot_data[['merge_data']]))), border = NA, ytop.column = 1, ybottom = 0, ...)
      n = 1
      gap <- length(unique(plot_data[['circle_4']]$plot_cell_type))
      for (i in seq(1,length(names(plot_data[['merge_data']])),1)){
        j = n + gap -1
        ylim = plot_data[['circle_4']][c(n:j),][get.cell.meta.data('sector.numeric.index'),4] + 0.03
        xlim = (plot_data[['circle_4']][c(n:j),][get.cell.meta.data('sector.numeric.index'),2] + plot_data[['circle_4']][c(n:j),][get.cell.meta.data('sector.numeric.index'),3])/2
        sector.name = plot_data[['circle_4']][c(n:j),][get.cell.meta.data('sector.numeric.index'),4]
        circos.text(xlim, ylim, sector.name, cex = 0.8, niceFacing = TRUE)
        n = j +1
      }
    })
  # add legend
  draw(merged_legend)
  dev.copy(which = copy)
  dev.off()
  dev.off()
  cat('MutiCircle Finish!\n')
}


#' @title Cell type ration visualization in multicircle plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_type reorder the cell type.
#' @param outdir output dir for multicircle_plot Eg. '/path/to/save/plot', default is './'
#'
#' @export
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import cols4all
#' @import circlize
#' @import gridBase
#' @import ComplexHeatmap
#' @import grid
#'
#' @examples
#'plot_multicircle(input, group_by = 'sample',cell_type = 'cluster.type')

plot_multicircle <- function(input, group_by=NULL, cell_type=NULL,sub_sample='all', sub_type='all', order_type = 'default', outdir='.'){
  prepared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  multicirlce_data <- prepare_multicircle(prepared)
  multicircle_plot(multicirlce_data, outdir)
}



