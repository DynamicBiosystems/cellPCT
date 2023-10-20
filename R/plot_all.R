#' @title Cell type ration visualization in all plot with default paraments and save as png and pdf.
#'
#' @param input the input data rds,txt,csv can be use
#' @param outdir the output files' path
#' @param sample the colname of sample in input data
#' @param cell_type the colname of cell in input data
#' @export
#' @import dplyr
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import ggalluvial
#' @import cols4all
#' @import treemapify
#' @import ggradar
#' @import circlize
#' @import ComplexHeatmap
#' @import gridBase
#' @import grid
#'
#' @examples
#' plot_all(input=input, outdir=outdir,group_by='sample',cell_type='cluster.type')

plot_all <- function(input,outdir, group_by=c('sample'), cell_type = c('cluster.type')){
  circle <- plot_circle(input, group_by, cell_type)
  for (i in names(circle)){
    ggsave(paste0(outdir,'/','Cell_Type_Circle_',i,'.png'), circle[[i]], width=15,height=10,units = c('cm'))
    ggsave(paste0(outdir,'/','Cell_Type_Circle_',i,'.pdf'), circle[[i]], width=15,height=10,units = c('cm'))
  }
  cat('Circle plot saved\n')
  cat('\n')

  sankey <- plot_sankey(input, group_by, cell_type)
  ggsave(paste0(outdir,'/','Cell_Type_Ratio_Sankey.png'),sankey,width=12,height=8,units = c('cm'))
  ggsave(paste0(outdir,'/','Cell_Type_Ratio_Sankey.pdf'),sankey,width=12,height=8,units = c('cm'))
  cat('Sankey plot saved\n')
  cat('\n')

  mutidonut <- plot_mutidonut(input, group_by, cell_type)
  ggsave(paste0(outdir,'/','Cell_Type_Ratio_mutidonut.png'),mutidonut,width=15,height=15,units = c('cm'))
  ggsave(paste0(outdir,'/','Cell_Type_Ratio_mutidonut.pdf'),mutidonut,width=15,height=15,units = c('cm'))
  cat('Mutidonut plot saved\n')
  cat('\n')

  singledonut <- plot_singledonut(input, group_by, cell_type)
  ggsave(paste0(outdir,'/','Cell_Type_singledonut','.png'), singledonut, width=15,height=10,units = c('cm'))
  ggsave(paste0(outdir,'/','Cell_Type_singledonut','.pdf'), singledonut, width=15,height=10,units = c('cm'))
  cat('Singledonut plot saved\n')
  cat('\n')

  stackbar <- plot_stackbar(input, group_by, cell_type)
  ggsave(paste0(outdir,'/','Cell_Type_Ratio_Stackbar.png'),stackbar,width=12,height=8,units = c('cm'))
  ggsave(paste0(outdir,'/','Cell_Type_Ratio_Stackbar.pdf'),stackbar,width=12,height=8,units = c('cm'))
  cat('Stackbar plot saved\n')
  cat('\n')

  stackline <- plot_stackline(input, group_by, cell_type)
  ggsave(paste0(outdir,'/','Cell_Type_Ratio_Stackline.png'),stackline,width=10,height=5,units = c('cm'))
  ggsave(paste0(outdir,'/','Cell_Type_Ratio_Stackline.pdf'),stackline,width=10,height=5,units = c('cm'))
  cat('Stackline plot saved\n')
  cat('\n')

  rose <- plot_rose(input, group_by, cell_type)
  for (i in names(rose)){
    ggsave(paste0(outdir,'/','Cell_Type_Rose_',i,'.png'), rose[[i]], width=15,height=10,units = c('cm'))
    ggsave(paste0(outdir,'/','Cell_Type_Rose_',i,'.pdf'), rose[[i]], width=15,height=10,units = c('cm'))
  }
  cat('Rose plot saved\n')
  cat('\n')

  treemap <- plot_treemap(input, group_by, cell_type)
  for (i in names(treemap)){
    ggsave(paste0(outdir,'/','Cell_Type_Treemap_',i,'.png'), treemap[[i]], width=15,height=10,units = c('cm'))
    ggsave(paste0(outdir,'/','Cell_Type_Treemap_',i,'.pdf'), treemap[[i]], width=15,height=10,units = c('cm'))
  }
  cat('Treemap plot saved\n')
  cat('\n')

  radar <- plot_radar(input, group_by, cell_type)
  ggsave(paste0(outdir,'/','Cell_Type_Radar','.png'), radar, width=15, height=10,units = c('cm'))
  ggsave(paste0(outdir,'/','Cell_Type_Radar','.pdf'), radar, width=15, height=10,units = c('cm'))
  cat('Radar plot saved\n')
  cat('\n')

  plot_multicircle(input, group_by = group_by, cell_type = cell_type, outdir = outdir)
  cat('Multicircle plot saved\n')
  cat('\n')
  cat('Finished all plot !')
  cat('\n')
}
