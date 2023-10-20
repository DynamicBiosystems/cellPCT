#' @title treemap_plot function
#'
#' @param data input data, Seurat obj or data.frame can be use.
#'
#' @return treemap_p A list contain ggplot2 objects
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import treemapify
#' @import cols4all
#'

# treemap_plot function
treemap_plot <- function(data){
  treemap_p <- list()
  for (sample in names(data)) {
    type_n <- length(unique(data[[sample]]$plot_cell_type))
    colours <- c4a('dynamic',type_n)
    sample_name <- unique(data[[sample]]$plot_sample)
    treemap_p[[sample]] <- ggplot(data[[sample]], aes(area = nums, fill = label, label = plot_cell_type, subgroup = plot_cell_type)) +
      geom_treemap() +
      scale_fill_manual(values = colours) +
      geom_treemap_subgroup_border(color = 'white') +
      geom_treemap_text(colour = "white",
                        place = "centre",
                        min.size = 0,
                        grow = T) +
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      labs(title = paste0('Cell type ratio in ',sample_name))
  }
  cat('Treemap Plot Finish!\n')
  return(treemap_p)
}

#' @title Cell type ration visualization in TreeMap plot
#'
#' @param input input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_type reorder the cell type.
#'
#' @return treemap A list contain ggplot2 objects
#'
#' @export
#'
#' @import SeuratObject
#' @import forcats
#' @import ggplot2
#' @import dplyr
#' @import ggalluvial
#' @import treemapify
#' @import cols4all
#'
#' @examples
#' plot_treemap(input, group_by = 'sample',cell_type = 'cluster.type')

# plot_treemap main function
plot_treemap <- function(input, group_by=NULL, cell_type=NULL,sub_sample='all', sub_type='all', order_type = 'default'){
  prerared <- prepare_data(input, group_by, cell_type, sub_sample, sub_type, order_type)
  treemap <- treemap_plot(prerared)
  return(treemap)
}




