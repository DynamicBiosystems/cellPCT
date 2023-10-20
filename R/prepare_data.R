#' @title Prepare plot data for multi-visualization functions.
#' @param data input data, Seurat obj or data.frame can be use.
#' @param group_by colname of Seurat metadata, which represent the groups or samples.
#' @param cell_type colname of Seurat metadata, which represent the cell type.
#' @param sub_sample select part of the groups for plot, default is 'all'.
#' @param sub_type select part of the cell types for plot, default is 'all'.
#' @param order_type reorder the cell type.
#'
#' @return sub_prepared_cell
#'
#' @import dplyr
#' @import SeuratObject
#' @import forcats

prepare_data <- function(data, group_by, cell_type, sub_sample, sub_type, order_type){
  if ( class(data) %in% c('Seurat','data.frame')){
    cat('Input data processing ...\n')
    if (class(data) == 'Seurat'){
      cat('Use Seurat object as input data\n')
      # prepare rds
      metadata <- data@meta.data
      cell_type_input <- cell_type
      sample_input <- group_by
      # check input group_by and cell_type colnames in Seurat metadata and rename
      if (cell_type_input %in% colnames(metadata)){
        cat('cell_type Check Finish\n')
        colnames(metadata)[colnames(metadata) %in% cell_type_input] <- c('plot_cell_type')
      } else {
        cat(paste0('Please check your input cell_type : ',cell_type_input,' Not Found!\n'))
      }
      if (sample_input %in% colnames(metadata)){
        cat('group_by Check Finish\n')
        colnames(metadata)[colnames(metadata) %in% sample_input] <- c('plot_sample')
      } else {
        cat(paste0('Please check your input group_by : ',sample_input,' Not Found!\n'))
      }
      sample_split <- list()
      sample_ratio <- list()
      for (i in as.character(unique(metadata$plot_sample))){
        sample_split[[i]] <- metadata[metadata$plot_sample == i,c('seurat_clusters','plot_cell_type')]
        sample_ratio[[i]] <- table(sample_split[[i]]$plot_cell_type) %>% as.data.frame() %>%
          mutate(plot_sample = i) %>% mutate(ratio = Freq/sum(Freq)) %>% rename( nums = Freq) %>% rename( plot_cell_type = Var1)
        sample_ratio[[i]] <- sample_ratio[[i]][order(-sample_ratio[[i]]$ratio),] %>%
          mutate(order = fct_reorder(plot_cell_type,nums)) %>% mutate(log = log(nums + 1))
        cell_type <- sample_ratio[[i]]$plot_cell_type
        cell_ratio <- round(sample_ratio[[i]]$ratio*100, digits = 2)
        sample_ratio[[i]]$label <- factor(sample_ratio[[i]]$plot_cell_type,
                                          levels = sample_ratio[[i]]$plot_cell_type,
                                          labels = paste0(cell_type,' (',cell_ratio,'%)'))
      }
    } else if (class(data) == 'data.frame'){
      cat('Use other file as input data\n')
      # prepare txt
      for (i in colnames(data)){
        if (class(data[,which(colnames(data) == i)]) == 'character'){
          data <- data %>% rename( plot_cell_type = i)
        }
      }
      sample_split <- list()
      sample_ratio <- list()
      for (i in colnames(data)[! colnames(data) %in% c('plot_cell_type')]){
        sample_split[[i]] <- data %>% select(plot_cell_type,i) %>% mutate(plot_sample = as.character(i))
        sample_ratio[[i]] <- sample_split[[i]] %>% mutate(ratio = !!sym(i)/colSums(across(where(is.numeric)))) %>% rename( nums = i)
        sample_ratio[[i]] <- sample_ratio[[i]][order(-sample_ratio[[i]]$ratio),] %>%
          mutate(order = fct_reorder(plot_cell_type,nums)) %>% mutate(log = log(nums + 1))
        cell_type <- sample_ratio[[i]]$plot_cell_type
        cell_ratio <- round(sample_ratio[[i]]$ratio*100, digits = 2)
        sample_ratio[[i]]$label <- factor(sample_ratio[[i]]$plot_cell_type,
                                          levels = sample_ratio[[i]]$plot_cell_type,
                                          labels = paste0(cell_type,' (',cell_ratio,'%)'))
      }
    }
    # sub_sample & sub_type
    # sub_sample select
    all_sample <- names(sample_ratio)
    sub_prepared <- list()
    if ('all' %in% sub_sample){
      sub_prepared <- sample_ratio
    } else {
      sub_sample = sub_sample
      for (i in names(sample_ratio)[which(names(sample_ratio) %in% sub_sample)]){
        sub_prepared[[i]] <- sample_ratio[[i]]
      }
    }
    # sub_type select
    if ('all' %in% sub_type) {
      #sub_type = all_type
      sub_prepared_cell <- sub_prepared
    } else {
      sub_type = sub_type
      sub_prepared_cell <- lapply(names(sub_prepared),function(sample){
        sub_prepared[[sample]] <- sub_prepared[[sample]] %>% filter(plot_cell_type %in% sub_type)
      })
      names(sub_prepared_cell) <- names(sub_prepared)
    }
    if ('default' %in% order_type) {
      sub_prepared_cell <- sub_prepared_cell
    } else {
      for (i in names(sub_prepared_cell)) {
        sub_prepared_cell[[i]] <- sub_prepared_cell[[i]] %>% mutate(plot_cell_type = fct_relevel(plot_cell_type, order_type)) %>% arrange(plot_cell_type)
      }
    }

    cat('Plot data already prepared\n')
  } else {
    cat('Error! Please check your input data type!\n')
  }
  return(sub_prepared_cell)
}

#' @title Select part of the Prepare plot data for multi-visualization functions.
#' @param sample_ratio the return value of prepare_data.
#' @param min_percent the minimum percent. use to combine some cell type with low ratio,default is 0.
#'
#' @return sample_new
#'
#' @import dplyr
#' @import SeuratObject
#' @import forcats

show_others <- function(sample_ratio, min_percent){
  sample_new <- list()
  if (min_percent >0) {
    cat(paste0('Use cell types that ratio upper ',min_percent,'% to plot\n'))
    for (i in names(sample_ratio)){
      sample_ratio[[i]]$plot_cell_type <- as.character(sample_ratio[[i]]$plot_cell_type)
      sample_ratio[[i]][sample_ratio[[i]]$ratio*100 < min_percent, c('plot_cell_type')] <- c('Others')
      other_ratio <- round(sum(sample_ratio[[i]][sample_ratio[[i]]$plot_cell_type == 'Others','ratio'])*100, digits = 2)
      sample_new[[i]] <- filter(sample_ratio[[i]], plot_cell_type != 'Others') %>%
        add_row(plot_cell_type = 'Others',
                nums = sum(sample_ratio[[i]][sample_ratio[[i]]$plot_cell_type == 'Others','nums']),
                plot_sample = unique(sample_ratio[[i]]$plot_sample),
                ratio = other_ratio/100,
                label = as.factor(paste0('Ohters',' (',other_ratio,'%)')))
    }
  } else if (min_percent == 0) {
    cat('Use all kind of cell types to plot\n')
    sample_new <- sample_ratio
  } else {
    stop('Error! min_percent must >= 0\n')
  }
  return(sample_new)
}


#' @title Prepare data for radar plot
#'
#' @param prepared the return value of prepare_data.
#'
#' @return radar_data
#'
#' @import dplyr
#' @import SeuratObject
#' @import forcats
#'
prepare_radar <- function(prepared){
  for (i in names(prepared)) {
    colnames(prepared[[i]]) <- gsub('nums',i,colnames(prepared[[i]]))
    prepared[[i]] <- prepared[[i]] %>% select('plot_cell_type',i) %>% .[order(prepared[[i]]$plot_cell_type),]
  }
  merge_data <- as.data.frame(prepared,col.names=NULL,row.names = 1)
  if (dim(merge_data)[2] >1) {
    merge_data <- merge_data[,-c(seq(2,dim(merge_data)[2],2))]
  } else {
    merge_data <- merge_data
  }
  pct_merge <- lapply(names(merge_data),function(col){
    as.data.frame(merge_data[,col]/sum(merge_data[,col]),row.names= rownames(merge_data))
  })
  pct_merge <- as.data.frame(pct_merge)
  colnames(pct_merge) <- names(merge_data)
  radar_data <- as.data.frame(t(pct_merge))
  radar_data$group <- row.names(radar_data)
  radar_data <- radar_data[,c(dim(radar_data)[2],2:dim(radar_data)[2]-1)]
  return(radar_data)
}


#' @title Prepare data for multicircle plot
#'
#' @param prepared the return value of prepare_data.
#'
#' @return multicircle_data
#'
#' @import dplyr
#' @import SeuratObject
#' @import forcats
#'
prepare_multicircle <- function(prepared){
  for (i in names(prepared)) {
    colnames(prepared[[i]]) <- gsub('nums',i,colnames(prepared[[i]]))
    prepared[[i]] <- prepared[[i]] %>% select('plot_cell_type',i) %>% .[order(prepared[[i]]$plot_cell_type),]
  }
  merge_data <- as.data.frame(prepared,col.names=NULL,row.names = 1)
  if (dim(merge_data)[2] >1) {
    merge_data <- merge_data[,-c(seq(2,dim(merge_data)[2],2))]
  } else {
    merge_data <- merge_data
  }
  circle_2 <- merge_data %>% mutate(start = 0) %>% mutate(end = rowSums(.)) %>% select('start','end')
  circle_1 <- circle_2 %>% mutate(end = ceiling(max(end)/100)*100)
  circle_1 <- circle_1 %>% mutate(type = rownames(.))
  circle_1 <- circle_1[,c(dim(circle_1)[2],1:dim(circle_1)[2]-1)]
  circle_2 <- circle_2 %>% mutate(type = rownames(.))
  circle_2 <- circle_2[,c(dim(circle_2)[2],1:dim(circle_2)[2]-1)]
  gap <- unique(circle_1$end)/4
  scale_circle_2 <- list()
  for (i in seq(2,4)){
    scale_circle_2[[i]] <- circle_2 %>% mutate(cal = end/10^i) %>% filter(cal >= 1 & cal <10) %>% mutate(site = gap*(i-1) + (end-10^i)/(10^(i+1)-10^i)*gap)
  }
  merge_scale_circle_2 <- do.call(rbind,scale_circle_2) %>% select('type','start','site','end') %>% .[c(circle_1$type),]
  circle_3_1 <- merge_data %>% mutate(sum = rowSums(.)) %>% mutate(total = sum(sum)) %>%
    mutate(total_pct = round(sum/total,digits = 2))
  if ( 0 %in% circle_3_1$total_pct){
    circle_3_1[circle_3_1$total_pct == 0,c('total_pct')] <- 0.01
  }
  circle_3_1 <- circle_3_1 %>% mutate(site1 = unique(circle_1$end)*total_pct) %>%
    mutate(type = rownames(.)) %>% mutate(start = 0) %>% select('type','start','site1','total_pct')
  circle_3_2 <- circle_3_1 %>% mutate(start = site1) %>% mutate(site1 = unique(circle_1$end)) %>% mutate(total_pct = 1-total_pct)
  circle_3 <- rbind(circle_3_1,circle_3_2)
  ratio_list <- lapply(prepared,function(sample){
    sample %>% mutate(ratio = round(sample[,2]/sum(sample[,2]),digits = 2))
  })
  n = 0
  for (i in names(ratio_list)){
    gap <- unique(circle_1$end)/dim(merge_data)[2]
    n = n + gap
    ratio_list[[i]] <- ratio_list[[i]] %>% mutate(start = n-gap) %>% mutate(end = n) %>% select('plot_cell_type','start','end','ratio')
  }
  circle_4 <- do.call(rbind,ratio_list)
  max_end <- unique(circle_1$end)
  end_gap <- max_end/4
  brk <- seq(0, max_end, end_gap)
  multicircle_data <- list(merge_data = merge_data, brk = brk, circle_1 = circle_1,circle_2 = circle_2,
                    merge_scale_circle_2 = merge_scale_circle_2, circle_3 = circle_3,
                    circle_3_1 = circle_3_1, circle_3_2 = circle_3_2, circle_4 = circle_4)
  return(multicircle_data)
}



