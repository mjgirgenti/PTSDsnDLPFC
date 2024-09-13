library(stringr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(glue)
library(data.table)

construct_metacell_new <- function(seurat_obj, name='agg', max_overlap = 0.8, K=50){
  nn_map <- seurat_obj@neighbors$RNA.nn@nn.idx
  cell_used = rep(0, dim(nn_map)[1]) #count how many cells are included
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1,
      replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0
  k2 <- K * 2
  get_shared <- function(other, this_choice) {
      k2 - length(union(cell_sample[other, ], this_choice))
  }
  while (length(good_choices) > 0 & it < length(cell_used)/((1-max_overlap)*K)) {
      it <- it + 1
      choice <- sample(seq_len(length(good_choices)), size = 1,
          replace = FALSE)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[good_choices != good_choices[choice]]
      cell_sample <- nn_map[new_chosen, ]
      others <- seq_len(nrow(cell_sample) - 1)
      this_choice <- cell_sample[nrow(cell_sample), ]
      shared <- sapply(others, get_shared, this_choice = this_choice)

      if (max(shared) < max_overlap * K) {
          chosen <- new_chosen
          cell_used[this_choice] <- 1 #note the cells used in current metacell
      }

      if (it %% 1000 == 0)
          message(paste(it, "meta cells tested"))
  }

  cell_sample <- nn_map[chosen, ]
  combs <- combn(nrow(cell_sample), 2)
  shared <- apply(combs, 2, function(x) {
      k2 - length(unique(as.vector(cell_sample[x, ])))
  })

  message(paste0("Overlap QC metrics:\nCells per bin: ",
                K, "\nMaximum shared cells bin-bin: ", max(shared),
                "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
                median(shared), "\nInput number of cells: ", length(cell_used),
                "\nNumber of cells captured: ", sum(cell_used),
                "\nNumber of metacells constructed: ", dim(cell_sample)[1]))
  if (mean(shared)/K > 0.1)
      warning("On average, more than 10% of cells are shared between paired bins.")

  exprs_old <- seurat_obj@assays$RNA@data

  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
      cell_sample[x, , drop = FALSE])
  mask <- Matrix::Matrix(mask)
  new_exprs <- (exprs_old %*% mask) / K
  colnames(new_exprs) <- paste0(name, '_', 1:ncol(new_exprs))
  rownames(cell_sample) <- paste0(name, '_', 1:ncol(new_exprs))
  colnames(cell_sample) <- paste0('knn_', 1:ncol(cell_sample))

  # make seurat obj:
  seurat_aggr <- CreateSeuratObject(
    counts = new_exprs
  )
  return(list(seurat_aggr, cell_sample))
}