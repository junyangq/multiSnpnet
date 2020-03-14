library(reshape2)
library(dplyr)
library(ggplot2)

###################################################

#  Main Parameters

# list of directories of the results
# result file is saved at, e.g. parent_dir[i]/${prefix_dir}[i]${ranks}[j]/results_rank_${rank}[j]/output_lambda_*.RData
# parent_dir[i] saves results of type fit_types[i] ("lazy"/"exact")
parent_dir <- c("/scratch/users/junyangq/multiresponse/results_mr_lazy_1_rank_35/",
                "/oak/stanford/groups/mrivas/users/mrivas/repos/multiresponse-ukbb/"
)
prefix_dir <- c("results_rank_", "results_biomarkers_unweightedresults_rank_")
fit_types <- c("lazy", "exact")

# directory of snpnet results (otherwise use NULL)
# each phenotype should have its own subdirectory,
# whose name is consistent with the one found in the multi-response results;
# result file is saved at ${snpnet_dir}/results/results/output_iter_*.RData
snpnet_dir <- "/oak/stanford/groups/mrivas/projects/biomarkers/snpnet/biomarkers/"

ranks <- c(5, 10, 20, 35)  # ranks to be plotted

# output directory of plots
plot_dir <- "/scratch/users/junyangq/multiresponse/plots_biomarkers/"

R2_lim <- c(0, NA)  # plot limit of validation R2
###################################################

dir.create(plot_dir)
is_multiplot <- FALSE
prefix_result_file <- "output_lambda_"
suffix_result_file <- ".RData"
train_name <- "metric_train"
val_name <- "metric_val"
multiplot_col <- 3  # number of columns in the multiplot

data_metric_full <- NULL

for (dir_idx in 1:length(prefix_dir)) {
  for (rank in ranks) {
    dir_rank <- file.path(parent_dir[dir_idx], paste0(prefix_dir[dir_idx], rank), paste0("results_rank_", rank))
    files_in_dir <- list.files(dir_rank)
    result_files <- files_in_dir[startsWith(files_in_dir, prefix_result_file)]
    max_iter <- max(as.numeric(gsub(suffix_result_file, "", gsub(pattern = prefix_result_file, "", result_files))))
    latest_result <- file.path(dir_rank, paste0(prefix_result_file, max_iter, suffix_result_file))
    
    myenv <- new.env()
    load(latest_result, envir = myenv)
    metric_train <- myenv[[train_name]]
    metric_val <- myenv[[val_name]]
    imax_train <- max(which(apply(metric_train, 1, function(x) sum(is.na(x))) == 0))
    imax_val <- max(which(apply(metric_val, 1, function(x) sum(is.na(x))) == 0))
    imax <- min(imax_train, imax_val)
    metric_train <- metric_train[1:imax, , drop = F]
    metric_val <- metric_val[1:imax, , drop = F]
    metric_train <- cbind(metric_train, lambda = 1:imax)
    metric_val <- cbind(metric_val, lambda = 1:imax)
    
    table_train <- melt(as.data.frame(metric_train), id.vars = "lambda", variable.name = "phenotype", value.name = "metric_train")
    table_val <- melt(as.data.frame(metric_val), id.vars = "lambda", variable.name = "phenotype", value.name = "metric_val")
    data_metric <- inner_join(table_train, table_val, by = c("phenotype", "lambda"))
    data_metric[["type"]] <- fit_types[dir_idx]
    data_metric[["rank"]] <- factor(rank, levels = as.character(ranks))
    
    data_metric_full <- rbind(data_metric_full, data_metric)
  }
}

if (!is.null(snpnet_dir)) {
  for (phe in as.character(unique(data_metric[["phenotype"]]))) {
    print(phe)
    phe_dir <- file.path(snpnet_dir, phe, "results", "results")
    files_in_dir <- list.files(phe_dir)
    result_files <- files_in_dir[startsWith(files_in_dir, "output_iter_") & endsWith(files_in_dir, ".RData")]
    max_iter <- max(as.numeric(gsub(suffix_result_file, "", gsub(pattern = "output_iter_", "", result_files))))
    latest_result <- file.path(phe_dir, paste0("output_iter_", max_iter, suffix_result_file))
    
    myenv <- new.env()
    load(latest_result, envir = myenv)
    metric_train <- myenv[["metric.train"]]
    metric_val <- myenv[["metric.val"]]
    
    imax_train <- max(which(!is.na(metric_train)))
    imax_val <- max(which(!is.na(metric_val)))
    imax <- min(imax_train, imax_val)
    
    table_snpnet <- data.frame(lambda = 1:imax, phenotype = rep(phe, imax), metric_train = metric_train[1:imax],
                               metric_val = metric_val[1:imax], type = "exact", rank = "snpnet")
    
    data_metric_full <- rbind(data_metric_full, table_snpnet)
  }
}


if (is_multiplot) {
  ggplot(data_metric_full, aes(x = metric_train, y = metric_val, colour = type)) +
    geom_line() + geom_point() +
    facet_wrap(. ~ phenotype, ncol = multiplot_col, scales = "free") +
    xlab("metric (train)") + ylab("metric (val)") +
    ylim(c(0, NA)) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
          legend.text=element_text(size=12), legend.title = element_text(size=12),
          legend.position = "bottom",
          strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
  save_path <- file.path(plot_dir, "metric_plot.pdf")
  ggsave(save_path)
} else {
  for (phe in as.character(unique(data_metric[["phenotype"]]))) {
    ggplot(dplyr::filter(data_metric_full, phenotype == phe), aes(x = metric_train, y = metric_val, shape = type, colour = rank)) +
      geom_path() + geom_point() +
      xlab("metric (train)") + ylab("metric (val)") +
      ylim(R2_lim) + 
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12),
            legend.text=element_text(size=12), legend.title = element_text(size=12),
            legend.position = "bottom",
            strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12)) +
      ggtitle(phe)
    save_path <- file.path(plot_dir, paste0("metric_plot_", phe, ".pdf"))
    ggsave(save_path)
  }
}