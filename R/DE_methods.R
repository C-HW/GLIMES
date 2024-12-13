
#' Two-sample t-test for gene expressions between two groups
#'
#' @param counts gene expression matrix with gene rows and cell columns
#' @param cellgroup1 the index of cells from group1
#' @param cellgroup2 the index of cells from group2
#' @return A dataframe with t-score, p-value, and BH p-value for each gene in a row
#' @examples
#' data(Bcells_sce)
#' simple_mean_DE(Bcells_sce@assays@data$counts, Bcells_sce$stim == "ctrl", Bcells_sce$stim == "stim")
#' @export
simple_mean_DE = function(counts, cellgroup1, cellgroup2){
  tmpcount1 = counts[,cellgroup1]
  tmpcount2 = counts[,cellgroup2]
  if(length(cellgroup1) == 1){
    tmpmean1 = mean(tmpcount1)
  }else{
    tmpmean1 = Matrix::rowMeans(tmpcount1)
  }
  if(length(cellgroup2) == 1){
    tmpmean2 = mean(tmpcount2)
  }else{
    tmpmean2 = Matrix::rowMeans(tmpcount2)
  }
  df = data.frame(genes = rownames(counts))
  df$t = sapply(1:nrow(counts), FUN = function(i){stats::t.test(tmpcount1[i,], tmpcount2[i,])$statistic})
  df$pval = 2*stats::pt(-abs(df$t), df = length(cellgroup1) + length(cellgroup2) -1)
  df$BHpval = stats::p.adjust(df$pval, method = "BH")
  return(df)
}

#' An implementation Poisson glmm DE method.
#'
#' @param sce a SingleCellExperiment object with raw count matrix.
#' @param comparison group comparison variable (e.g., conditions, celltype). The variable consists of exact two different groups.
#' @param replicates donor variable.
#' @param exp_batch experimental batches.
#' @param other_fixed donor variable.
#' @param freq_expressed a threshold for gene detection rate.
#' @return A dataframe of Poisson glmm DE results with each row for a gene.
#' @examples
#' data(Bcells_sce)
#' poisson_glmm_DE(Bcells_sce[1:10,], comparison = "stim", replicates = "ind")
#' @export
poisson_glmm_DE = function(sce,
                           comparison,        # Mandatory: group comparison variable (e.g., conditions, celltype)
                           replicates,          # Mandatory: random effect for donors
                           exp_batch = NULL,     # Optional: random effect for experimental batches
                           other_fixed = NULL, # Optional: other fixed effects (e.g., sex, age)
                           freq_expressed = 0.05) {
  # Prepare the data frame
  countdf <- data.frame(comparison = as.factor(SummarizedExperiment::colData(sce)[, comparison]),
                        replicates = as.factor(SummarizedExperiment::colData(sce)[, replicates]))
  # Build fixed effects formula
  fixed_effects <- paste("comparison",
                         if (!is.null(other_fixed)) paste("+", paste(other_fixed, collapse = " + ")),
                         sep = " ")

  # Build random effects formula
  random_effects_list <- list(replicates = ~1)  # First random effect for replicates

  # Add other fixed effects if provided
  if (!is.null(other_fixed)) {
    countdf = cbind(countdf, data.frame(SummarizedExperiment::colData(sce)[, other_fixed, drop = FALSE]))
  }

  # Add experimental batch if provided
  if (!is.null(exp_batch)) {
    countdf$exp_batch <- SummarizedExperiment::colData(sce)[, exp_batch]
    random_effects_list$exp_batch <- ~1  # Add second random effect for exp_batch
  }

  # Initialize the output data frame
  df = data.frame(genes = rownames(sce), mu = NA, beta_comparison = NA,
                  log2FC = NA, sigma_square = NA, status = "done", pval = NA, BH = NA,
                  log2mean = NA, log2meandiff = NA)

  start_time <- Sys.time()  # Start time for tracking progress
  # Loop through each gene
  for(i in 1:nrow(sce)){
    # Progress tracking: Estimated time remaining
    if (i %% round(nrow(sce)/10) == 0) {  # Update every 10% iterations
      current_time <- Sys.time()
      elapsed_time <- as.numeric(difftime(current_time, start_time, units = "secs"))
      avg_time_per_iter <- elapsed_time / i
      remaining_time <- (nrow(sce) - i) * avg_time_per_iter
      cat(sprintf("Progress: %d/%d genes, Estimated remaining time: %.2f minutes\n",
                  i, nrow(sce), remaining_time / 60))
    }

    countdf$count = as.numeric(round(pmax(sce@assays@data$counts[i,],0)))

    # Compute gene mean and mean difference for each comparison group
    genemean = stats::aggregate(count ~ comparison, data = countdf, FUN = mean, na.rm = TRUE)
    genemean = genemean[order(genemean$comparison), ]
    genemean1 = genemean[1,2]
    genemean2 = genemean[2,2]

    df$log2mean[i] = log2(genemean1*genemean2)/2
    df$log2meandiff[i] = log2(abs(genemean1-genemean2))

    # Skip lowly expressed genes
    if (mean(countdf$count != 0, na.rm = TRUE) <= freq_expressed) {
      df$status[i] <- ifelse(mean(countdf$count != 0, na.rm = TRUE) == 0, "zero mean", "lowly expressed")
      next
    }

    # Fit Poisson GLMM model using glmmPQL
    gm = tryCatch(
      summary(MASS::glmmPQL(as.formula(paste("count ~", fixed_effects)),
                            random = random_effects_list,
                            family = stats::poisson, data = countdf,
                            verbose = FALSE)),
      error = function(e){NULL}
    )

    if (is.null(gm)){
      df$status[i] = "not converge"
      next
    }
    df$pval[i] = gm$tTable[2 ,"p-value"]
    df$sigma_square[i] = gm$sigma^2
    df$mu[i] = gm$coefficients$fixed[1]
    df$beta_comparison[i] = gm$coefficients$fixed[2]
  }
  df$log2FC = log2(exp(df$beta_comparison))
  df$BH = stats::p.adjust(df$pval, method = "BH")
  return(df)
}

#' An implementation of Binomial glmm DE method.
#'
#' @param sce a SingleCellExperiment object with raw count matrix.
#' @param comparison group comparison variable (e.g., conditions, celltype). The variable consists of exact two different groups.
#' @param replicates donor variable.
#' @param exp_batch experimental batches.
#' @param other_fixed donor variable.
#' @param freq_expressed a threshold for gene detection rate.
#' @return A dataframe of Poisson glmm DE results with each row for a gene.
#' @examples
#' data(Bcells_sce)
#' binomial_glmm_DE(Bcells_sce[1:10,], comparison = "stim", replicates = "ind")
#' @export
binomial_glmm_DE = function(sce,
                            comparison,        # Mandatory: group comparison variable (e.g., conditions, cell type)
                            replicates,        # Mandatory: random effect for donors (patients)
                            exp_batch = NULL,  # Optional: random effect for experimental batches
                            other_fixed = NULL, # Optional: other fixed effects (e.g., sex, age)
                            freq_expressed = 0.05) {

  # Prepare the data frame
  countdf <- data.frame(comparison = as.factor(SummarizedExperiment::colData(sce)[, comparison]),
                        replicates = as.factor(SummarizedExperiment::colData(sce)[, replicates]))
  # Build fixed effects formula
  fixed_effects <- paste("comparison",
                         if (!is.null(other_fixed)) paste("+", paste(other_fixed, collapse = " + ")),
                         sep = " ")

  # Build random effects formula
  random_effects_list <- list(replicates = ~1)  # First random effect for replicates

  # Add other fixed effects if provided
  if (!is.null(other_fixed)) {
    countdf = cbind(countdf, data.frame(SummarizedExperiment::colData(sce)[, other_fixed, drop = FALSE]))
  }
  # Add experimental batch if provided
  if (!is.null(exp_batch)) {
    countdf$exp_batch <- SummarizedExperiment::colData(sce)[, exp_batch]
    random_effects_list$exp_batch <- ~1  # Add second random effect for exp_batch
  }
  # Initialize the output data frame
  df = data.frame(genes = rownames(sce), mu = NA, beta_comparison = NA,
                  log2FC = NA, sigma_square = NA, status = "done", pval = NA, BH = NA,
                  log2mean = NA, log2meandiff = NA)

  start_time <- Sys.time()  # Start time for tracking progress
  # Loop through each gene
  for (i in 1:nrow(sce)) {
    # Progress tracking: Estimated time remaining
    if (i %% round(nrow(sce)/10) == 0) {  # Update every 10% iterations
      current_time <- Sys.time()
      elapsed_time <- as.numeric(difftime(current_time, start_time, units = "secs"))
      avg_time_per_iter <- elapsed_time / i
      remaining_time <- (nrow(sce) - i) * avg_time_per_iter
      cat(sprintf("Progress: %d/%d genes, Estimated remaining time: %.2f minutes\n",
                  i, nrow(sce), remaining_time / 60))
    }

    countdf$count = as.numeric(1*(sce@assays@data$counts[i,]>0))

    # Compute gene mean and mean difference for each comparison group
    genemean = stats::aggregate(count ~ comparison, data = countdf, FUN = mean, na.rm = TRUE)
    genemean = genemean[order(genemean$comparison), ]
    genemean1 = genemean[1,2]
    genemean2 = genemean[2,2]
    df$log2mean[i] = log2(genemean1*genemean2)/2
    df$log2meandiff[i] = log2(abs(genemean1-genemean2))

    # Skip lowly expressed genes
    if (mean(countdf$count) <= freq_expressed) {
      df$status[i] <- ifelse(mean(countdf$count != 0, na.rm = TRUE) == 0, "zero mean", "lowly expressed")
      next
    }
    # genemean_bygroup = aggregate(count ~ cellgroups, data = countdf, FUN = mean)
    #  if (genemean_bygroup$count[1] == genemean_bygroup$count[2]){
    #    df$status[i] = "no difference between groups"
    #    next
    #  }

    # Fit Binomial GLMM model using glmmPQL
    gm = tryCatch(
      summary(MASS::glmmPQL(as.formula(paste("count ~", fixed_effects)),
                            random = random_effects_list,
                            family = stats::poisson, data = countdf,
                            verbose = FALSE,
                            niter = 50)),
      error = function(e){NULL}
    )

    if (is.null(gm)){
      df$status[i] = "not converge"
      next
    }
    df$pval[i] = gm$tTable[2 ,"p-value"]
    df$sigma_square[i] = gm$sigma^2
    df$mu[i] = gm$coefficients$fixed[1]
    df$beta_comparison[i] = gm$coefficients$fixed[2]
  }
  df$log2FC = log2(exp(df$beta_comparison))
  df$BH = stats::p.adjust(df$pval, method = "BH")
  return(df)
}

#' Identify DEGs for a list of genes after performing DE analysis.
#'
#' An old framework select genes with adjusted p-values smaller than
#' a threshold and absolute log2 fold change greater than a threshold.
#' A new framework filters out genes with small average log2 gene means,
#' but genes showing large difference in mean would be considered as a
#' candidate for DEGs.
#'
#' @param adj_pval a vector of adjusted p-values obtained from a DE analysis
#' @param log2FC a vector of log2 fold change obtained from a DE analysis
#' @param log2mean a vector of log2(genemean1*genemean2)/2 with genemean1 and genemean2 representing the gene mean from raw counts
#' @param log2meandiff a vector of log2(abs(genemean1-genemean2)) with genemean1 and genemean2 representing the gene mean from raw counts
#' @param pvalcutoff the p-value threshold to determine DEGs
#' @param log2FCcutoff the log2 fold change threshold to determine DEGs
#' @param log2meancutoff the log2 mean threshold to determine DEGs
#' @param log2meandiffcutoff the log2 difference of gene mean threshold to determine DEGs
#' @param newcriteria logical. Whether the gene mean and difference of mean should be included in the criteria
#' @return A logical vector indicating DEGs
#' @examples
#' identifyDEGs(runif(1000), rnorm(1000), newcriteria = FALSE)
#' @export
identifyDEGs = function(adj_pval, log2FC, log2mean = NA, log2meandiff = -Inf,
                        pvalcutoff = 0.05, log2FCcutoff = log2(1.5),
                        log2meancutoff = -2.25, log2meandiffcutoff = -1, newcriteria = T){
  if(newcriteria){
    DEGs = adj_pval<pvalcutoff & abs(log2FC)>log2FCcutoff & (log2mean > log2meancutoff | log2meandiff > log2meandiffcutoff)
  }else{
    DEGs = adj_pval<pvalcutoff & abs(log2FC)>log2FCcutoff
  }
  DEGs = ifelse(is.na(adj_pval), NA, DEGs)
  return(DEGs)
}

