DA.mascp = function (data, predictor, paired = NULL, covars = NULL, out.all = NULL, 
                    p.adj = "fdr", coeff = 2, coeff.ref = 1, allResults = FALSE, 
                    ...) 
{
  ok <- tryCatch({
    loadNamespace("Maaslin2")
    TRUE
  }, error = function(...) FALSE)
  if (ok) {
    if (is(data, "phyloseq")) {
      DAdata <- DA.phyloseq(data, predictor, paired, covars)
      count_table <- DAdata$count_table
      predictor <- DAdata$predictor
      paired <- DAdata$paired
      covars <- DAdata$covars
    }
    else {
      count_table <- data
    }
    predictor <- as.factor(predictor)
    if (coeff == coeff.ref) 
      stop("coeff and coeff.ref cannot be the same")
    if (!coeff %in% seq_along(unique(predictor)) | !coeff.ref %in% 
        seq_along(unique(predictor))) 
      stop(paste("coeff and coeff.ref should be integers between 1 and", 
                 length(unique(predictor))))
    if (is.null(out.all)) {
      if (length(unique(predictor)) == 2) 
        out.all <- FALSE
      if (length(unique(predictor)) > 2) 
        out.all <- TRUE
      if (is.numeric(predictor)) 
        out.all <- FALSE
    }
    if (is.null(paired)) {
      if (is.null(covars)) {
        predictordf <- data.frame(predictor = factor(predictor))
        row.names(predictordf) <- colnames(count_table)
        x <-  as.data.frame(t(count_table))
      }
      else {
        predictordf <- as.data.frame(c(list(predictor = factor(predictor)), 
                                       covars))
        row.names(predictordf) <- colnames(count_table)
        x <- as.data.frame(t(count_table))
      }
    }
    else {
      if (is.null(covars)) {
        predictordf <- data.frame(predictor = factor(predictor), 
                                  paired = factor(paired))
        row.names(predictordf) <- colnames(count_table)
        x <- as.data.frame(t(count_table))
      }
      else {
        predictordf <- as.data.frame(c(list(predictor = factor(predictor), 
                                            paired = factor(paired)), covars))
        row.names(predictordf) <- colnames(count_table)
        x <- as.data.frame(t(count_table))
      }
    }
    x <- Maaslin2::Maaslin2(
      input_data = x,
      input_metadata = predictordf,
      output = "./",
      #min_abundance = 0.0001,
      #min_prevalence = 0.1,
      normalization = "NONE",
      transform = "NONE",
      analysis_method = "CPLM",
      max_significance = 0.1,
      fixed_effects = colnames(predictordf),
      #random_effects = ,
      correction = "BH",
      cores=1,
      standardize = FALSE,
      plot_scatter = FALSE,
      plot_heatmap = FALSE)
    
    res <- x$results[x$results$metadata == "predictor",]
    res$ordering <- NA
    res[!is.na(res$coef) & res$coef > 
          0, "ordering"] <- paste0(levels(as.factor(predictor))[coeff], 
                                   ">", levels(as.factor(predictor))[coeff.ref])
    res[!is.na(res$coef) & res$coef < 
          0, "ordering"] <- paste0(levels(as.factor(predictor))[coeff.ref], 
                                   ">", levels(as.factor(predictor))[coeff])
    colnames(res)[6] <- "pval"
    
    res <- res[, -which(names(res) %in% c("qval",",metadata"))]
    res$pval.adj <- p.adjust(res$pval, method = p.adj)
    names(res)[names(res) == 'feature'] <- 'Feature'
    
    if (file.exists("all_results.tsv")) {
      unlink("all_results.tsv")
    }
    if (file.exists("significant_results.tsv")) {
      unlink("all_results.tsv")
    }
    
    
    res$Method <- "Maaslin2 (CPLM) "
    if (is(data, "phyloseq")) 
      res <- addTax(data, res)
    if (allResults) 
      return(x)
    else return(res)
  }
  else {
    stop("Maaslin2 package required")
  }
}
