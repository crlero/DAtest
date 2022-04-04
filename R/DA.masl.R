#' Maaslin2 with LM (no normalisation / no transformation)
#'
#' Implementation of DESeq2 for \code{DAtest}
#' Manual geometric means calculated to avoid errors, see https://github.com/joey711/phyloseq/issues/387
#' @param data Either a matrix with counts/abundances, OR a \code{phyloseq} object. If a matrix/data.frame is provided rows should be taxa/genes/proteins and columns samples
#' @param predictor The predictor of interest. Either a Factor or Numeric, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param paired For paired/blocked experimental designs. Either a Factor with Subject/Block ID for running paired/blocked analysis, OR if \code{data} is a \code{phyloseq} object the name of the variable in \code{sample_data(data)} in quotation
#' @param covars Either a named list with covariables, OR if \code{data} is a \code{phyloseq} object a character vector with names of the variables in \code{sample_data(data)}
#' @param out.all If TRUE, will run "LRT" which will produce one p-value for the \code{predictor}. If FALSE will run "Wald" test and will output p-value from one level of the predictor specified by \code{coeff}. If NULL (default) set as TRUE for multi-class \code{predictor} and FALSE otherwise
#' @param p.adj Character. P-value adjustment. Default "fdr". See \code{p.adjust} for details
#' @param coeff Integer. The log2FoldChange (and p-value if test="Wald") will be associated with this coefficient. This coefficient is by default compared to the intercept (1. level of \code{predictor}), change this with \code{coeff.ref}. Default 2, i.e. the 2. level of the \code{predictor}.
#' @param coeff.ref Integer. Reference level of the \code{predictor}. Default the intercept, = 1 
#' @param allResults If TRUE will return raw results from the \code{DESeq} function
#' @param ... Additional arguments for the \code{DESeq} function
#' @return A data.frame with with results.
#' @examples
#' # Creating random count_table and predictor
#' set.seed(4)
#' mat <- matrix(rnbinom(200, size = 0.1, mu = 500), nrow = 20, ncol = 10)
#' rownames(mat) <- 1:20
#' pred <- c(rep("Control", 5), rep("Treatment", 5))
#' 
#' # Running masl
#' res <- DA.masl(data = mat, predictor = pred)
#' @export

DA.masl = function (data, predictor, paired = NULL, covars = NULL, out.all = NULL, 
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

    # normalize and transform data
    count_table <- apply(count_table, 2, function(x) {x/max(sum(x), 1e-32)}) # TSS normalization
    count_table <- apply(count_table, 2, function(x) log(x + min(x[x > 0])/2)) # LOG transformation (add pseudocount of min/2)

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
            min_abundance = 0,
            min_prevalence = 0,
            normalization = "NONE",
            transform = "NONE",
            analysis_method = "LM",
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
    
    
    res$Method <- "Maaslin2 (LM (tss + log)) "
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
