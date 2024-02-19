#' Computes Cell Capture Score for a given annotation vector
#'
#' This function computes a cell capture score for a given annotation by
#' testing with a Welch's t-test the difference between the values of
#' a pattern associated to a given annotation against all the others.
#' Hence, it gives an idea of how much of the annotations information
#' is present in that given pattern, compared to all the others.
#'
#' @param sce <SingleCellExperiment object> SCE object with at least one
#' reducedDim representation saved.
#' @param decomp <character> specifying the reducedDim representation to use
#' @param alternative <character> specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less".
#' You can specify just the initial letter.
#' @param adj_pvalue_method <character> correction method, a character string.
#' Can be abbreviated. The adjustment methods include the Bonferroni correction
#' ("bonferroni") in which the p-values are multiplied by the number of
#' comparisons. Less conservative corrections are also included by Holm
#' (1979) ("holm"), Hochberg (1988) ("hochberg"), Hommel (1988) ("hommel"),
#'Benjamini & Hochberg (1995) ("BH" or its alias "fdr"), and Benjamini &
#' Yekutieli (2001) ("BY"), respectively. A pass-through option ("none")
#' is also included. The set of methods are contained in the p.adjust.methods
#' vector for the benefit of methods that need to have the method as an option
#' and pass it on to p.adjust. (from the p.adjusted function)
#' @param adj_pvalue <bool> either TRUE or FALSE, to specify whether to
#' correct the pvalue for multiple testing or not, respectively.
#' @param min_log_10 <bool> whether to compute the -log10(p.value).
#' @param cap_multiplier <double or integer> number which multiplies the maximum
#' non infinite pvalue to impute infinite values.
#' @param remove_infs <bool> whether to impute Infinite values by the
#' max(pvalue) * cap_value + 10
#' @return matrix metadata x prototypes with cell capture scores
#' @examples
#' cell_capture(sce, sce$cluster, 'pca')
#' @export
cell_capture <- function(sce,
                         annotations,
                         decomp,
                         alternative='greater',
                         adj_pvalue_method='BH',
                         adj_pvalue=TRUE,
                         min_log_10=TRUE,
                         cap_multiplier=2,
                         remove_infs=FALSE){

  ### Description ###
  # Computing the cell capture score
  # via a Welch's t-test considering
  # the cells annotated to a given Identity
  # and comparing one prototype / partition
  # to all the others.

  # in case anything brakes… Ask Carlo Manenti!

  ## check for normality of the data!


  is_package_installed('SingleCellExperiment')


  sorted_annotation <- sort(unique(annotations))

  #sorted_partition <- sort(unique(partitions))

  # getting the H matrix for a given decomposition level
  H <- SingleCellExperiment::reducedDim(sce, decomp)
  sorted_partition <- seq(ncol(H))

  # getting the patterns for each annotations
  #H_split <- split_sparse(H, annotations, by_col = T, sorted = T)

  # computing
  res <- sapply(sorted_annotation, # looping over annotations
                function(annotation_i){
                  sapply(sorted_partition, # looping over partitions
                         function(partition_i){

                           if(length(H[annotations == annotation_i , partition_i]) == 0) return(10e7)
                           res_i <- t.test(H[annotations == annotation_i , partition_i],
                                           H[annotations == annotation_i , -partition_i],
                                           alternative = alternative)


                           # extracting the p value
                           return(res_i$p.value)})
                }) # closing the first sapply

  # Do we want to adjust for multiple testing?
  if(adj_pvalue){
    dimensions <- dim(res)
    adj_pvalues <- p.adjust(as.vector(res), method = adj_pvalue_method)
    res <- matrix(adj_pvalues, nrow = dimensions[1], ncol = dimensions[2])

  }

  # Do we want to compute the -log10 p value?
  if(min_log_10){
    # computing the -log10
    res <- -log10(res)
    # capping the Inf values
    if(remove_infs){
    res[is.infinite(res)] <- max(res[!is.infinite(res)]) + 10 * cap_multiplier
  }}

  # giving names
  colnames(res) <- sorted_annotation
  rownames(res) <- sorted_partition
  # transposing to have identity x prototypes
  res <- t(res)

  return(res)}

