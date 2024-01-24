cell_capture <- function(ace,
                         annotations,
                         #partitions,
                         decomp,
                         alternative = 'greater',
                         adj_pvalue_method = 'BH',
                         adj_pvalue = T,
                         min_log_10 = T,
                         cap_multiplier = 2,
                         remove_infs = F){

  ### Description ###
  # Computing the cell capture score
  # via a Welch's t-test considering
  # the cells annotated to a given Identity
  # and comparing one prototype / partition
  # to all the others.

  # in case anything brakes… Ask Carlo Manenti!

  ## check for normality of the data!



  sorted_annotation <- sort(unique(annotations))

  #sorted_partition <- sort(unique(partitions))

  # getting the H matrix for a given decomposition level
  H <- colMaps(ace)[[decomp]]
  sorted_partition <- seq(ncol(H))

  # getting the patterns for each annotations
  #H_split <- split_sparse(H, annotations, by_col = T, sorted = T)

  # computing
  res <- sapply(sorted_annotation, # looping over annotations
                function(annotation_i){
                  sapply(sorted_partition, # looping over partitions
                         function(partition_i){

                           # Welch’s paired t-test
                           # assumption : normal data (to be tested)

                           # if you do not have any values, than is not specific!
                           #if(length(H_split[[annotation_i]][, partition_i]) == 0) return(10e7)

                           # testing the values for a given prototype
                           #res_i <- t.test(H_split[[annotation_i]][, partition_i],
                                           # against all others.
                           #               H_split[[annotation_i]][, -partition_i],
                           #                alternative = alternative)


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

