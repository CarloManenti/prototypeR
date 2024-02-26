#' Partition Approach : Recursive Clustering
#'
#' This function performs Recursive Clustering as implemented via the
#' scrattch.hicat package of the Allen Brain Atlas. The core idea is to define
#' sub communities which show a binary gene expression pattern.
#' The analysis has X key steps:
#' - Iterative Clustering
#' - High Variance Gene Selection
#' - Dimensionality reduction using PCA or WGCNA
#' - Jaccard-Louvain or Hierarchical (Ward) Clustering
#' - Cluster Merging Based on Presence of Differentially Expressed Genes.
#' These steps are repeated for each community, till it is no longer possible
#' differentiate sub-communities.
#' Furthermore, all the process is iterated multiple times to ensure robustness
#' via bootstrapping on %80 of all available data.
#' -----------------------------------------------------------------------------
#' Please refer to [What is hicat](https://taxonomy.shinyapps.io/scrattch_tutorial/)
#' for a more comprehensive explanation
#' -----------------------------------------------------------------------------
#' Single Cell RNA-seq Analysis for Transcriptomic Type CHaracterization -
#' Hierarchical, Iterative Clustering for Analysis of Transcriptomics
#' (SCRATCH - HICAT)
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param de.score.th <integer> default 150;
#' Differential Expression Score (DES) to state that a cluster can be split into
#' two smaller ones. The DES is computed based on the significance of the
#' differential expression. Each gene can contribute at most 20 'points'
#' to the score, regardless of how significant it is.
#' For more than 1.000 cells the suggested value is 150.
#' For less than 40 cells the suggested value is 40.
#' Adjusting this score will result more or less communities; but it is also
#' important to set the other parameters accordingly.
#' @param q.diff.th <double, from 0 to 1> default 0.7; Fraction of
#' cells with differentially expressed genes between communities
#' it is defined as abs(q1.th - q2.th)/max(q1.th, q2.th) and
#' should be greater than q.diff.th. to split a community into two.
#' Lower values will fit best a continuous of cell states-types, while higher
#' values will lead to desecrate communities.
#' @param mc.cores <integer or NULL> default 2; Number of cores to be used for
#' parallelization. If 0 or NULL, it will fallback to estimating by default the
#' number of cores to use, basically all the available cores - 2.
#' @param n_iter <integer> default 100; Number of iterations used to define the
#' final consensus.
#' @param min.cells <integer> default 10; Minimum number of cells to define
#' a community.
#' @param low.th <double> default 1; Threshold of logFC to define a gene as
#' detected in a given community.
#' @param q1.th <double, from 0 to 1 or NULL> default 0.4; Fraction of
#' cells in the  Background set with an expression > low.th.
#' @param q2.th <double, from 0 to 1 or NULL> default 0.4; Fraction of
#' cells in the Foreground set with an expression > low.th.
#' @param padj.th <double> default 0.01; Adjusted p value threshold to
#' consider a gene differentially expressed.
#' @param lfc.th <double> default 1; Log Fold Change threshold to
#' consider a gene differentially expressed.
#' @param dim.method <character> default pca; Type of dimensionality reduction
#' used in the analysis. Can be on of:
#' *pca*    : Principal Component Analysis, suggested for medium to large size
#' data sets.
#' *WGCNA* : Weighted Correlation Gene Network Analysis, suggested for small
#' size data sets with high amount of genes detected (more than 5000).
#' @param output_dir <character> default subsample_PCA;
#' Directory used to store intermediate results of each single run.
#' If the processes is interrupted, it will still be able to use the computed
#' runs and re-start from there.
#' note : to re-run completely the method one must change the output_dir or
#' delete complete it.
#' @param result_name <character> default 'recursive_clustering';
#' Name used to store the result in the SingleCellExperiment object.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only
#' the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @return either a SingleCellExperiment object with a membership vector
#' representations or the SingleCellExperiment object and the model
#' used to perform Recursive Clustering
#' @examples
#' # it can not be run on such small data set
#' @export
partition_recursive<- function(sce,
                               de.score.th=150,
                               q.diff.th=0.7,
                               mc.cores=2,
                               n_iter=100,
                               min.cells=10,
                               low.th=1,
                               q1.th=0.4,
                               q2.th=NULL,
                               padj.th=0.01,
                               lfc.th=1,
                               dim.method="pca",
                               output_dir="subsample_PCA",
                               result_name='recursive_clustering',
                               return_model=FALSE,
                               seed=42,
                               verbose=FALSE
                               ){
    ### Description ###
    # Runs Recursive Clustering as implemented by the Allen Brain Institute
    # via scrattch::hicat

    # note: It would be best to re-implement it via scrattch::bigcat!
    #       It would be significantly faster on big data sets.


    if(verbose){
        message('--- Checking packages ---')
    }
    is_package_installed('scrattch.hicat')
    is_package_installed('foreach')
    is_package_installed('doParallel')
    is_package_installed('iterators')
    is_package_installed('parallel')
    # setting the number of parallel processes by default
    if(is.null(mc.cores) | mc.cores == 0){
        mc.cores <- parallel::detectCores() - 2
        if(mc.cores < 1){
            mc.cores <- 1
        }
        if(verbose){
            message(paste0('note: ',mc.cores ,' cores used'))
        }
    }



    if(verbose){
        message('--- Pre-processing ---')
    }
    counts.dgCMatrix <- SingleCellExperiment::counts(sce)
    # might be quite slow!
    cpm.dgCMatrix <- scrattch.hicat::cpm(counts.dgCMatrix)
    logcounts.dgCMatrix <- methods::as(log2(cpm.dgCMatrix + 1), 'sparseMatrix')



    if(verbose){
        message('--- Consensus Definition ---')
    }
    de.param <- scrattch.hicat::de_param(padj.th     = padj.th,
                                         lfc.th      = lfc.th,
                                         low.th      = low.th,
                                         q1.th       = q1.th,
                                         q2.th       = q2.th,
                                         q.diff.th   = q.diff.th,
                                         de.score.th = de.score.th,
                                         min.cells   = min.cells)
    # WARNING
    # it uses plenty of RAM even with few processes.
    # Working locally may not be a good idea.
    # Reserve a node on a HPC computer with plenty of RAM and processors.
    recursive.model <- scrattch.hicat::run_consensus_clust(logcounts.dgCMatrix,
                                                           niter      = n_iter,
                                                           de.param   = de.param,
                                                           dim.method = dim.method,
                                                           # dim.method to be removed
                                                           # when using bigcat
                                                           output_dir = output_dir,
                                                           mc.cores   = mc.cores)



    if(verbose){
        message('--- Storing results ---')
    }
    result_name <- change_default_name(result_name, colnames(SummarizedExperiment::colData(sce)))
    names(recursive.model[['cl.result']])[1] <- result_name
    recursive_clustering <- recursive.model[['cl.result']][1]
    SummarizedExperiment::colData(sce) <- cbind(SummarizedExperiment::colData(sce), result_name = recursive_clustering)



    # this is the return
    return_model(sce = sce, model = recursive.model, return_model = return_model)
}




