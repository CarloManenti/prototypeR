partition_recursive<- function(sce,
                               result_name='recursive_clustering',
                               return_model=FALSE,
                               seed=42,
                               ## Main parameters
                               mc.cores=10, # number of parallel processes (RAM depends on the data used)
                               q.diff.th=0.7, # absolute difference of expression.(discrete) 1-0 (continous)
                               q1.th=0.4, # Background set : proportion of cells with expression > low.th
                               min.cells=10,  # minimal number of cells of a cluster
                               de.score.th=150, # 150 for more than 1000 cells | 40 for less than 1000
                               ### Selecting DEGs (with limma)
                               padj.th=0.05,  # adjusted p value threshold for DE genes
                               lfc.th=1,     # log2 fold change threshold for DE genes.
                               ## Splitting or merging clusters
                               low.th=1, # minimum log2 expression to determine if a gene is detected
                               q2.th=NULL, # Foreground set : proportion of cells with expression > low.th
                               ## Consensus parameters
                               niter=100, # number of iterations to define the consensus
                               dim.method="pca", # type of clustering
                               output_dir="subsample_PCA", # output directory to store intermediate results
                               single_run=TRUE){
    ### Description ###
    # Runs Recursive Clustering as implemented by the Allen Brain Institute
    # via scrattch::hicat

    # note: It would be best to re-implement it via scrattch::bigcat!
    #       It would be significantly faster on big data sets.

    # example usage
    #sce <- partition_recursive(sce)


    message('--- Checking packages ---')
    is_package_installed('scrattch.bigcat')
    is_package_installed('foreach')
    is_package_installed('doParallel')
    is_package_installed('iterators')
    is_package_installed('parallel')




    message('--- Pre-processing ---')
    counts.dgCMatrix <- counts(sce)

    # might be quite slow!
    cpm.dgCMatrix <- scrattch.hicat::cpm(counts.dgCMatrix)
    logcounts.dgCMatrix <- as(log2(cpm.dgCMatrix + 1), 'sparseMatrix')

    # make the matrix non dgC


    message('--- Consensus Definition ---')
    de.param <- scrattch.hicat::de_param(padj.th     = padj.th,
                                         lfc.th      = lfc.th,
                                         low.th      = low.th,
                                         q1.th       = q1.th,
                                         q2.th       = q2.th,
                                         q.diff.th   = q.diff.th,
                                         de.score.th = de.score.th,
                                         min.cells   = min.cells)

    ## WARNING
    # it uses plenty of RAM even with few processes.
    # Working locally may not be a good idea.
    # Reserve a node on a HPC computer with plenty of RAM and processors.
    recursive.model <- scrattch.hicat::run_consensus_clust(logcounts.dgCMatrix,
                                                           niter      = niter,
                                                           de.param   = de.param,
                                                           dim.method = dim.method,
                                                           # dim.method to be removed
                                                           # when using bigcat
                                                           output_dir = output_dir,
                                                           mc.cores   = mc.cores)



    message('--- Storing results ---')
    result_name <- change_default_name(result_name, colnames(colData(sce)))
    names(recursive.model[['cl.result']])[1] <- result_name
    recursive_clustering <- recursive.model[['cl.result']][1]
    colData(sce) <- cbind(colData(sce), result_name = recursive_clustering)

    if(!return_model) return(sce)
    return(list(model = result, obj = sce))
}




