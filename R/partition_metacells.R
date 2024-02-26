#' Partition Approach : MetaCells
#'
#' This function performs MetaCell analysis on a SingleCellExperiment object
#' and stores the results directly in it. This function leverages the python
#' implementation of MetaCell to achieve higher performance.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param target_number_of_metacells <integer or NULL> default 100;
#' Number of approximate number of MetaCells (latents) to be defined.
#' If NULL or 0, it will fallback to use approximately 100 cells for each
#' MetaCells.
#' @param num_parallel_piles <NULL or integer> Number of parallel processes.
#' If set to NULL, the method will estimate the appropriate number of
#' parallel processes that the machine can run.
#' @param quality_filters <bool> default TRUE; Whether to perform quality
#' filtering steps prior to the estimation of MetaCells. It is strongly
#' suggested to perform appropriate quality filtering prior to the estimation
#' of MetaCells.
#' @param min_umi <integer> default 800; Minimum number of Unique Molecular
#' Identifier (UMI) to consider a cell.
#' @param max_umi <integer> default 20000; Maximum number of Unique Molecular
#' Identifier (UMI) to consider a cell.
#' @param excluded_gene_names <list of character> default c('XIST, MALAT1');
#' List of blacklisted genes. These genes  should not be considered during the
#' analysis. example : XIST is related to inactivation of the X chromosome;
#' MALAT1 is associated with metastatic processes.
#' @param excluded_gene_patterns <list of characters> default c('MT-.*');
#' List of blacklisted gene patterns, which should not be considered during the
#' analysis. example : MT- (all the genes starting with MT-; basically all
#' the mitochondrial genes).
#' @param max_excluded_gene_fraction <double> default 0.25; maximum fraction
#' of genes excluded for a given cell. If a cell has more than 25% of its
#' expressed genes excluded, then is not considered for the following analysis.
#' @param lateral_gene_name <list of characters> default c("AURKA", "MCM3",
#' "MCM4", "MCM7", "MKI67", "PCNA", "RRM2", "SMC4", "TPX2", "FOS",
#' "HSP90AB1", "TXN"); (Cell-Cycle and Stress related genes)
#' List of genes which are not considered during the
#' community detection step, thus they do not contribute into defining
#' MetaCells, even though they are still used to build the KNN graph.
#' @param lateral_gene_patterns <list of characters> default c("RP[LS].*");
#' (Ribosomial Genes). List of gene patterns which are not considered during
#' the community detection step, thus they do not contribute into defining
#' MetaCells, even though they are still used to build the KNN graph.
#' @param noisy_gene_names <list of characters> default c();
#' List of genes which are expected to be highly variable between cells without
#' bearing any information in their variability. This genes will not be
#' considered Highly Variable Genes even if they are Highly Variable.
#' @param noisy_gene_patterns <list of characters> default c();
#' List of gene patterns which are expected to be highly variable between
#' cells without bearing any information in their variability. This genes
#' will not be considered Highly Variable Genes even if they are Highly
#' Variable.
#' @param assay_name <character> default 'full_run';
#' Name assigned to the metacell assay inside the AnnData Object.
#' It is useful when you want to retrieve also the full model, thus enabling you
#' to set a specific name for the results of the MetaCell run.
#' @param result_name <character> default 'metacells';
#' Name used to store the result in the SingleCellExperiment object.
#' @param ignore_warnings <bool> default TRUE; Whether to ignore warnigns
#' @param envname <character> default 'r-decomp';
#' Specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only
#' the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @param ... <extra arguments for the divide_and_conquer_pipeline function>
#' @return either a SingleCellExperiment object with W matrix
#' representation or the SingleCellExperiment object and the model
#' used to perform MetaCell analysis.
#' @examples
#' #partition_metacells(sce, target_number_of_metacells = 2, min_umi = 5)
#' @export
partition_metacells <- function(sce,
                               target_number_of_metacells=100,
                               num_parallel_piles=NULL,
                               quality_filters=TRUE,
                               min_umi=800,
                               max_umi=20000,
                               excluded_gene_names=c('XIST, MALAT1'),
                               excluded_gene_patterns=c('MT-.*'),
                               max_excluded_gene_fraction=0.25,
                               lateral_gene_name=c(# Cell-cycle
                                 "AURKA", "MCM3", "MCM4",
                                 "MCM7", "MKI67", "PCNA",
                                 "RRM2", "SMC4", "TPX2",
                                 # Stress
                                 "FOS", "HSP90AB1", "TXN"),
                               lateral_gene_patterns=c("RP[LS].*"), # Ribosomal
                               noisy_gene_names=c(),
                               noisy_gene_patterns=c(),
                               assay_name='full_run',
                               ignore_warnings=TRUE,
                               result_name='metacells',
                               envname='r-partition',
                               return_model=FALSE,
                               seed=42,
                               verbose=FALSE,
                               ...){

    ### Description ###
    # Computes Metacells from a dense or dgCMatrix! (NOT from dgTMatrix)
    # Many of the parameters are task-specific and the results of the method
    # should be evaluated carefully!
    # In case, ask Carlo

    # To do : • remove loggings for all the functions

    # example usage
    #load('data/sce.rda')
    # partition_metacells(sce)

    if(verbose){
        message('--- Checking packages ---')
    }
    is_package_installed('reticulate')
    is_python_package_installed(packages.vec = 'metacells', envname = envname)
    #.rs.restartR() in case there are problems with the loading metacells
    # enforcing the use of the correct environment

    # avoiding silly warning (where possibile)
    ignore_warnings(ignore_warnings = ignore_warnings,
                    envname = envname,
                    verbose = verbose)

    # loading the python packages
    mc <- reticulate::import('metacells', delay_load = TRUE) # performing MetaCell analysis
    ad <- reticulate::import('anndata', delay_load = TRUE)   # from SCE to AnnData
    np <- reticulate::import('numpy', delay_load = TRUE)     # data handling

    # enforcing the variables type to avoid crashes in python
    if(!is.null(min_umi)){
        min_umi <- as.integer(min_umi)
    }
    if(!is.null(max_umi)){
        max_umi <- as.integer(max_umi)
    }
    if(!is.null(max_excluded_gene_fraction)){
        max_excluded_gene_fraction <- as.double(max_excluded_gene_fraction)
    }
    seed <- as.integer(seed)

    # IMPORTANT
    # float32 type required by MetaCell will be deprecated in late 2024
    matrix.dgCMatrix <- SingleCellExperiment::counts(sce)
    data.h5ad <- ad$AnnData(Matrix::t(matrix.dgCMatrix), dtype = 'float32')
    data.h5ad$obs_names <- colnames(matrix.dgCMatrix)
    data.h5ad$var_names <- rownames(matrix.dgCMatrix)
    # enforcing unique feature names
    data.h5ad$var_names_make_unique()

    mc$ut$top_level(data.h5ad)
    mc$ut$set_name(data.h5ad, name = assay_name)



    if(quality_filters & verbose){
        message('--- Quality Filtering Steps ---')
    }
    # overwriting all the quality filters
    if(!quality_filters){
        min_umi                    <- NULL
        max_umi                    <- NULL
        excluded_gene_names        <- NULL
        excluded_gene_patterns     <- NULL
        lateral_gene_name          <- NULL
        lateral_gene_patterns      <- NULL
        noisy_gene_names           <- NULL
        noisy_gene_patterns        <- NULL
        max_excluded_gene_fraction <- NULL
    }

    # General messages to acknowledge the complete lack of quality filters
    if(is.null(min_umi) & is.null(max_umi) & verbose){
        message('note : Avoiding filtering cells by total UMI counts')
    }

    if(is.null(excluded_gene_names) & is.null(excluded_gene_patterns) & verbose){
        message('note : Considering all gene names and patterns')
    }

    # it also excludes highly variant genes which are also uncorrelated
    # with any other gene.
    if(quality_filters & verbose){
        message('--- 1. excluding genes ---')
    }
    mc$pl$exclude_genes(data.h5ad,
                        excluded_gene_names = excluded_gene_names,
                        excluded_gene_patterns = excluded_gene_patterns,
                        random_seed = seed);

    if(quality_filters & verbose){
        message('--- 2. excluding cells ----')
    }
    mc$tl$compute_excluded_gene_umis(data.h5ad) # for each cell
    mc$pl$exclude_cells(data.h5ad,
                        properly_sampled_min_cell_total = min_umi,
                        properly_sampled_max_cell_total = max_umi,
                        properly_sampled_max_excluded_genes_fraction = max_excluded_gene_fraction)

    if(quality_filters & verbose){
        message('--- 3. cleaning the data ---')
    }
    data.h5ad <- mc$pl$extract_clean_data(data.h5ad, name = paste0(assay_name, '_clean'))


    if(quality_filters & verbose){
        message('--- 4. marking lateral genes ---')
    }
    mc$pl$mark_lateral_genes(data.h5ad,
                             lateral_gene_names = lateral_gene_name,
                             lateral_gene_patterns = lateral_gene_patterns,
                             op = 'set') # op = 'add' adds gene names

    if(quality_filters & verbose){
        message('--- 5. marking noisy genes ---')
    }
    mc$pl$mark_noisy_genes(data.h5ad,
                           noisy_gene_names = noisy_gene_names,
                           noisy_gene_patterns = noisy_gene_patterns,
                           op = 'set')

    # setting the proper number of parallel processes
    if(is.null(num_parallel_piles)){
        max_parallel_piles <- mc$pl$guess_max_parallel_piles(data.h5ad)
        mc$pl$set_max_parallel_piles(max_parallel_piles)
    }else{
        mc$pl$set_max_parallel_piles(num_parallel_piles)
    }

    # setting the correct number of latents (MetaCells)
    if(is.null(target_number_of_metacells) | target_number_of_metacells == 0){
        target_metacell_size <- as.integer(100)
    }else{
        metacell2cell_ratio <- ncol(matrix.dgCMatrix) / target_number_of_metacells
        target_metacell_size <- as.integer(ceiling(metacell2cell_ratio))
    }



    if(verbose){
        message('--- Computing MetaCells ---')
    }
    mc$pl$divide_and_conquer_pipeline(data.h5ad,
                                      target_metacell_size = target_metacell_size,
                                      random_seed = seed,
                                      ...)



    if(verbose){
      message('--- Storing Results ---')
    }

    # gathering MetaCells
    metacells.h5ad <- mc$pl$collect_metacells(data.h5ad, random_seed = seed)

    # gene view
    metacells_w.dgCMatrix <- Matrix::t(metacells.h5ad$X)
    rownames(metacells_w.dgCMatrix) <- np$array(data.h5ad$var_names)
    colnames(metacells_w.dgCMatrix) <- paste0(rep('MetaCell_'),
                                         seq(1, ncol(metacells_w.dgCMatrix)))

    # also introducing missing genes, with null contribution, in the w matrix
    missing_genes <- sapply(unique(data.h5ad$obs$metacell), function(metacell_i){
          # cells <- which(data.h5ad$obs$metacell == metacell_i)
          missing <- !(rownames(matrix.dgCMatrix) %in% rownames(metacells_w.dgCMatrix))
          missing_genes <- rownames(matrix.dgCMatrix)[missing]
          # this might be useful down the road;
          # but not now since we are interested in how it behaves metacells…
          # rowSums(matrix.dgCMatrix[missing_genes, cells])
          enforced_zeros <- rep(0, length(missing_genes))
          names(enforced_zeros) <- missing_genes
          enforced_zeros
    })

    metacells_w.dgCMatrix <- rbind(metacells_w.dgCMatrix, missing_genes)
    result_name <- change_default_name(result_name, names(S4Vectors::metadata(sce)))
    S4Vectors::metadata(sce)[[result_name]] <- metacells_w.dgCMatrix



    # this is the return
    return_model(sce = sce, model = metacells.h5ad, return_model = return_model)
}










