partition_metacell <- function(sce,
                               target_number_of_metacells=100,
                               envname='r-metacells',
                               assay_name='full_run',
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
                               num_parallel_piles=NULL,
                               return_model=FALSE,
                               seed=42,
                               result_name='metacells',
                               ...){

    ### Description ###
    # Computes Metacells from a dense or dgCMatrix! (NOT from dgTMatrix)
    # Many of the parameters are task-specific and the results of the method
    # should be evaluated carefully!
    # In case, ask Carlo

    # To do : • remove loggings for all the functions
    #         • flag to keep all the genes after the definition of metacells

    # example usage
    #load('data/sce.rda')
    # partition_metacell(sce)


    message('--- Checking packages ---')
    is_python_package_installed(packages.vec = 'metacells', envname = envname)
    #.rs.restartR() in case there are problems with the loading metacells
    # enforcing the use of the correct environment
    reticulate::use_virtualenv(envname)

    mc <- reticulate::import('metacells')
    ad <- reticulate::import('anndata')
    np <- reticulate::import('numpy')

    # enforcing the variables type to avoid crashes in python
    if(!is.null(min_umi)) min_umi <- as.integer(min_umi)
    if(!is.null(max_umi)) max_umi <- as.integer(max_umi)
    if(!is.null(max_excluded_gene_fraction)){
      max_excluded_gene_fraction <- as.double(max_excluded_gene_fraction)}
    seed <-  as.integer(seed)

    ## IMPORTANT
    # float32 type required by metacells will be deprecated in late 2024
    matrix.dgCMatrix <- counts(sce)
    data.h5ad <- ad$AnnData(t(matrix.dgCMatrix), dtype = 'float32')
    data.h5ad$obs_names <- colnames(matrix.dgCMatrix)
    data.h5ad$var_names <- rownames(matrix.dgCMatrix)
    # enforcing unique feature names
    data.h5ad$var_names_make_unique()

    mc$ut$top_level(data.h5ad)
    mc$ut$set_name(data.h5ad, name = assay_name)



    if(quality_filters) message('--- Quality Filtering Steps ---')

    if(!quality_filters){
      min_umi <- NULL
      max_umi <- NULL
      excluded_gene_names <- NULL
      excluded_gene_patterns <- NULL
      lateral_gene_name <- NULL
      lateral_gene_patterns <- NULL
      noisy_gene_names <- NULL
      noisy_gene_patterns <- NULL
      max_excluded_gene_fraction <- NULL
    }

    if(is.null(min_umi) & is.null(max_umi)){
      message('note : Avoiding filtering cells by total UMI counts')
    }

    if(is.null(excluded_gene_names) & is.null(excluded_gene_patterns)){
      message('note : Considering all gene names and patterns')
    }

    # it also excludes highly variant genes which are also uncorrelated
    # with any other gene
    if(quality_filters) message('--- 1. excluding genes ---')
    mc$pl$exclude_genes(data.h5ad,
                        excluded_gene_names = excluded_gene_names,
                        excluded_gene_patterns = excluded_gene_patterns,
                        random_seed = seed);

    if(quality_filters) message('--- 2. excluding cells ----')
    mc$tl$compute_excluded_gene_umis(data.h5ad) # for each cell
    mc$pl$exclude_cells(data.h5ad,
                        properly_sampled_min_cell_total = min_umi,
                        properly_sampled_max_cell_total = max_umi,
                        properly_sampled_max_excluded_genes_fraction = max_excluded_gene_fraction)

    if(quality_filters) message('--- 3. cleaning the data ---')
    data.h5ad <- mc$pl$extract_clean_data(data.h5ad, name = paste0(assay_name, '_clean'))


    if(quality_filters) message('--- 4. marking lateral genes ---')
    mc$pl$mark_lateral_genes(data.h5ad,
                             lateral_gene_names = lateral_gene_name,
                             lateral_gene_patterns = lateral_gene_patterns,
                             op = 'set') # op = 'add' adds gene names

    if(quality_filters) message('--- 5. marking noisy genes ---')
    mc$pl$mark_noisy_genes(data.h5ad,
                           noisy_gene_names = noisy_gene_names,
                           noisy_gene_patterns = noisy_gene_patterns,
                           op = 'set')

    if(is.null(num_parallel_piles)){
      max_parallel_piles <- mc$pl$guess_max_parallel_piles(data.h5ad)
      mc$pl$set_max_parallel_piles(max_parallel_piles)
    }else{mc$pl$set_max_parallel_piles(num_parallel_piles)
    }

    if(is.null(target_number_of_metacells) | target_number_of_metacells == 0){
      target_metacell_size <- as.integer(100)
    }else{
    metacell2cell_ratio <- ncol(matrix.dgCMatrix) / target_number_of_metacells
    target_metacell_size <- as.integer(ceiling(metacell2cell_ratio))
    }



    message('--- Computing MetaCells ---')
    mc$pl$divide_and_conquer_pipeline(data.h5ad,
                                      target_metacell_size = target_metacell_size,
                                      random_seed = seed,
                                      ...)

    metacells.h5ad <- mc$pl$collect_metacells(data.h5ad, random_seed = seed)

    metacells_w.dgCMatrix <- t(metacells.h5ad$X)
    rownames(metacells_w.dgCMatrix) <- np$array(data.h5ad$var_names)
    colnames(metacells_w.dgCMatrix) <- paste0(rep('MetaCell_'),
                                         seq(1, ncol(metacells_w.dgCMatrix)))

    # also introducing missing genes, with null contribution, in the w matrix
    missing_genes <- sapply(unique(data.h5ad$obs$metacell), function(metacell_i){
        #cells <- which(data.h5ad$obs$metacell == metacell_i)
        missing <- !(rownames(matrix.dgCMatrix) %in% rownames(metacells_w.dgCMatrix))
        missing_genes <- rownames(matrix.dgCMatrix)[missing]
        # this might be useful down the road;
        # but not now since we are interested in how it behaves metacells…
        #rowSums(matrix.dgCMatrix[missing_genes, cells])
        enforced_zeros <- rep(0, length(missing_genes))
        names(enforced_zeros) <- missing_genes
        enforced_zeros
    })

    metacells_w.dgCMatrix <- rbind(metacells_w.dgCMatrix, missing_genes)
    result_name <- change_default_name(result_name, names(metadata(sce)))
    metadata(sce)[[result_name]] <- metacells_w.dgCMatrix

    if(return_model){
      return(list(model = metacells.h5ad, obj =  sce))
    }else{return(sce)}
}










