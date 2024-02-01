decomp_kb <- function(){}

load('data/sce.rda')
n_components=NULL # number of gene sets + 1 (thums rule)
envname='r-decomp'
# nested list of lists gene_set_dict | cell type_i | gene_set_i
# gene sets are not cell specific
# Note that one key in the dictionary must be 'global' with the corresponding
# value being a dictionary of gene sets which apply to all cells

# The cell type labels have to match with the cell type labels in the gene set dictionary

# Note that if you have a cell type in your adata for which you do not have any gene sets
# in your gene set annotation dictionary you must include an empty dictionary under that cell type key.
prior_knowledge = list('0' = list('gene_set_1' = list("ZC3H12A", "MEAF6", "SNIP1", "DNALI1"),
                                  'gene_set_2' = list( "GNL2", "RSPO1", "C1orf109")),
                       '1' = list('gene_set_1' = list("ZC3H12A", "MEAF6", "SNIP1", "DNALI1"),
                                  'gene_set_3' = list("HCRTR1", "PEF1", "COL16A1", "ADGRB2")),
                       'global' = list('gene_set_4' = list("NOL9", "TAS1R1", "ZBTB48")))
metadata = 'is.doublet'




message('--- Checking Packages ---')
is_package_installed('reticulate')
is_package_installed('sceasy')     # moving from and to sce and anndata

packages.vec = c('scSpectra', 'pandas')

is_python_package_installed(envname = envname,
                            packages.vec = packages.vec)

spectra <- reticulate::import('Spectra')
pd <- reticulate::import('pandas')

# Setting the number of components with guidelines…
if(is.null(n_components)){
    class_patterns <- lapply(prior_knowledge, names)
    unique_patterns <- unique(unlist(class_patterns))
    n_components <- length(n_unique_patterns) + 1
    # as suggested by the original paper
}
n_components <- as.integer(n_components)


colData(sce)[[metadata]] <- as.character(colData(sce)[[metadata]])


# converting SCE into AnnData (with sparse matrix)
adata <- sce2adata_sparse(sce = sce, envname = envname)

# this could be a function on its own
#####
# to match correctly the metadata with prior knowledge we need to
# transform into characters the chosen filed for the metadata
if(as.logical(toupper(adata$obs[[metadata]]$dtype != 'S'))){ # is for string
    adata$obs[[metadata]] <- adata$obs[[metadata]]$astype('double')
    adata$obs[[metadata]] <- adata$obs[[metadata]]$astype('string')
}

# checking if prior_knowledge is found in the data
prior_knowledge = spectra$Spectra_util$check_gene_set_dictionary(adata = adata,
                                                                 annotations = prior_knowledge,
                                                                 obs_key     = metadata,
                                                                 global_key  = 'global')
prior_knowledge
#####



message('--- Running Knowledge Based Decomposition ---')


