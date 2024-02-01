decomp_kb <- function(){}

load('data/sce.rda')
n_components=NULL # number of gene sets + 1 (thums rule)
                  # or a dictionary that maps cell type
                  # to an integer per cell type
# FIX THIS TO WORK WITH A GIVEN NUMBER OF PARAMETERS!!!!
envname='r-decomp'
ignore_warnings=TRUE
use_highly_variable=FALSE
metadata_field='is.doublet'
cell_type_key=metadata_field
result_name='kb'

# nested list of lists gene_set_dict | cell type_i | gene_set_i
# gene sets are not cell specific
# Note that one key in the dictionary must be 'global' with the corresponding
# value being a dictionary of gene sets which apply to all cells

# The cell type labels have to match with the cell type
# labels in the gene set dictionary


# if use_cell_types == False then maps gene set names to gene sets ;
# must contain "global" key in addition to every unique cell type under .obs.
# Note that if you have a cell type in your adata for which
# you do not have any gene sets
# in your gene set annotation dictionary you must include an
# empty dictionary under that cell type key.
prior_knowledge = list('0' = list('gene_set_1' = list("KLHL17",
                                                      "CCNL2",
                                                      "ATAD3B",
                                                      "NOL9"),
                                  'gene_set_2' = list( "CENPS-CORT",
                                                       "FBXO2",
                                                       "KLHL21")),
                       '1' = list('gene_set_1' = list("NPPA",
                                                      "FHAD1",
                                                      "FBXO42",
                                                      "ATP13A2"),
                                  'gene_set_3' = list("PADI2",
                                                      "AHDC1",
                                                      "PPP1R8",
                                                      "LAPTM5")),
                       # set of genes which will be evauated for all the
                       # cells in data…
                       'global' = list('gene_set_4' = list("AZIN2",
                                                           "ZSCAN20",
                                                           "ZC3H12A")))

n_markers = 50
num_epochs = 2
min_gs_num = 0
gene_set_confidence_assignment = 0.2 # fraction of overlapping genes
use_cell_types = TRUE
use_weights = TRUE # edge weights are estimated based on graph
# structure and used throughout training
lam = 0.1 # trust in the prior (the higher, the higher the trust in it!)
delta = 0.001 # scaling to make the values comparable
# varies depending on data and gene sets, try between 0.5 and 0.001



message('--- Checking Packages ---')
is_package_installed('reticulate')
is_package_installed('sceasy')     # moving from and to sce and anndata

packages.vec = c('scSpectra', 'pandas')

is_python_package_installed(envname      = envname,
                            packages.vec = packages.vec)

spectra <- reticulate::import('Spectra')
pd      <- reticulate::import('pandas')

if(ignore_warnings){
    warnings <- reticulate::import('warnings')
    warnings$filterwarnings("ignore")
}

if(!is.null(n_components)){
    if(is.character(n_components)){
        n_components <- as.integer(n_components)
    }
}
n_markers = as.integer(n_markers)
num_epochs = as.integer(num_epochs)

if(min_gs_num == 0){
    min_gs_num <- 1
    warning('> Forcing the min_gs_num to 1; must be a positive integer <')
}
min_gs_num = as.integer(min_gs_num)


colData(sce)[[metadata_field]] <- as.character(colData(sce)[[metadata_field]])


# converting SCE into AnnData (with sparse matrix)
adata <- sce2adata_sparse(sce = sce, envname = envname)

# this could be a function on its own
#####
# to match correctly the metadata_field with prior knowledge we need to
# transform into characters the chosen filed for the metadata_field

# S is for string
if(as.logical(toupper(adata$obs[[metadata_field]]$dtype != 'S'))){
    adata$obs[[metadata_field]] <- adata$obs[[metadata_field]]$astype('int')
    adata$obs[[metadata_field]] <- adata$obs[[metadata_field]]$astype('string')
}

# checking if prior_knowledge is found in the data
check_gene_set_dictionary <- spectra$Spectra_util$check_gene_set_dictionary
prior_knowledge <- reticulate::r_to_py(prior_knowledge)
prior_knowledge = check_gene_set_dictionary(adata       = adata,
                                            annotations = prior_knowledge,
                                            obs_key     = metadata_field,
                                            global_key  = 'global')

#####



message('--- Running Knowledge Based Decomposition ---')
model = spectra$est_spectra(adata               = adata,
                            L                   = n_components,
                            gene_set_dictionary = prior_knowledge,
                            use_highly_variable = use_highly_variable,
                            use_cell_types      = use_cell_types,
                            cell_type_key       = cell_type_key,
                            use_weights         = use_weights,
                            lam                 = lam,
                            delta               = delta,
                            n_top_vals          = n_markers,
                            overlap_threshold   = gene_set_confidence_assignment,
                            min_gs_num          = min_gs_num,
                            num_epochs          = num_epochs)



message('--- Storing Results ---')
# General Parameters
result_name <- change_default_name(result_name, reducedDimNames(sce))
# pattern names will be defined later…

# Gene View
# Spectra factors the patterns of the W matrix
kb_h.dense <- t(as.matrix(adata$uns['SPECTRA_factors']))
kb_h.dgCmatrix <- as(kb_h.dense, 'sparseMatrix')

# Setting the number of patterns if they were selected by default
if(is.null(n_components)){
    n_components <- ncol(kb_h.dgCmatrix)
}
patterns_names <- paste0(rep('KB_', n_components), seq_len(n_components))

colnames(kb_h.dgCmatrix) <- patterns_names
rownames(kb_h.dgCmatrix) <- rownames(sce)
metadata(sce)[[result_name]] <- nmf_w.dgCmatrix

# Cell View







#explore eta parameter to detect new factors
model$return_eta_diag()


