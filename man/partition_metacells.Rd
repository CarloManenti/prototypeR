% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition_metacells.R
\name{partition_metacells}
\alias{partition_metacells}
\title{Partition Approach : MetaCells}
\usage{
partition_metacells(
  sce,
  target_number_of_metacells = 100,
  num_parallel_piles = NULL,
  quality_filters = TRUE,
  min_umi = 800,
  max_umi = 20000,
  excluded_gene_names = c("XIST, MALAT1"),
  excluded_gene_patterns = c("MT-.*"),
  max_excluded_gene_fraction = 0.25,
  lateral_gene_name = c("AURKA", "MCM3", "MCM4", "MCM7", "MKI67", "PCNA", "RRM2", "SMC4",
    "TPX2", "FOS", "HSP90AB1", "TXN"),
  lateral_gene_patterns = c("RP[LS].*"),
  noisy_gene_names = c(),
  noisy_gene_patterns = c(),
  assay_name = "full_run",
  ignore_warnings = TRUE,
  result_name = "metacells",
  envname = "r-decomp",
  return_model = FALSE,
  seed = 42,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{sce}{<SingleCellExperiment object> SCE object}

\item{target_number_of_metacells}{<integer or NULL> default 100;
Number of approximate number of MetaCells (latents) to be defined.
If NULL or 0, it will fallback to use approximately 100 cells for each
MetaCells.}

\item{num_parallel_piles}{<NULL or integer> Number of parallel processes.
If set to NULL, the method will estimate the appropriate number of
parallel processes that the machine can run.}

\item{quality_filters}{<bool> default TRUE; Whether to perform quality
filtering steps prior to the estimation of MetaCells. It is strongly
suggested to perform appropriate quality filtering prior to the estimation
of MetaCells.}

\item{min_umi}{<integer> default 800; Minimum number of Unique Molecular
Identifier (UMI) to consider a cell.}

\item{max_umi}{<integer> default 20000; Maximum number of Unique Molecular
Identifier (UMI) to consider a cell.}

\item{excluded_gene_names}{<list of character> default c('XIST, MALAT1');
List of blacklisted genes. These genes  should not be considered during the
analysis. example : XIST is related to inactivation of the X chromosome;
MALAT1 is associated with metastatic processes.}

\item{excluded_gene_patterns}{<list of characters> default c('MT-.*');
List of blacklisted gene patterns, which should not be considered during the
analysis. example : MT- (all the genes starting with MT-; basically all
the mitochondrial genes).}

\item{max_excluded_gene_fraction}{<double> default 0.25; maximum fraction
of genes excluded for a given cell. If a cell has more than 25% of its
expressed genes excluded, then is not considered for the following analysis.}

\item{lateral_gene_name}{<list of characters> default c("AURKA", "MCM3",
"MCM4", "MCM7", "MKI67", "PCNA", "RRM2", "SMC4", "TPX2", "FOS",
"HSP90AB1", "TXN"); (Cell-Cycle and Stress related genes)
List of genes which are not considered during the
community detection step, thus they do not contribute into defining
MetaCells, even though they are still used to build the KNN graph.}

\item{lateral_gene_patterns}{<list of characters> default c("RP[LS].*");
(Ribosomial Genes). List of gene patterns which are not considered during
the community detection step, thus they do not contribute into defining
MetaCells, even though they are still used to build the KNN graph.}

\item{noisy_gene_names}{<list of characters> default c();
List of genes which are expected to be highly variable between cells without
bearing any information in their variability. This genes will not be
considered Highly Variable Genes even if they are Highly Variable.}

\item{noisy_gene_patterns}{<list of characters> default c();
List of gene patterns which are expected to be highly variable between
cells without bearing any information in their variability. This genes
will not be considered Highly Variable Genes even if they are Highly
Variable.}

\item{assay_name}{<character> default 'full_run';
Name assigned to the metacell assay inside the AnnData Object.
It is useful when you want to retrieve also the full model, thus enabling you
to set a specific name for the results of the MetaCell run.}

\item{ignore_warnings}{<bool> default TRUE; Whether to ignore warnigns}

\item{result_name}{<character> default 'metacells';
Name used to store the result in the SingleCellExperiment object.}

\item{envname}{<character> default 'r-decomp';
Specify the name of the python virtual
environment to be used. If it does not exists it will create one and use it.}

\item{return_model}{<bool> default FALSE; Whether to return also
the model and not only
the SingleCellExperiment object.}

\item{seed}{<integer> default 42; to set the seed for reproducibility.}

\item{verbose}{<bool> default FALSE; Whether to be prompted with message
for each step of the analysis.}

\item{...}{<extra arguments for the divide_and_conquer_pipeline function>}
}
\value{
either a SingleCellExperiment object with W matrix
representation or the SingleCellExperiment object and the model
used to perform MetaCell analysis.
}
\description{
This function performs MetaCell analysis on a SingleCellExperiment object
and stores the results directly in it. This function leverages the python
implementation of MetaCell to achieve higher performance.
}
\examples{
library(packageX)
data(sce)
partition_metacells(sce, target_number_of_metacells = 2, min_umi = 5)
}
