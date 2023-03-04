#' Population doubling example dataset
#'
#' A dataset containing the raw counts for bulk RNA-seq of fibroblasts 
#' at 16, 24, 46, 64, and 72 population doublings, with 3 samples per condition
#'
#' @format A data frame with 61541 rows and 16 columns:
#' \describe{
#'   \item{Geneid}{Ensembl IDs of genes}
#'   \item{SRR1660543}{Raw counts for 1st sample}
#'   ...
#' }
#'
"pd_raw"

#' Growing, quiescent, senescent example dataset
#'
#' A dataset containing the raw counts for bulk RNA-seq of fibroblasts 
#' under growing, quiescent, senescent, and deep senescent conditions, with 2 samples
#' per condition
#'
#' @format A data frame with 62703 rows and 9 columns:
#' \describe{
#'   \item{Geneid}{Ensembl IDs of genes}
#'   \item{SRR9601022}{Raw counts for 1st sample}
#'   ...
#' }
#'
"gqs_raw"

#' Neural stem cell (NSC) example dataset
#'
#' A dataset containing the raw counts for scRNA-seq of 
#' dormant, primed, early-stage activated, and late-stage activated neural stem cells (NSC),
#' with 2785 total cells
#' 
#'
#' @format A data frame with 27998 rows and 2785 columns:
#' \describe{
#'   \item{AAACCTGCACCAGTTA-1}{Raw counts for 1st sample}
#'   \item{AAACCTGCATGGTCAT-1}{Raw counts for 2nd sample}
#'   ...
#' }
#'
"nsc_raw"

#' Rat embryonic fibroblast (REF) example dataset
#'
#' A dataset containing the raw counts for bulk RNA-seq of 
#' rat embryonic fibroblasts (REF) under 2, 3, 4, 6, 8, 10, 12, 14, and 16 
#' days of serum starvation, with 3 samples per group
#' 
#' 
#'
#' @format A data frame with 30560 rows and 28 columns:
#' \describe{
#'   \item{Geneid}{Ensembl IDs of genes}
#'   \item{SRR8353399}{Raw counts for 1st sample}
#'   ...
#' }
#'
"ref_raw"