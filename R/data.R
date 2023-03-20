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

#' 16-74 population doublings dataset annotations
#'
#' The annotations for the 16-74 population doublings bulk RNA-seq dataset
#' 
#' @format A data frame with 15 rows and 2 columns:
#' \describe{
#'   \item{sample_id}{IDs of the cells}
#'   \item{population_doublings}{The number of population doublings the sample has gone through}
#' }
#'
"pd_ann"

#' Growing, quiescent, senescent dataset annotations
#'
#' The annotations for the growing, quiescent, senescent bulk RNA-seq dataset
#' 
#' @format A data frame with 8 rows and 2 columns:
#' \describe{
#'   \item{sample_id}{IDs of the cells}
#'   \item{Condition}{The state of the samples (growing, quiescence, senescence, deep senescence)}
#' }
#'
"gqs_ann"

#' Neural stem cell (NSC) dataset annotations
#'
#' The annotations for the NSC scRNA-seq dataset
#' 
#' @format A data frame with 2785 rows and 3 columns:
#' \describe{
#'   \item{sample_id}{IDs of the cells}
#'   \item{celltype}{The dormancy state (dormant, primed, activated (early), or activted (late) )}
#'   \item{age}{The age (old or young)}
#' }
#'
"nsc_ann"

#' Rat embryonic fibroblast (REF) dataset annotations
#'
#' The annotations for the bulk RNA-seq of 
#' rat embryonic fibroblasts (REF) under 2, 3, 4, 6, 8, 10, 12, 14, and 16 
#' days of serum starvation
#' 
#' @format A data frame with 27 rows and 3 columns:
#' \describe{
#'   \item{sample_id}{IDs of the cells}
#'   \item{Condition}{The dormancy state (all Q, for quiescent))}
#'   \item{Dataset}{The dataset (all REF)}
#'   \item{Day}{The number of serum starvation days (2-16)}
#' }
#'
"ref_ann"

#' Gene map for neural stem cell (NSC) dataset
#'
#' The Ensembl IDs corresponding to the gene names in the NSC dataset
#' 
#' 
#' @format A data frame with 27998 rows and 2 columns:
#' \describe{
#'   \item{initial_alias}{gene names}
#'   \item{converted_alias}{Ensembl IDs}
#' }
#'
"nsc_genemap"

#' Homology table for mouse and rat Ensembl IDs
#'
#' The corresponding Ensembl IDs for mouse and rat, as well as the homology type and confidence
#' 
#' @format A data frame with 128614 rows and 2 columns:
#' \describe{
#'   \item{Gene stable ID}{Mouse Ensembl gene IDs}
#'   \item{Gene stable ID version}{Mouse Ensembl gene IDs, with version}
#'   \item{Transcript stable ID}{Mouse Ensembl transcript IDs}
#'   \item{Transcript stable ID version}{Mouse Ensembl transcript IDs, with version}
#'   \item{Rat gene stable ID}{Rat Ensembl gene IDs}
#'   \item{Rat orthology confidence}{Confidence in homology relationship (0=low, 1=high)}
#'   \item{Rat homology type}{Homology type (one-to-one, one-to-many, many-to-many)}
#' }
#'
"rat_mouse_hm"

#' Homology table for human and rat Ensembl IDs
#'
#' The corresponding Ensembl IDs for human and rat, as well as the homology type and confidence
#' 
#' @format A data frame with 179981 rows and 2 columns:
#' \describe{
#'   \item{Gene stable ID}{Human Ensembl gene IDs}
#'   \item{Gene stable ID version}{Human Ensembl gene IDs, with version}
#'   \item{Transcript stable ID}{Human Ensembl transcript IDs}
#'   \item{Transcript stable ID version}{Human Ensembl transcript IDs, with version}
#'   \item{Rat gene stable ID}{Rat Ensembl gene IDs}
#'   \item{Rat orthology confidence}{Confidence in homology relationship (0=low, 1=high)}
#'   \item{Rat homology type}{Homology type (one-to-one, one-to-many, many-to-many)}
#' }
#'
"rat_human_hm"

