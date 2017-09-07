#' External Annotation.
#'
#' A dataset containing external information that can be used to tune the weights in
#' the polygenic gene score using boosted regression tree models.
#'
#' @format A data frame with 927 rows and 48 variables:
#' \describe{
#'   \item{SNP_ID}{SNP rs ID}
#'   \item{BMI_beta_adj}{Standardized regression coefficient for association of a SNP and BMI obtained from external consortium GIANT}
#'   \item{CARDIOGRAM_beta_adj}{Standardized regression coefficient for association of a SNP and CAD obtained from external consortium CARDIOGRAM}
#'   \item{DIAGRAM_beta_adj}{Standardized regression coefficient for association of a SNP and T2D obtained from external consortium DIAGRAM}
#'   \item{HDL_beta_adj}{Standardized regression coefficient for association of a SNP and HDL obtained from external consortium Global Lipids Genetics Consortium Results}
#' }
#' @source
#' \url{http://www.diagram-consortium.org}
#' \url{http://csg.sph.umich.edu/abecasis/public/lipids2013/}
#' \url{http://www.cardiogramplusc4d.org/}
#' \url{http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium}
"annotation_data"
