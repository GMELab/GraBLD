#' BMI Univariate beta regression coefficient.
#'
#' A dataset containing the univariate regression coefficient of 927 SNPs for association with BMI from the simulated data.
#'
#' @format A data frame with 927 rows and 2 variables:
#' \describe{
#' The univariate regression coefficients were computed from the simulated
#' datasets \code{BMI} and \code{geno} using: \cr
#'
#'    \code{coef(summary(lm(pheno_data ~ geno[,j])))[2,1]}.
#'
#' Both genotype data and phenotype data over individuals need to be standardized with mean = 0 and variance = 1.
#'   \item{V1}{chromosome}
#'   \item{V2}{BMI univariate regression coefficient}
#'   ...
#' }
"univariate_beta"
