#'
#' Gradient Boosted and LD adjusted Prediction.
#'
#' The function returns the prediction r-squared on the target population using polygenic gene score
#'   based on the GraBLD heuristic.
#'
#' @param source_data the name of the file which the genotype data are to be read from.
#'    Each row of the matrix appears as one line of the file. Could be an absolute path to the
#'    file or the name of the file assuming in the present directory \code{\link{getwd}()}.
#'
#' @param PLINK a logic indicating whether the supplied file is of the .raw format from PLINK, if not, the first six
#'    columns will not be removed and all columns will be read in.
#'
#' @param chr an integer indicating the maximum number of chromosomes to be read in, this option is used in combination with
#'    \code{source_data}, to perform analysis by each chromosome. In this case, the file name should follow:
#'    ``source_data_i.raw'' for all \eqn{i \le chr}. For example, \code{Geno_Matrix_23.raw}.
#'
#' @param SNPnames a vector of characters for the names of SNPs used, these are only used in the output of GraBLD weights.
#'
#' @param geno_raw the genotype matrix assuming each row is the individual and each column is the SNP. Can be skipped
#'    if \code{source_data} was provided.
#'
#' @param max_size an integer for the maximum size of SNPs to be standardized at a time, the default is 100000.
#'    This option can be ignored for data with less than 1 million SNPs.
#'
#' @param NAval the missing value used in the output data matrix. The default is NA, for use with
#'    PLINK, set \code{NAval} = -9.
#'
#' @param Pheno a numeric vector of quantitative traits with the same length as the number of rows in the genotype matrix.
#'
#' @param trait_name a character indicating the name of the quantitative trait.
#'
#' @param WRITE a logic indicating whether the results of the GraBLD weights should be written to a file with
#'    file name \code{trait_name_gbm_beta.txt}.
#'
#' @param LDadjVal a numeric vector of LD adjusted scores with length matching the number of SNPs in the
#'    genotype matrix. This can be taken from the output of \code{\link{LDadj}()}.
#'
#' @param gbmVal a numeric vector of gradient boosted weights with length matching the number of SNPs in
#'    the genotype matrix. This can be taken from the output of \code{\link{GraB}()}.
#'
#' @return if \code{Pheno} is supplied, both the polygenic gene score as well as the prediction R-squared (adjusted)
#'    are returned, otherwise only the polygenic gene score is returned.
#'
#' @importFrom data.table fread
#' @importFrom gbm gbm
#' @importFrom stats lm, sd, as.formula
#' @importFrom utils write.table
#'
#' @examples
#' data(geno)
#' data(univariate_beta)
#' data(annotation_data)
#' LD_val <- LDadj(geno_raw = geno, chr = 1, size = 200)
#' gbm_val <- list()
#' for (j in 1:5){
#' gbm_val[[j]] <- GraB(betas = univariate_beta, annotations = annotation_data,
#' trait_name = 'BMI', steps = j, validation = 5)
#' }
#' data(BMI)
#' gs <- GraBLD.score(geno_raw = geno, LDadjVal=LD_val, gbmVal = unlist(gbm_val),
#'  trait_name='BMI', Pheno = BMI[,1])
#'
#' @references Guillaume Pare, Shihong Mao, Wei Q Deng (2017)
#' A machine-learning heuristic to improve gene score prediction of
#' polygenic traits Short title: Machine-learning boosted gene scores,
#' \emph{bioRxiv 107409}; doi: \url{https://doi.org/10.1101/107409};
#' \url{http://www.biorxiv.org/content/early/2017/02/09/107409}
#'
#'

GraBLD.score <- function(source_data = NULL, chr = NULL, geno_raw,
    PLINK = TRUE, SNPnames = NULL, max_size = 1e+05, NAval = NA,
    Pheno = NULL, LDadjVal, gbmVal, trait_name = NULL, WRITE = FALSE) {

    if (is.null(LDadjVal) | is.null(gbmVal)) {
        stop("Please provide the LD adjustments and the gbm tuned weights.")
    }

    if (length(gbmVal) != dim(LDadjVal)[1]) {
        stop("Please ensure the length of the LD adjustments matches that of the gbm tuned weights.")
    } else {
        gbm_adj = gbmVal/LDadjVal[, 2]
    }

    geno_norm = as.numeric()

    if (is.null(chr) & is.null(source_data)) {

        geno_norm = full_normal_geno(geno_raw, NAval = NAval,
            max_size = max_size)

    } else {

        for (i in 1:chr) {
            geno_pre_chr = load_geno(source_data, PLINK = PLINK,
                chr = chr)
            geno_norm_chr = full_normal_geno(geno_pre_chr,
                NAval = NAval, max_size = max_size)
            geno_norm = cbind(geno_norm, geno_norm_chr)
        }
    }

    if (dim(geno_norm)[2] != length(gbm_adj)) {
        stop("Please ensure the number of SNPs matches the length of the GraBLD weights.")
    }

    gene_score = as.matrix(geno_norm) %*% gbm_adj

    if (WRITE == TRUE) {

        gbm_adj_title = c("gbm_beta", gbm_adj)

        if (is.null(SNPnames)) {
            if (length(SNPnames) != length(gbm_adj_title)) {
                warnings("the number of SNPs supplied does not match the length of the GraBLD weights,
               only the GraBLD weights are outputted.")
                utils::write.table(gbm_adj_title, paste(trait_name,
                  "_gbm_beta.txt", sep = ""), col.names = F,
                  row.names = F, quote = F, sep = "\t")
            } else {
                utils::write.table(cbind(SNPnames, gbm_adj_title),
                  paste(trait_name, "_gbm_beta.txt", sep = ""),
                  col.names = F, row.names = F, quote = F,
                  sep = "\t")
            }
        } else {
            utils::write.table(gbm_adj_title, paste(trait_name,
                "_gbm_beta.txt", sep = ""), col.names = F,
                row.names = F, quote = F, sep = "\t")
        }
    }

    if (is.null(Pheno)) {

        return(GraBLD.gs = gene_score)

    } else {
        Pheno_norm = (Pheno - mean(Pheno, na.rm = T))/stats::sd(Pheno,
            na.rm = T)
        adjR2 = summary(stats::lm(Pheno_norm ~ gene_score))$adj.r.squared
        return(list(GraBLD.gs = gene_score, adj.R2 = adjR2))
    }
}

