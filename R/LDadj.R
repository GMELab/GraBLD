#' LD adjustment.
#'
#' The function calculates LD adjustment for each individual SNP based on the number of
#'   up/down stream SNPs for genome-wide SNPs.
#'
#' @param source_data the name of the file which the genotype data are to be read from.
#'    Each row of the matrix appears as one line of the file. Could be an absolute path to the
#'    file or the name of the file assuming in the present directory \code{\link{getwd}()}.
#'
#' @param geno_raw a matrix of the raw genotype data, contains \code{n} individuals (by row) and \code{m} SNPs (by column).
#'    Either \code{source_data} or \code{geno_raw} can be provided.
#'
#' @param size an integer for the number of SNPs that should be included in a block. Usually, for genome-wide
#'    datasets, there are about 2 million SNPs and +/- 300 SNPs is roughly equivalent to 1Mb
#'    physical distance, thus 300 is set as the default size.
#'
#' @param max_size an integer for the maximum size of SNPs to be standardized at a time, the default is 100000.
#'    This option can be ignored for data with less than 1 million SNPs.
#'
#' @param chr an integer for the chromosome of genotype data supplied,
#'    is not required and only used to name the output file. It is recommended to compute LD adjustments by
#'    chromosome to save memory.
#'
#' @param NAval the missing value used in the output data matrix. The default is NA, for use with
#'    PLINK, set \code{NAval} = -9.
#'
#' @param write a logic indicating whether the results should be written in a text file. If TRUE,
#'    the user should also provide the output file name by specifying \code{outname}; otherwise the default is
#'    `LDadj_\code{chr}_size_\code{size}_SNPs.txt'.
#'
#' @param outname a character giving the name of the output file if \code{write} is TRUE. Otherwise, the
#'    default name `LDadj_\code{chr}_size_\code{size}_SNPs.txt' is used.
#'
#' @return a numeric vector of LD adjustments with length matching the
#'    number of SNPs in the genotype data provided.
#'
#' @importFrom utils write.table
#'
#' @details
#' For large datasets, it is recommended to run from the command line with \cr
#' \preformatted{for((i = 1; i <= chr; i++))
#' do
#' Rscript PerformLDadj.R size data_name ${i} &
#' done}
#'
#' where the R script \code{PerformLDadj.R} might look something like this, while additional options can be
#' added to the argument list:  \cr
#'
#' \preformatted{#!/bin/sh
#' rm(list = ls())
#' library('GraBLD')
#' args = (commandArgs(TRUE))
#' size = eval(parse(text=args[1]))
#' source_data = args[2]
#' chr = eval(parse(text=args[3]))
#' geno_data = load_geno(source_data = source_data, PLINK = TRUE, chr = chr)
#' geno_norm = full_normal_geno(geno_data)
#' LD_OUT <- LDadj(geno_raw = geno_norm, chr = chr, size = size, write = TRUE)}
#'
#'
#' @examples
#' data(geno)
#' LDadj(geno_raw = geno, chr = 1, size = 200)
#'
#' @references Guillaume Pare, Shihong Mao, Wei Q Deng (2017)
#' A machine-learning heuristic to improve gene score prediction of
#' polygenic traits Short title: Machine-learning boosted gene scores,
#' \emph{bioRxiv 107409}; doi: \url{https://doi.org/10.1101/107409};
#' \url{http://www.biorxiv.org/content/early/2017/02/09/107409}
#'
#' Pare, Guillaume, Shihong Mao, and Wei Q. Deng.
#' A method to estimate the contribution of regional genetic
#' associations to complex traits from summary association statistics.
#' \emph{Scientific reports} \strong{6} (2016): 27644.
#'
#'
#'

LDadj <- function(source_data = NULL, geno_raw, size = 300,
    chr = NULL, max_size = 1e+05, write = FALSE, outname = NULL,
    NAval = NA) {

    if (is.null(geno_raw) & is.null(source_data)) {
        stop("Please either provide the source of the genotype data or genotype data")
    }

    if (is.null(geno_raw)) {
        try(geno_data <- load_geno(source_data))
        if (is.null(geno_data)) {
            stop("Please ensure the data source is correct. No file found or readable.")
        } else {
            print("Genotype data loaded.")
            try(geno_norm <- full_normal_geno(geno_data, max_size = max_size,
                NAval = NAval))
        }
    } else {
        try(geno_norm <- full_normal_geno(geno_raw, max_size = max_size,
            NAval = NAval))
        if (is.null(geno_norm)) {
            stop("Please ensure the data format is correct.")
        } else {
            print("Genotype data loaded and normalized.")
        }
    }

    LDadj <- genomewideLD(geno_norm, size = size)
    LDadj[which(LDadj < 1)] <- 1
    chr_out <- rep(chr, length(LDadj))
    LDadjs <- cbind(chr_out, LDadj)

    if (write == TRUE) {
        if (is.null(outname)) {
            utils::write.table(LDadjs, paste("LDadj_", chr,
                "_size_", size, "_SNPs.txt", sep = ""), col.names = F,
                row.names = T, quote = F, sep = "\t")
        } else {
            utils::write.table(LDadjs, paste(outname), col.names = F,
                row.names = T, quote = F, sep = "\t")
        }
    } else {
        return(LDadjs)
    }
}
