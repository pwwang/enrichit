#' Build the numbers for the tests
#'
#' @param userlist Vector of user-provided genes
#' @param geneset Vector of genes in the gene set
#' @param background Vector of all genes in the universe or a number of genes
#' in the universe
#' @return A list with the following elements:
#' * k: Number of genes in the user list that are also in the gene set
#' * M: Number of genes in the gene set
#' * n: Number of genes in the user list
#' * N: Total number of genes in the universe
#' @keywords internal
BuildNumbers <- function(userlist, geneset, background) {
  k <- length(intersect(userlist, geneset))
  M <- length(geneset)
  n <- length(userlist)
  N <- if (is.numeric(background)) background else length(background)
  return(list(k = k, M = M, n = n, N = N))
}


#' Compute statistic p-value by given numbers
#'
#' @param k Number of genes in both user list and gene set
#' @param M Number of genes in the gene set
#' @param n Number of genes in the user list
#' @param N Total number of genes in the universe
#' @param method Method for computing p-value, either "fisher" or "hypergeometric"
#' @return p-value from the specified test
#' @keywords internal
#' @importFrom stats fisher.test phyper
ComputePValue <- function(k, M, n, N, method = c("fisher", "hypergeometric")) {
    method <- match.arg(method)
    if (method == "fisher") {
        mat <- matrix(c(k, n - k, M - k, N - M - n + k), nrow = 2)
        return(fisher.test(mat, alternative = "greater")$p.value)
    } else {  # method == "hypergeometric")
        return(phyper(k - 1, M, N - M, n, lower.tail = FALSE))
    }
}


#' Check user list against gene sets
#'
#' @param userlist Vector of user-provided genes
#' @param genesets List of gene sets
#' @param use_matched_only Logical, whether to return only matched genes
#' Otherwise the genes in the user list will be all used.
#' This will affect the number of genes in the user list when computing the
#' p-value.
#' @param show_all_unmatched Logical, whether to show all unmatched genes. Or instead
#' show only the first 10 unmatched genes.
#' @return A vector of user-provided genes that are also in the gene sets
#' @keywords internal
CheckUserList <- function(userlist, genesets, use_matched_only, show_all_unmatched = FALSE) {
    gs <- unique(unlist(genesets))
    if (length(userlist) == 0) {
        stop("[EnrichIt] User list is empty.")
    }
    userlist <- unique(userlist)
    unmatched <- setdiff(userlist, gs)
    if (length(unmatched) > 0) {
        if (isFALSE(show_all_unmatched)) {
            n_unmatched <- length(unmatched)
            if (n_unmatched > 10) {
                unmatched <- unmatched[1:10]
                warning("User list contains ", n_unmatched, " unmatched genes. Showing first 10: ", paste(unmatched, collapse = ", "))
            } else {
                warning("User list contains ", n_unmatched, " unmatched genes: ", paste(unmatched, collapse = ", "))
            }
        } else {
            warning("User list contains unmatched genes: ", paste(unmatched, collapse = ", "))
        }
    }
    if (isTRUE(use_matched_only)) {
        userlist <- intersect(userlist, gs)
    }
    if (length(userlist) == 0) {
        stop("[EnrichIt] User list contains no genes that are in the gene sets.")
    }
    return(userlist)
}


#' Perform enrichment analysis
#'
#' @param userlist Vector of user-provided genes
#' @param dbs List of gene sets or paths to GMT files
#' It can be a vector of database names. You can set the names of the vector, which
#' will be used as the database names.
#' Otherwise a list is expected, where each element is a database (e.g. parsed from a gmt file).
#' If a database is given directly (without a name), the expression of this argument will
#' be used as the database name.
#' @param method Method for computing p-value, either "fisher" or "hypergeometric"
#' When `style` is "enrichr", the method defaults to "fisher".
#' When `style` is "clusterProfiler", the method defaults to "hypergeometric".
#' @param use_matched_only Logical, whether to use only matched genes against the gene sets.
#' This will affect the number of genes in the user list when computing the
#' p-value.
#' By default, when `style` is "enrichr", this is set to FALSE.
#' When `style` is "clusterProfiler", this is set to TRUE.
#' @param padjust_method Method for adjusting p-values, either "BH", "bonferroni",
#' "holm", "hochberg", "hommel", "BY", "fdr"
#' @param background Vector of all genes in the universe or a number of genes
#' in the universe. If NULL, the number of genes in the gene set will be used.
#' For "enrichr", the default is 20,000.
#' For "clusterProfiler", the default is the number of unique genes in the gene set.
#' Note that for "enrichr", if a vector is given, the length of it will be used, no
#' checking will be done to see if userlist and genes from dbs are in the vector.
#' @param style Style of the output, either "enrichr" or "clusterProfiler"
#' @param return_all Logical, whether to return all results (all gene sets in dbs)
#' or only those with at least one gene in the user list.
#' @returns A data frame with the results.
#' When `style` is "enrichr", the columns are:
#' * Database: Name of the database
#' * Term: Name of the term
#' * Overlap: Number of genes in the user list and the gene set
#' * P.value: p-value from the enrichment test
#' * Adjusted.P.value: Adjusted p-value from the enrichment test
#' * Odds.Ratio: Odds ratio from the enrichment test
#' * Combined.Score: Combined score from the enrichment test
#' * Genes: Genes in the user list that are also in the gene set
#' * Rank: Rank of the term based on the combined score
#'
#' When `style` is "clusterProfiler", the columns are:
#' * ID: Name of the term
#' * Description: Description of the term
#' * GeneRatio: Ratio of genes in the user list and the gene set
#' * BgRatio: Ratio of genes in the gene set and the universe
#' * Count: Number of genes in the user list that are also in the gene set
#' * pvalue: p-value from the enrichment test
#' * p.adjust: Adjusted p-value from the enrichment test
#' * qvalue: Q-value from the enrichment test
#' * geneID: Genes in the user list that are also in the gene set
#' * Database: Name of the database
#' @export
#' @importFrom rlang %||%
#' @importFrom stats p.adjust
#' @importFrom tools file_path_sans_ext
#' @importFrom qvalue qvalue
EnrichIt <- function(
    userlist, dbs,
    method = ifelse(tolower(style[1]) == "enrichr", "fisher", "hypergeometric"),
    use_matched_only = ifelse(tolower(style[1]) == "enrichr", FALSE, TRUE),
    padjust_method = c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY", "fdr"),
    background = NULL,
    style = c("enrichr", "Enrichr", "clusterProfiler", "ClusterProfiler", "clusterprofiler"),
    return_all = FALSE
) {
    method <- match.arg(method, c("fisher", "hypergeometric"))
    padjust_method <- match.arg(padjust_method)
    style <- match.arg(style)
    style <- tolower(style)
    dbs_name <- deparse(substitute(dbs))

    if (is.character((dbs))) {
        if (is.null(names(dbs))) {
            names(dbs) <- file_path_sans_ext(gsub("\\.gz$", "", basename(dbs)))
        }
        dbs <- lapply(dbs, FetchGMT)
        dbs <- lapply(dbs, ParseGMT)
    } else if (!is.list(dbs)) {
        stop("[EnrichIt] Invalid input for dbs.")
    } else {
        if (all(sapply(dbs, is.character))) {
            dbs <- list(dbs)
            names(dbs) <- dbs_name
        }
    }

    results <- lapply(names(dbs), function(db_name) {
        db <- dbs[[db_name]]
        db_gene_count <- length(unique(unlist(db)))
        background <- background %||% ifelse(style == "enrichr", max(20000, db_gene_count), db_gene_count)
        background <- 10922
        ul <- CheckUserList(userlist, db, use_matched_only = use_matched_only)
        res <- lapply(seq_along(db), function(i) {
            term_name <- names(db)[i]
            genes <- db[[term_name]]
            numbers <- BuildNumbers(ul, genes, background)
            pvalue <- ComputePValue(numbers$k, numbers$M, numbers$n, numbers$N, method)
            if (style == "enrichr") {
                # oddsRatio = (1.0 * a * d) / Math.max(1.0 * b * c, 1)
                #  where: a are the overlapping genes,
                #  b are the genes in the annotated set - overlapping genes,
                #  c are the genes in the input set - overlapping genes, and
                #  d are the 20,000 genes (or total genes in the background) - genes in the annotated set - genes in the input set + overlapping genes
                a <- numbers$k
                b <- numbers$M - numbers$k
                c <- numbers$n - numbers$k
                d <- numbers$N - numbers$M - numbers$n + numbers$k
                odds_ratio <- (a * d) / max(1, b * c)
                return(data.frame(
                    Database = db_name,
                    Term = term_name,
                    Overlap = paste(numbers$k, numbers$M, sep = "/"),
                    P.value = pvalue,
                    Odds.Ratio = odds_ratio,
                    Combined.Score = -log(pvalue) * odds_ratio,
                    Genes = paste(intersect(genes, ul), collapse = ";"),
                    stringsAsFactors = FALSE
                ))
            } else {
                # For other output styles, you can modify this part
                return(data.frame(
                    Database = db_name,
                    ID = paste0(db_name, "_", i),
                    Description = term_name,
                    GeneRatio = paste(numbers$k, numbers$n, sep = "/"),
                    BgRatio = paste(numbers$M, numbers$N, sep = "/"),
                    Count = numbers$k,
                    pvalue = pvalue,
                    geneID = paste(intersect(genes, ul), collapse = "/"),
                    stringsAsFactors = FALSE
                ))
            }
        })
        res <- do.call(rbind, res)
        if (style == "enrichr") {
            if (isFALSE(return_all)) {
                res <- res[res$Odds.Ratio > 0, , drop = FALSE]
            }
            res$k <- NULL
            res$Adjusted.P.value <- p.adjust(res$P.value, method = padjust_method)
            res$Rank <- rank(res$P.value, ties.method = "min")
            res <- res[
                order(res$Rank),
                c("Term", "Overlap", "P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes", "Rank", "Database"),
                drop = FALSE
            ]
        } else {  # style == "clusterprofiler"
            if (isFALSE(return_all)) {
                res <- res[res$Count > 0, , drop = FALSE]
            }
            res$p.adjust <- p.adjust(res$pvalue, method = padjust_method)
            res$qvalue <- tryCatch(qvalue(p=res$pvalue, lambda=0.05, pi0.method="bootstrap")$qvalues, error = function(e) NA)
            res <- res[
                order(res$pvalue),
                c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Database"),
                drop = FALSE
            ]
        }

        res
    })

    return(do.call(rbind, results))
}
