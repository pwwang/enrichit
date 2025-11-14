#' Fetch GMT file
#'
#' @param gmtpath Path to the GMT file.
#' Can be one of the following:
#' * A local file path (will be checked for existence)
#' * A URL starting with http:// or https:// (will be downloaded and cached)
#' * A built-in library name (e.g., "KEGG", "GO_Biological_Process", "Reactome")
#' * An Enrichr library name (e.g., "ChEA_2016", "GWAS_Catalog_2019")
#'
#' Built-in libraries include:
#' * "BioCarta" or "BioCarta_2016"
#' * "GO_Biological_Process" or "GO_Biological_Process_2025"
#' * "GO_Cellular_Component" or "GO_Cellular_Component_2025"
#' * "GO_Molecular_Function" or "GO_Molecular_Function_2025"
#' * "KEGG", "KEGG_Human", "KEGG_2021", or "KEGG_2021_Human"
#' * "Hallmark", "MSigDB_Hallmark", or "MSigDB_Hallmark_2020"
#' * "Reactome", "Reactome_Pathways", or "Reactome_Pathways_2024"
#' * "WikiPathways", "WikiPathways_2024", "WikiPathways_Human", or "WikiPathways_2024_Human"
#'
#' For Enrichr libraries, see https://maayanlab.cloud/Enrichr/#libraries
#'
#' @return The path to the GMT file (local, built-in, or cached).
#' If `gmtpath` is a URL or Enrichr library, the file is downloaded to a cache
#' directory if not already present and the path to the cached file is returned.
#' @importFrom digest digest
#' @importFrom utils download.file
#' @importFrom httr parse_url
#' @importFrom rlang %||%
#' @importFrom tools file_path_sans_ext
#' @keywords internal
#' @examples
#' \donttest{
#' if (FALSE) {
#' # Use a built-in library
#' FetchGMT("KEGG")
#' FetchGMT("GO_Biological_Process")
#'
#' # Fetch from Enrichr by library name
#' FetchGMT("ChEA_2016")
#' FetchGMT("GWAS_Catalog_2019")
#'
#' # Fetch a GMT file from a URL (will be cached)
#' gmtpath_url <- paste0(
#'   "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=",
#'   "Data_Acquisition_Method_Most_Popular_Genes"
#' )
#' gmt_file_url <- FetchGMT(gmtpath_url)
#' # Example output: ~/.cache/enrichit/xxxxxxxx/Data_Acquisition_Method_Most_Popular_Genes.gmt
#'
#' # Fetching the same URL again uses the cache
#' gmt_file_url_cached <- FetchGMT(gmtpath_url)
#' print(identical(gmt_file_url, gmt_file_url_cached)) # TRUE
#'
#' # Use a local GMT file
#' local_gmt_path <- tempfile(fileext = ".gmt")
#' writeLines(c(
#'   "GENESET1\tdesc1\tGENE1\tGENE2",
#'   "GENESET2\tdesc2\tGENE3\tGENE4"
#' ), local_gmt_path)
#' FetchGMT(local_gmt_path)
#' # Returns: /tmp/RtmpXXXXXX/filexxxx.gmt (original path)
#' }
#' }
FetchGMT <- function(gmtpath) {
    if (gmtpath %in% c("BioCarta", "BioCarta_2016")) {
        gmtpath <- system.file("extdata", "BioCarta_2016.gmt.gz", package = "enrichit")
    } else if (gmtpath %in% c("GO_Biological_Process", "GO_Biological_Process_2025")) {
        gmtpath <- system.file("extdata", "GO_Biological_Process_2025.gmt.gz", package = "enrichit")
    } else if (gmtpath %in% c("GO_Cellular_Component", "GO_Cellular_Component_2025")) {
        gmtpath <- system.file("extdata", "GO_Cellular_Component_2025.gmt.gz", package = "enrichit")
    } else if (gmtpath %in% c("GO_Molecular_Function", "GO_Molecular_Function_2025")) {
        gmtpath <- system.file("extdata", "GO_Molecular_Function_2025.gmt.gz", package = "enrichit")
    } else if (gmtpath %in% c("KEGG", "KEGG_Human", "KEGG_2021", "KEGG_2021_Human")) {
        gmtpath <- system.file("extdata", "KEGG_2021_Human.gmt.gz", package = "enrichit")
    } else if (gmtpath %in% c("Hallmark", "MSigDB_Hallmark", "MSigDB_Hallmark_2020")) {
        gmtpath <- system.file("extdata", "MSigDB_Hallmark_2020.gmt.gz", package = "enrichit")
    } else if (gmtpath %in% c("Reactome", "Reactome_Pathways", "Reactome_Pathways_2024")) {
        gmtpath <- system.file("extdata", "Reactome_Pathways_2024.gmt.gz", package = "enrichit")
    } else if (gmtpath %in% c("WikiPathways", "WikiPathways_2024", "WikiPathways_Human", "WikiPathways_2024_Human")) {
        gmtpath <- system.file("extdata", "WikiPathways_2024_Human.gmt.gz", package = "enrichit")
    } else if (grepl("^[a-zA-Z0-9_-]+$", gmtpath)) {
        # Check if the GMT file is available in Enrichr
        gmtpath <- paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", gmtpath)
    }
    # Check if the path is a URL
    if (grepl("^https?://", gmtpath)) {
        # Use a persistent user cache directory
        cache_dir <- tools::R_user_dir("enrichit", which = "cache")
        # Generate a unique filename based on the URL hash
        file_hash <- substr(digest(gmtpath, algo = "md5"), 1, 8)
        cache_dir <- file.path(cache_dir, file_hash)
        if (!dir.exists(cache_dir)) {
            dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        }
        # Parse the URL to get the libraryName
        parsed_url <- parse_url(gmtpath)
        library_name <- parsed_url$query$libraryName %||% file_path_sans_ext(sub("\\.gz$", "", basename(gmtpath)))

        # Create the target file path
        target_file <- file.path(cache_dir, paste0(library_name, ".gmt"))

        # Check if the file exists in cache and is not empty
        if (file.exists(target_file) && file.info(target_file)$size > 0) {
            message("Using cached file: ", target_file)
            return(target_file)
        } else {
            message("Downloading GMT file from URL: ", gmtpath)
            # Download the file
            tryCatch({
                download_status <- download.file(gmtpath, target_file, mode = "wb", quiet = TRUE)
                if (download_status != 0) {
                   stop("[EnrichIt] Failed to download file. Status code: ", download_status)
                }
                # Verify download integrity (basic check: file exists and is not empty)
                if (!file.exists(target_file) || file.info(target_file)$size == 0) {
                    stop("Downloaded file is empty or does not exist: ", target_file)
                }
                message("File downloaded and cached successfully: ", target_file)
                return(target_file)
            }, error = function(e) {
                # Clean up potentially incomplete file on error
                if (file.exists(target_file)) {
                    unlink(target_file)
                }
                stop("[EnrichIt] Error downloading ", gmtpath, ": ", e$message)
            })
        }
    } else {
        # If it's a local file, check existence and return the path
        if (file.exists(gmtpath)) {
            return(gmtpath)
        } else {
            stop("[EnrichIt] The specified local GMT file does not exist: ", gmtpath)
        }
    }
}


#' Read GMT file
#'
#' @param gmtfile Path to the GMT file
#' @param swap_name_desc_if_needed Logical, whether to swap name and description fields.
#' They will be swapped only if:
#' * `swap_name_desc_if_needed` is `TRUE`; and
#' * The descriptions are not empty; and
#' * The descriptions are shorter than the names; and
#' * The descriptions are not ID-like (i.e., hsa00001, or 123456).
#' @return A list of gene sets, the names of the list are the gene set names,
#' and each element is a vector of gene names.
#' @export
#' @examples
#' \donttest{
#' if (FALSE) {
#' # Example GMT file content
#' gmtfile <- system.file("extdata", "KEGG_2021_Human.gmt", package = "enrichit")
#' # Read the GMT file
#' ParseGMT(gmtfile)
#' }
#' }
ParseGMT <- function(gmtfile, swap_name_desc_if_needed = TRUE) {
    gmtfile <- FetchGMT(gmtfile)
    # Support gzipped files
    con <- if (grepl("\\.gz$", gmtfile)) gzfile(gmtfile) else gmtfile
    on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    # Parse lines into gene sets
    libraries <- lapply(lines, function(line) {
        fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
        list(name = fields[1], desc = fields[2], genes = fields[-(1:2)])
    })
    # Remove empty gene sets and empty genes
    libraries <- Filter(function(lib) length(lib$genes) > 0, libraries)
    libraries <- lapply(libraries, function(lib) {
        lib$genes <- lib$genes[nzchar(lib$genes)]
        lib
    })
    # Swap name and desc if needed
    name_nchars <- vapply(libraries, function(lib) nchar(lib$name), integer(1))
    desc_nchars <- vapply(libraries, function(lib) nchar(lib$desc), integer(1))
    desc_prefix <- sub("[0-9]+$", "", libraries[[1]]$desc)
    desc_are_ids <- vapply(libraries, function(lib) startsWith(lib$desc, desc_prefix) | grepl("^[0-9]+$", lib$desc), logical(1))
    if (
        swap_name_desc_if_needed &&
        all(desc_nchars > 0) &&
        all(desc_nchars < name_nchars) &&
        !all(desc_are_ids)
    ) {
        warning("Swapping name and desc fields in GMT file, as desc is shorter. Set 'swap_name_desc_if_needed' to FALSE to disable this.")
        libraries <- lapply(libraries, function(lib) {
            tmp <- lib$name
            lib$name <- lib$desc
            lib$desc <- tmp
            lib
        })
    }
    # Remove duplicates by name
    unique_names <- !duplicated(vapply(libraries, `[[`, character(1), "name"))
    libraries <- libraries[unique_names]
    # Set names and remove name/desc field
    names(libraries) <- vapply(libraries, `[[`, character(1), "name")
    lapply(libraries, `[[`, "genes")
}
