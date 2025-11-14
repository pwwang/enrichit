# Fetch GMT file

Fetch GMT file

## Usage

``` r
FetchGMT(gmtpath)
```

## Arguments

- gmtpath:

  Path to the GMT file. Can be one of the following:

  - A local file path (will be checked for existence)

  - A URL starting with http:// or https:// (will be downloaded and
    cached)

  - A built-in library name (e.g., "KEGG", "GO_Biological_Process",
    "Reactome")

  - An Enrichr library name (e.g., "ChEA_2016", "GWAS_Catalog_2019")

  Built-in libraries include:

  - "BioCarta" or "BioCarta_2016"

  - "GO_Biological_Process" or "GO_Biological_Process_2025"

  - "GO_Cellular_Component" or "GO_Cellular_Component_2025"

  - "GO_Molecular_Function" or "GO_Molecular_Function_2025"

  - "KEGG", "KEGG_Human", "KEGG_2021", or "KEGG_2021_Human"

  - "Hallmark", "MSigDB_Hallmark", or "MSigDB_Hallmark_2020"

  - "Reactome", "Reactome_Pathways", or "Reactome_Pathways_2024"

  - "WikiPathways", "WikiPathways_2024", "WikiPathways_Human", or
    "WikiPathways_2024_Human"

  For Enrichr libraries, see https://maayanlab.cloud/Enrichr/#libraries

## Value

The path to the GMT file (local, built-in, or cached). If `gmtpath` is a
URL or Enrichr library, the file is downloaded to a cache directory if
not already present and the path to the cached file is returned.

## Examples

``` r
# \donttest{
if (FALSE) {
# Use a built-in library
FetchGMT("KEGG")
FetchGMT("GO_Biological_Process")

# Fetch from Enrichr by library name
FetchGMT("ChEA_2016")
FetchGMT("GWAS_Catalog_2019")

# Fetch a GMT file from a URL (will be cached)
gmtpath_url <- paste0(
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=",
  "Data_Acquisition_Method_Most_Popular_Genes"
)
gmt_file_url <- FetchGMT(gmtpath_url)
# Example output: ~/.cache/enrichit/xxxxxxxx/Data_Acquisition_Method_Most_Popular_Genes.gmt

# Fetching the same URL again uses the cache
gmt_file_url_cached <- FetchGMT(gmtpath_url)
print(identical(gmt_file_url, gmt_file_url_cached)) # TRUE

# Use a local GMT file
local_gmt_path <- tempfile(fileext = ".gmt")
writeLines(c(
  "GENESET1\tdesc1\tGENE1\tGENE2",
  "GENESET2\tdesc2\tGENE3\tGENE4"
), local_gmt_path)
FetchGMT(local_gmt_path)
# Returns: /tmp/RtmpXXXXXX/filexxxx.gmt (original path)
}
# }
```
