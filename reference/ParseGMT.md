# Read GMT file

Read GMT file

## Usage

``` r
ParseGMT(gmtfile, swap_name_desc_if_needed = TRUE)
```

## Arguments

- gmtfile:

  Path to the GMT file

- swap_name_desc_if_needed:

  Logical, whether to swap name and description fields. They will be
  swapped only if:

  - `swap_name_desc_if_needed` is `TRUE`; and

  - The descriptions are not empty; and

  - The descriptions are shorter than the names; and

  - The descriptions are not ID-like (i.e., hsa00001, or 123456).

## Value

A list of gene sets, the names of the list are the gene set names, and
each element is a vector of gene names.

## Examples

``` r
# \donttest{
if (FALSE) {
# Example GMT file content
gmtfile <- system.file("extdata", "KEGG_2021_Human.gmt", package = "enrichit")
# Read the GMT file
ParseGMT(gmtfile)
}
# }
```
