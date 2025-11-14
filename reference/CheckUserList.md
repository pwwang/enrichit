# Check user list against gene sets

Check user list against gene sets

## Usage

``` r
CheckUserList(userlist, genesets, use_matched_only, show_all_unmatched = FALSE)
```

## Arguments

- userlist:

  Vector of user-provided genes

- genesets:

  List of gene sets

- use_matched_only:

  Logical, whether to return only matched genes Otherwise the genes in
  the user list will be all used. This will affect the number of genes
  in the user list when computing the p-value.

- show_all_unmatched:

  Logical, whether to show all unmatched genes. Or instead show only the
  first 10 unmatched genes.

## Value

A vector of user-provided genes that are also in the gene sets
