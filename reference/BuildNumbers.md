# Build the numbers for the tests

Build the numbers for the tests

## Usage

``` r
BuildNumbers(userlist, geneset, background)
```

## Arguments

- userlist:

  Vector of user-provided genes

- geneset:

  Vector of genes in the gene set

- background:

  Vector of all genes in the universe or a number of genes in the
  universe

## Value

A list with the following elements:

- k: Number of genes in the user list that are also in the gene set

- M: Number of genes in the gene set

- n: Number of genes in the user list

- N: Total number of genes in the universe
