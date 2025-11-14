# Compute statistic p-value by given numbers

Compute statistic p-value by given numbers

## Usage

``` r
ComputePValue(k, M, n, N, method = c("fisher", "hypergeometric"))
```

## Arguments

- k:

  Number of genes in both user list and gene set

- M:

  Number of genes in the gene set

- n:

  Number of genes in the user list

- N:

  Total number of genes in the universe

- method:

  Method for computing p-value, either "fisher" or "hypergeometric"

## Value

p-value from the specified test
