# Perform enrichment analysis

Perform enrichment analysis

## Usage

``` r
EnrichIt(
  userlist,
  dbs,
  method = ifelse(tolower(style[1]) == "enrichr", "fisher", "hypergeometric"),
  use_matched_only = ifelse(tolower(style[1]) == "enrichr", FALSE, TRUE),
  padjust_method = c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY", "fdr"),
  background = NULL,
  style = c("enrichr", "Enrichr", "clusterProfiler", "ClusterProfiler",
    "clusterprofiler"),
  return_all = FALSE
)
```

## Arguments

- userlist:

  Vector of user-provided genes

- dbs:

  List of gene sets or paths to GMT files It can be a vector of database
  names. You can set the names of the vector, which will be used as the
  database names. Otherwise a list is expected, where each element is a
  database (e.g. parsed from a gmt file). If a database is given
  directly (without a name), the expression of this argument will be
  used as the database name.

- method:

  Method for computing p-value, either "fisher" or "hypergeometric" When
  `style` is "enrichr", the method defaults to "fisher". When `style` is
  "clusterProfiler", the method defaults to "hypergeometric".

- use_matched_only:

  Logical, whether to use only matched genes against the gene sets. This
  will affect the number of genes in the user list when computing the
  p-value. By default, when `style` is "enrichr", this is set to FALSE.
  When `style` is "clusterProfiler", this is set to TRUE.

- padjust_method:

  Method for adjusting p-values, either "BH", "bonferroni", "holm",
  "hochberg", "hommel", "BY", "fdr"

- background:

  Vector of all genes in the universe or a number of genes in the
  universe. If NULL, the number of genes in the gene set will be used.
  For "enrichr", the default is 20,000. For "clusterProfiler", the
  default is the number of unique genes in the gene set. Note that for
  "enrichr", if a vector is given, the length of it will be used, no
  checking will be done to see if userlist and genes from dbs are in the
  vector.

- style:

  Style of the output, either "enrichr" or "clusterProfiler"

- return_all:

  Logical, whether to return all results (all gene sets in dbs) or only
  those with at least one gene in the user list.

## Value

A data frame with the results. When `style` is "enrichr", the columns
are:

- Database: Name of the database

- Term: Name of the term

- Overlap: Number of genes in the user list and the gene set

- P.value: p-value from the enrichment test

- Adjusted.P.value: Adjusted p-value from the enrichment test

- Odds.Ratio: Odds ratio from the enrichment test

- Combined.Score: Combined score from the enrichment test

- Genes: Genes in the user list that are also in the gene set

- Rank: Rank of the term based on the combined score

When `style` is "clusterProfiler", the columns are:

- ID: Name of the term

- Description: Description of the term

- GeneRatio: Ratio of genes in the user list and the gene set

- BgRatio: Ratio of genes in the gene set and the universe

- Count: Number of genes in the user list that are also in the gene set

- pvalue: p-value from the enrichment test

- p.adjust: Adjusted p-value from the enrichment test

- qvalue: Q-value from the enrichment test

- geneID: Genes in the user list that are also in the gene set

- Database: Name of the database
