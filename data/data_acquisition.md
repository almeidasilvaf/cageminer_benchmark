Data acquisition
================

# Overview

Here, we describe the code to create the input data for
*[cageminer](https://bioconductor.org/packages/3.15/cageminer)* and
[Camoco](https://camoco.readthedocs.io/en/latest/index.html), which are
stored in `cageminer_input/` and `camoco_input/`, respectively.

``` r
# Load required packages
library(here)
library(BioNERO)
library(tidyverse)
library(SummarizedExperiment)
```

## cageminer_input

To mine candidate genes with
*[cageminer](https://bioconductor.org/packages/3.15/cageminer)*, we will
need:

1.  Gene coordinates in a `GRanges` object - `gene_ranges.rda`;
2.  SNP positions in a `GRanges` object - `snp_positions.rda`;
3.  Gene expression matrix and sample metadata in a
    `SummarizedExperiment` object - `gmax_se.rda`;
4.  Guide (a.k.a. “reference”) genes as a character vector -
    `guides.rda`.
5.  A gene coexpression network inferred with
    *[BioNERO](https://bioconductor.org/packages/3.15/BioNERO)* -
    inferred on-the-fly.

Below, you can find the code to create each of these objects.

### gene_ranges.rda

``` r
gene_ranges <- rtracklayer::import(
    "https://raw.githubusercontent.com/almeidasilvaf/SoyFungi_GWAS_GCN/main/data/PLAZA_selected.transcripts.gff.gz"
)

# Clean file to remove unnecessary columns and keep only genes
gene_ranges <- gene_ranges[gene_ranges$type == "gene"]
gene_ranges$source <- NULL
gene_ranges$score <- NULL
gene_ranges$phase <- NULL
gene_ranges$pacid <- NULL
gene_ranges$pid <- NULL
gene_ranges$tid <- NULL
gene_ranges$old_id <- NULL
gene_ranges$old_tid <- NULL
gene_ranges$UniProtKB <- NULL
gene_ranges$Parent <- NULL

# Save file
save(
    gene_ranges,
    file = here("data", "cageminer_input", "gene_ranges.rda"),
    compress = "xz"
)
```

### snp_positions.rda

``` r
# Download file
snp_file <- file.path(tempdir(), "snp_file.rda")
download.file(
    "https://raw.githubusercontent.com/almeidasilvaf/SoyFungi_GWAS_GCN/main/products/result_files/snp_granges.rda",
    destfile = snp_file
)

# Load file
load(snp_file)
snp_positions <- snp_grangeslist$Fgraminearum

# Save file
save(
    snp_positions,
    file = here("data", "cageminer_input", "snp_positions.rda")
)
```

### gmax_se.rda

``` r
# Load whole expression data and filter to only contain fungi-infected soybean
load("~/Dropbox/Atlas/atlasv2_tpm.rda")
exp <- atlas_tpm[, atlas_tpm$Stress_info == "fungus" &
                   !is.na(atlas_tpm$Stress_info)]
rm(atlas_tpm)

#----Filter and preprocess the SE object----
gmax_se <- exp_preprocess(exp, min_exp = 5, Zk_filtering = FALSE)
metadata_fungi <- colData(gmax_se) %>%
  as.data.frame() %>%
  mutate(Pathogen = str_sub(Pathogen, 1, 3)) %>%
  mutate(annot = paste(Pathogen, Sample_description, sep = "_")) %>%
  select(annot)
colData(gmax_se) <- DataFrame(metadata_fungi)

# Save object
save(
    gmax_se,
    file = here::here("data", "cageminer_input", "gmax_se.rda"),
    compress = "xz"
)
```

### guides.rda

``` r
# Get data as a character vector
guides <- read.csv("https://raw.githubusercontent.com/almeidasilvaf/SoyFungi_GWAS_GCN/main/products/tables/sup_table3.tsv", header = TRUE, sep = "\t")
guides <- guides$Gene

# Save
save(
    guides,
    file = here("data", "cageminer_input", "guides.rda")
)
```

## camoco_input

To mine candidate genes with
[Camoco](https://camoco.readthedocs.io/en/latest/index.html), we will
need:

1.  Gene annotation in a .GFF file - `annotation.gff`;
2.  Gene expression in a tab-separated file with genes in row names and
    samples in column names - `exp_tpm.tsv`;
3.  SNP positions in a tab-separated file - `gmax_fungi_snps.tsv` - with
    the fields:
    -   **Trait**: Trait name
    -   **CHR:** Chromosome number
    -   **POS:** SNP position

Below, you can find the code to create each of these files.

### annotation.gff

``` r
# Load annotationa and filter it
annotation <- rtracklayer::import(
    "https://raw.githubusercontent.com/almeidasilvaf/SoyFungi_GWAS_GCN/main/data/PLAZA_selected.transcripts.gff.gz"
)
annotation$old_id <- NULL
annotation$old_tid <- NULL
annotation <- annotation[annotation$type == "gene"]
annotation$UniProtKB <- NULL

# Export as GFF
rtracklayer::export.gff3(
    annotation, 
    con = here("data", "camoco_input", "annotation.gff3")
)
```

### exp_tpm.tsv

``` r
# Load whole expression data and filter to only contain fungi-infected soybean
load("~/Dropbox/Atlas/atlasv2_tpm.rda")
exp_tpm <- atlas_tpm[, atlas_tpm$Stress_info == "fungus" &
                   !is.na(atlas_tpm$Stress_info)]
rm(atlas_tpm)

# Get expression matrix
exp_tpm <- assay(exp_tpm)

# Remove genes with median expression < 5 (Camoco default)
exp_tpm_final <- exp_tpm[matrixStats::rowMedians(as.matrix(exp_tpm)) >= 5, ]

# Save as a tab-delimited file
write.table(
    exp_tpm_final, quote = FALSE, sep = "\t",
    file = here::here("data", "camoco_input", "exp_tpm.tsv")
)
```

### gmax_fungi_snps.tsv

``` r
# Load SNP positions 
load(here("data", "cageminer_input", "snp_positions.rda"))

# Create a table with columns "Trait", "CHR", and "POS"
gmax_fungi_snps <- as.data.frame(snp_positions)
gmax_fungi_snps <- gmax_fungi_snps[, c("seqnames", "start", "Trait")]
names(gmax_fungi_snps) <- c("CHR", "POS", "Trait")

# Save as a tab-delimited file
readr::write_tsv(
    gmax_fungi_snps,
    file = here("data", "camoco_input", "gmax_fungi_snps.tsv")
)
```

# Session information

This document was created under the following conditions:

``` r
sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.2.1 (2022-06-23)
    ##  os       Ubuntu 20.04.4 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       Europe/Brussels
    ##  date     2022-07-11
    ##  pandoc   2.17.1.1 @ /usr/lib/rstudio/bin/quarto/bin/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package     * version date (UTC) lib source
    ##  BiocManager   1.30.18 2022-05-18 [1] CRAN (R 4.2.0)
    ##  BiocStyle     2.25.0  2022-06-15 [1] Github (Bioconductor/BiocStyle@7150c28)
    ##  cli           3.3.0   2022-04-25 [1] CRAN (R 4.2.0)
    ##  digest        0.6.29  2021-12-01 [1] CRAN (R 4.2.0)
    ##  evaluate      0.15    2022-02-18 [1] CRAN (R 4.2.0)
    ##  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
    ##  htmltools     0.5.2   2021-08-25 [1] CRAN (R 4.2.0)
    ##  knitr         1.39    2022-04-26 [1] CRAN (R 4.2.0)
    ##  magrittr      2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
    ##  rlang         1.0.3   2022-06-27 [1] CRAN (R 4.2.1)
    ##  rmarkdown     2.14    2022-04-25 [1] CRAN (R 4.2.0)
    ##  rstudioapi    0.13    2020-11-12 [1] CRAN (R 4.2.0)
    ##  sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
    ##  stringi       1.7.6   2021-11-29 [1] CRAN (R 4.2.0)
    ##  stringr       1.4.0   2019-02-10 [1] CRAN (R 4.2.0)
    ##  xfun          0.31    2022-05-10 [1] CRAN (R 4.2.0)
    ##  yaml          2.3.5   2022-02-21 [1] CRAN (R 4.2.0)
    ## 
    ##  [1] /home/faalm/R/x86_64-pc-linux-gnu-library/4.2
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
