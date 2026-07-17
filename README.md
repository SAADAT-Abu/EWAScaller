# EWAScaller

EWAScaller is an R client for the [EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/) database. It queries CpG probes, gene symbols, and genomic regions for epigenome-wide association study (EWAS) results, runs trait/GO/KEGG/genomic-location enrichment analyses on probe sets, and provides plotting helpers including word clouds of enriched trait terms.

## Installation

```r
# install.packages("remotes")
remotes::install_github("SAADAT-Abu/EWAScaller")
```

## Usage

### Query by CpG probe, gene, or genomic region

```r
library(EWAScaller)

cpg_res <- query_cpg(c("cg05575921", "cg11903855", "cg00240195"))
summary(cpg_res)
head(cpg_res$associations)

gene_res <- query_gene("AHRR")
region_res <- query_region(chr = "5", start = 373000, end = 374000)
```

All three query functions dispatch requests across a small worker pool
(`workers`) and pause `delay` seconds between batches, so large probe lists
can be resolved faster without overwhelming the remote service.

### Enrichment analysis

```r
data(example_cpgs)
probes <- unique(example_cpgs$query_cpg)

enr <- ewas_enrichment(probes, background = "850K",
                        types = c("trait", "genomic_location", "gene_ontology", "kegg"))
summary(enr)
top_traits(enr)
```

### Plots

```r
plot_wordcloud(cpg_res)
plot_wordcloud(enr, type = "trait")

autoplot(cpg_res, type = "traits")
autoplot(cpg_res, type = "chromosome")
autoplot(cpg_res, type = "direction")
autoplot(enr, type = "trait")
```

## Data source

EWAScaller queries the EWAS Atlas REST and toolkit web services
(<https://ngdc.cncb.ac.cn/ewas/api>). Please cite the EWAS Atlas when using
data retrieved through this package.

## License

MIT
