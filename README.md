# EWAScaller

EWAScaller is an R client for the [EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/) database. It queries CpG probes, gene symbols, and genomic regions for epigenome-wide association study (EWAS) results, runs trait/GO/KEGG/genomic-location enrichment analyses on probe sets, finds expression quantitative trait methylation (eQTM) pairs, and provides customisable plots including word clouds of enriched trait terms.

Supported arrays: Illumina 450K, EPIC (850K), and EPIC v2, including
non-CpG `ch.*` probes.

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
# enrichment needs 20-5000 distinct probes
probes <- c(
  "cg05575921", "cg00240195", "cg16419584", "cg11903855", "cg21566642",
  "cg23771366", "cg19859270", "cg00382138", "cg25648203", "cg01940273",
  "cg23576855", "cg00401753", "cg00300637", "cg12803068", "cg09935388",
  "cg01899089", "cg04180046", "cg00050873", "cg02228185", "cg18147296"
)

enr <- ewas_enrichment(probes, background = "850K",
                        types = c("trait", "genomic_location", "gene_ontology", "kegg"))
summary(enr)
top_traits(enr)
```

Which background matters depends on the analysis: in our tests, trait
enrichment gave identical results with the 450K and 850K backgrounds,
while genomic location, Gene Ontology, and KEGG results changed.

### EPIC v2 arrays

EPIC v2 probe names such as `cg05575921_BC21` are accepted everywhere and
converted to the 450K/EPIC v1 names EWAS Atlas uses (replicate probes are
merged). For enrichment and eQTM, use `background = "EPICv2"`:

```r
enr <- ewas_enrichment(my_epicv2_probes, background = "EPICv2")
```

Things to know:

- EWAS Atlas has no EPIC v2 background, so `"EPICv2"` uses the EPIC/850K
  background. That background also contains EPIC v1 probes that are not on
  EPIC v2, so genomic location, GO, and KEGG p-values are approximate
  (trait enrichment is unaffected).
- Probes that are only on EPIC v2 (about 200,000 of them) have no EWAS
  Atlas data, and the toolkit would reject the whole job if the list
  contained one. EWAScaller therefore drops them automatically before
  submitting, using a bundled list of the 721,378 probes on both EPIC v1
  and EPIC v2, and tells you how many it dropped. They are kept in
  `enr$excluded` (or `eq$excluded`). At least 20 probes must remain.

### eQTM: CpG methylation vs gene expression

`ewas_eqtm()` finds expression quantitative trait methylation (eQTM)
pairs: CpG-gene pairs whose methylation correlates with gene expression.
The correlations were computed by EWAS Atlas from matched methylation and
RNA-seq data in The Cancer Genome Atlas (TCGA), i.e. tumour tissue, for
six tissues (brain, colon, kidney, liver, stomach, testis). Keep that in
mind when interpreting results for blood or healthy tissue. It takes the
same 20-5000 probe set as `ewas_enrichment()`.

```r
eq <- ewas_eqtm(probes, background = "850K")
summary(eq)
head(eq$pairs)

# per-sample methylation and expression behind one pair
eqtm_samples("cg05575921", "ENSG00000063438", "brain")
```

### Plots

```r
plot_wordcloud(cpg_res)
plot_wordcloud(enr, type = "trait")
plot_wordcloud(eq)

autoplot(cpg_res, type = "traits")
autoplot(cpg_res, type = "chromosome")
autoplot(cpg_res, type = "direction")
autoplot(enr, type = "trait")
autoplot(enr, type = "kegg", style = "bar")

autoplot(eq, type = "heatmap")   # probes x genes, coloured by correlation
autoplot(eq, type = "volcano")   # correlation vs -log10(p)
autoplot(eq, type = "tissue")    # pairs per tissue
autoplot(eq, type = "scatter", probe = "cg05575921", gene = "AHRR", tissue = "brain")
```

Every plot accepts the same customisation arguments and returns a regular
`ggplot` object, so you can keep adding layers with `+`:

```r
autoplot(enr,
  colors = "#6d597a",              # colours (meaning depends on the plot)
  title = "Smoking-related CpGs",  # NA removes the title
  subtitle = "850K background",
  theme = ggplot2::theme_bw,       # any theme function or theme object
  base_size = 14,                  # font size
  wrap = 30                        # wrap long labels at 30 characters
)
```

Word clouds take `title`, `subtitle`, `base_size`, `theme`, `max_size`
(largest font; scaled automatically when there are few words), and `seed`
(layout seed).

## Data source

EWAScaller queries the EWAS Atlas REST and toolkit web services
(<https://ngdc.cncb.ac.cn/ewas/api>). Please cite the EWAS Atlas when using
data retrieved through this package.

## License

MIT
