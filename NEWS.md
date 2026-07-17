# EWAScaller 0.1.0

* Initial release.
* `query_cpg()`, `query_gene()`, `query_region()` query the EWAS Atlas REST
  API with bounded-concurrency, rate-limited requests.
* `ewas_enrichment()` runs trait, genomic location, Gene Ontology, and KEGG
  pathway enrichment analyses via the EWAS Atlas toolkit.
* `plot_wordcloud()` and `autoplot()` methods for visualising query and
  enrichment results.
