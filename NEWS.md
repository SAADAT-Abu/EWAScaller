# EWAScaller 0.2.0

* New `ewas_eqtm()` finds expression quantitative trait methylation (eQTM)
  pairs: CpG-gene correlations between methylation and expression in TCGA
  tissue data, via the EWAS Atlas toolkit's "Expression Regulation"
  analysis. `eqtm_samples()` downloads the per-sample data behind one pair.
* New eQTM plots: `autoplot()` with `type = "heatmap"`, `"volcano"`,
  `"tissue"`, or `"scatter"`, and `plot_wordcloud()` of eQTM genes.
* EPIC v2 support: EPIC v2 probe names (e.g. `cg05575921_BC21`) are
  converted to their 450K/EPIC v1 names everywhere, and enrichment/eQTM
  accept `background = "EPICv2"`, which uses the EPIC/850K background as an
  approximation because EWAS Atlas has no EPIC v2 background. Probes that
  are only on EPIC v2 (which the toolkit would reject) are dropped
  automatically using a bundled list of probes on both EPIC v1 and v2, and
  reported in the result's new `excluded` element. Input probes missing
  from a custom background are dropped the same way.
* Non-CpG `ch.*` probe names (e.g. `ch.4.37949012R`) are no longer
  rejected as invalid.
* Word clouds with few words now use a larger font by default.
* `autoplot()` for enrichment results gains `style = "bar"`.
* All plots can be customised with `colors`, `title`, `subtitle`,
  `base_size`, `theme`, and `wrap` (word clouds: `max_size` and `seed`).
  Defaults keep the previous appearance.
* `ewas_enrichment()` now uses the 450K background when `background` is
  left at its default (previously this errored).
* `plot_wordcloud()` no longer changes the global random seed.

# EWAScaller 0.1.0

* Initial release.
* `query_cpg()`, `query_gene()`, `query_region()` query the EWAS Atlas REST
  API with bounded-concurrency, rate-limited requests.
* `ewas_enrichment()` runs trait, genomic location, Gene Ontology, and KEGG
  pathway enrichment analyses via the EWAS Atlas toolkit.
* `plot_wordcloud()` and `autoplot()` methods for visualising query and
  enrichment results.
