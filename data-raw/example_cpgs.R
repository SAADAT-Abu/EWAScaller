## Build the bundled example dataset from a small set of test CpG probes.
example_cpgs <- data.frame(
  query_cpg = c(
    rep("cg11903855", 6L),
    rep("cg00240195", 16L),
    rep("cg16419584", 8L)
  ),
  stringsAsFactors = FALSE
)

usethis::use_data(example_cpgs, overwrite = TRUE)
