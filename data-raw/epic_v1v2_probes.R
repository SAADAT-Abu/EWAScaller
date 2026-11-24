## Build the internal list of probes present on both the Illumina EPIC v1
## (850K) and EPIC v2 arrays, stored in R/sysdata.rda.
##
## EWAS Atlas has no EPIC v2 background, so `background = "EPICv2"` runs
## against the EPIC/850K background. The toolkit rejects the whole job if
## any input probe is not in that background, so probes that are only on
## EPIC v2 are dropped client-side using this list.
##
## Requires the Bioconductor annotation packages below (not package
## dependencies; only needed to regenerate the data).

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

v1 <- rownames(minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, what = "Locations"))
v2 <- rownames(minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38, what = "Locations"))
# EPIC v2 names carry a strand/design/replicate suffix, e.g. cg05575921_BC21.
v2 <- unique(sub("_(BC|TC|BO|TO)[0-9]+$", "", v2))

epic_v1v2_probes <- sort(intersect(v1, v2))
message(length(epic_v1v2_probes), " probes on both EPIC v1 and EPIC v2")

usethis::use_data(epic_v1v2_probes, internal = TRUE, overwrite = TRUE, compress = "xz")
