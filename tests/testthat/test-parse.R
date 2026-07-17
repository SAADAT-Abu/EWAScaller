single_probe_json <- '{
  "probeId": "cg05575921",
  "chrHg19": 5,
  "posHg19": 373378,
  "cpgIsland": "Shelf",
  "relatedTranscription": [
    {"geneName": "AHRR", "posToTss": 69087, "ensemblTranscriptId": "ENST1"},
    {"geneName": "AHRR", "posToTss": 69079, "ensemblTranscriptId": "ENST2"}
  ],
  "associationList": [
    {"studyId": "ES00580", "trait": "maternal smoking", "correlation": "neg", "rank": 1, "pmid": 29618728},
    {"studyId": "ES00084", "trait": "smoking", "correlation": "neg", "rank": null, "pmid": 23070629}
  ]
}'

multi_probe_json <- '[
  {
    "probeId": "cg00300637",
    "chrHg19": 5,
    "posHg19": 319433,
    "cpgIsland": "Shore",
    "relatedTranscription": [{"geneName": "AHRR", "posToTss": 15142, "ensemblTranscriptId": "ENST1"}],
    "associationList": [{"studyId": "ES00118", "trait": "smoking", "correlation": "pos", "rank": 205, "pmid": 27651444}]
  },
  {
    "probeId": "cg00401753",
    "chrHg19": 5,
    "posHg19": 405986,
    "cpgIsland": "Other",
    "relatedTranscription": [],
    "associationList": []
  }
]'

test_that(".rows_from_probe_data handles a single-probe object payload", {
  data <- jsonlite::fromJSON(single_probe_json, simplifyVector = TRUE)
  parsed <- EWAScaller:::.rows_from_probe_data(data)

  expect_s3_class(parsed$probes, "data.frame")
  expect_equal(nrow(parsed$probes), 1)
  expect_equal(parsed$probes$probeId, "cg05575921")
  expect_equal(parsed$probes$chrHg19, "5")
  expect_equal(parsed$probes$gene, "AHRR")

  expect_s3_class(parsed$associations, "data.frame")
  expect_equal(nrow(parsed$associations), 2)
  expect_setequal(parsed$associations$trait, c("maternal smoking", "smoking"))
  expect_true(is.na(parsed$associations$rank[parsed$associations$trait == "smoking"]))
})

test_that(".rows_from_probe_data handles a multi-probe array payload", {
  data <- jsonlite::fromJSON(multi_probe_json, simplifyVector = TRUE)
  parsed <- EWAScaller:::.rows_from_probe_data(data)

  expect_equal(nrow(parsed$probes), 2)
  expect_setequal(parsed$probes$probeId, c("cg00300637", "cg00401753"))
  expect_equal(nrow(parsed$associations), 1)
  expect_equal(parsed$associations$probeId, "cg00300637")
})

test_that(".rows_from_probe_data returns NULLs for empty input", {
  parsed <- EWAScaller:::.rows_from_probe_data(NULL)
  expect_null(parsed$probes)
  expect_null(parsed$associations)
})
