## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The only NOTE is the expected "New submission" note for a first-time
package submission.

## Test environments

* Fedora Linux 44, R 4.6.0 (local, `R CMD check --as-cran`)

## Notes for reviewers

* Functions that contact the EWAS Atlas web service
  (<https://ngdc.cncb.ac.cn/ewas/>) are wrapped in `\donttest{}` in
  examples and skipped via `testthat::skip_on_cran()` /
  `skip_if_not(reachable)` in tests, so the check suite does not depend on
  network access or CRAN's test machines being able to reach the service.
* The vignette's live-API code chunks are set to `eval = FALSE` for the
  same reason; the vignette documents usage without requiring network
  access to build.
