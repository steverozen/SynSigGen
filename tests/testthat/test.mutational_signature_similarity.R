context("Mutational signature matching")

test_that("MatchSigs2Directions", {
  load("data.for.matchsigs.Rdata")
  expect_equal(
    MatchSigs2Directions(.test.extracted.sigs, .test.ground.truth.sigs),
    .test.match.out)
})
