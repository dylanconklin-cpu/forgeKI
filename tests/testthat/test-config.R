test_that("hdr_config creates an inert validated object", {
  cfg <- hdr_config(gene = "ACTB", project_dir = tempdir())
  expect_s3_class(cfg, "hdr_config")
  expect_equal(cfg$gene, "ACTB")
  expect_equal(cfg$cassette_id, "HiBiT_dTAG_GFP_EF1a_BSD_P2A_mKATE")
})

test_that("run_hdr_pipeline is the controlled package orchestrator", {
  expect_true(is.function(run_hdr_pipeline))
})
