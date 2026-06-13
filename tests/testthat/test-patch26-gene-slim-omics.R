test_that("patch26 creates gene-slim omics bundles and maps local long-table schemas", {
  td <- tempfile("forgeki_patch26_gene_slim_"); dir.create(td)
  global <- data.frame(
    depmap_id = c("ACH-1", "ACH-2", "ACH-3"),
    cell_line_name = c("A", "B", "C"),
    lineage = c("Lung", "Breast", "CNS"),
    GeneContext_Score = c(95, 80, 70),
    HDR_Recommendation_Rank = c(1L, 2L, 3L)
  )
  expr <- data.frame(depmap_id = c("ACH-1", "ACH-2", "ACH-3", "ACH-1"), gene = c("ACTB", "ACTB", "TP53", "ACTB"), rna_expression = c(10, 0.5, 3, 12), entrez_id = c(60, 60, 7157, 60), gene_name = c("ACTB", "ACTB", "TP53", "ACTB"), cell_line = c("A", "B", "C", "A"))
  cn <- data.frame(depmap_id = c("ACH-1", "ACH-2", "ACH-3"), gene = c("ACTB", "ACTB", "TP53"), log_copy_number = c(2.1, 1.2, 3.0), entrez_id = c(60, 60, 7157), gene_name = c("ACTB", "ACTB", "TP53"), cell_line = c("A", "B", "C"))
  crispr <- data.frame(depmap_id = c("ACH-1", "ACH-2", "ACH-3"), gene = c("ACTB", "ACTB", "TP53"), dependency = c(-0.9, -0.1, -0.5), entrez_id = c(60, 60, 7157), gene_name = c("ACTB", "ACTB", "TP53"), cell_line = c("A", "B", "C"))
  mut <- data.frame(depmap_id = c("ACH-2", "ACH-3"), gene_name = c("ACTB", "TP53"), var_annotation = c("missense", "nonsense"), is_deleterious = c(FALSE, TRUE))
  fus <- data.frame(ModelID = c("ACH-1", "ACH-3"), Gene1 = c("X", "TP53"), Gene2 = c("ACTB", "Y"), CanonicalFusionName = c("X--ACTB", "TP53--Y"), FFPM = c(1.2, 0.2))
  rrbs_tss <- data.frame(gene = c("ACTB", "TP53"), DMS53_LUNG = c(0.2, 0.8), SW1116_LARGE_INTESTINE = c(0.4, 0.1), check.names = FALSE)
  rrbs_cpg <- data.frame(gene_name = c("ACTB", "ACTB", "TP53"), DMS53_LUNG = c(0.3, 0.5, 0.8), SW1116_LARGE_INTESTINE = c(0.4, 0.6, 0.2), check.names = FALSE)
  designs <- data.frame(Design_ID = "DESIGN_001_g001", Guide_ID = "g001", Design_Rank = 1L, Final_Design_Score = 90, Recommendation_Status = "PASS_recommended_for_production")

  gp <- file.path(td, "global.csv"); ep <- file.path(td, "expr.csv"); cp <- file.path(td, "cn.csv"); crp <- file.path(td, "crispr.csv"); mp <- file.path(td, "mut.csv"); fp <- file.path(td, "fusion.csv"); tp <- file.path(td, "rrbs_tss.csv"); cpgp <- file.path(td, "rrbs_cpg.csv"); dp <- file.path(td, "designs.csv")
  write.csv(global, gp, row.names = FALSE); write.csv(expr, ep, row.names = FALSE); write.csv(cn, cp, row.names = FALSE); write.csv(crispr, crp, row.names = FALSE); write.csv(mut, mp, row.names = FALSE); write.csv(fus, fp, row.names = FALSE); write.csv(rrbs_tss, tp, row.names = FALSE); write.csv(rrbs_cpg, cpgp, row.names = FALSE); write.csv(designs, dp, row.names = FALSE)

  bundle_path <- file.path(td, "full_bundle.rds")
  forgeki_compile_stage10_omics_bundle(output_rds = bundle_path, global_ranking_path = gp, expression_path = ep, copy_number_path = cp, crispr_dependency_path = crp, mutation_path = mp, fusion_path = fp, rrbs_tss_path = tp, rrbs_cpg_path = cpgp, release_label = "patch26_toy", compress = "gzip")
  slim_path <- forgeki_make_gene_slim_stage10_omics_bundle(bundle_path, gene = "ACTB", output_dir = file.path(td, "slim"), verbose = FALSE)
  slim <- forgeki_load_stage10_omics_bundle(slim_path)
  expect_lt(nrow(slim$tables$expression_path), nrow(expr))
  expect_equal(unique(slim$tables$expression_path$gene), "ACTB")
  expect_true(file.exists(file.path(td, "slim", "ACTB_gene_slim_bundle_table_summary.csv")))

  out <- forgeki_build_stage10_reference(gene = "ACTB", output_dir = file.path(td, "out"), omics_bundle_path = bundle_path, design_table_path = dp, module_label = "test_modules", mode = "internal", build_10a = TRUE, build_10b = TRUE, build_10c = TRUE, build_10d = TRUE, build_10e = TRUE, top_n = 2L)
  expect_true(out$builder_qc$Stage10A_Context_Constructed[[1]])
  expect_true(all(c("Target_Gene_Expression", "Target_Gene_Copy_Number", "Target_Gene_Dependency", "Target_Gene_Mutation_Status", "Target_Gene_Fusion_Status") %in% names(out$stage10a_context)))
  expect_true(all(c("RNA_expression", "copy_number", "CRISPR_dependency", "mutation", "fusion") %in% out$stage10a_gene_feature_schema_audit$Feature_Source))
  expect_true(any(out$stage10a_gene_feature_schema_audit$Schema_Status == "PASS_gene_feature_schema_mapped"))
  expect_true(any(out$stage10a_feature_status$Feature_Source == "RNA_expression" & out$stage10a_feature_status$Feature_Status == "PASS_feature_loaded"))
  expect_true(file.exists(out$output_paths$stage10a_gene_feature_schema_audit))
  expect_true(out$builder_qc$Stage10E_Final_Ranking_Constructed[[1]])
})
