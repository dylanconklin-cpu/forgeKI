hdr_stage10_mock_stage9 <- function(cfg, score = 84) {
  locus <- list(gene_symbol = cfg$gene, transcript_id = "tx1")
  rec <- tibble::tibble(
    Design_Rank = 1L,
    Guide_ID = "g001",
    Stage2_Rank = 1L,
    Final_Design_Score = score,
    Recommendation_Tier = "BACKUP_candidate",
    Recommendation_Status = "WARN_backup_candidate"
  )
  summ <- tibble::tibble(
    N_Designs_Scored = 1L,
    N_Recommended_Primary = 0L,
    N_Backup_Candidates = 1L,
    N_Failed_Candidates = 0L,
    Top_Guide_ID = "g001",
    Top_Final_Design_Score = score,
    Stage7_QC_Status = "PASS_virtual_allele_validated",
    Stage8_QC_Status = "PASS_donor_modules_constructed",
    N_Total_Donor_Edits = 0L,
    Stage9_QC_Status = "WARN_no_primary_recommendation"
  )
  out <- list(stage = "stage9_design_scoring", schema_version = 1L, cfg = cfg, locus = locus, design_recommendations = rec, recommendation_summary = summ)
  class(out) <- c("hdr_stage9_result", "list")
  out
}
