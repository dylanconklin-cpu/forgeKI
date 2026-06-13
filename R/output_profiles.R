# Output profile handling.

redacted_columns_for_profile <- function(profile) {
  switch(profile, full_internal = character(), user_facing = c("Internal_Model_Notes", "Private_Scoring_Equation"), collaborator = c("Internal_Model_Notes", "Private_Scoring_Equation"), public = c("Internal_Model_Notes", "Private_Scoring_Equation", "Private_Feature_Matrix_Path", "DepMap_Private_Feature_Detail"), character())
}

forgeki_user_facing_output_names <- function() {
  c(
    detailed_html = "forgeki_report.html",
    executive_summary_html = "forgeki_executive_summary.html",
    order_csv = "forgeki_order_sheet.csv",
    report_model_json = "report_model.json",
    report_model_rds = "report_model.rds"
  )
}
