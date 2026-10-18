# Run individual lines after run_analysis.R.
# These commands display existing results without refitting any model.

# Main tables
Main_tables$Baseline
Main_tables$Categorical_OR

# Supplementary tables (NHANES only)
Supplement_tables$Variable_definitions
Supplement_tables$Threshold
Supplement_tables$Joint16_OR
Supplement_tables$Components_OR
Supplement_tables$Mediation

# Main figures: standalone objects
print(p_PA_ST_categories)
print(p_PA_RCS)
print(p_ST_RCS)
print(p_Fig2d_interaction)
print(p_Fig3a_joint_share)
print(p_Fig3b_joint_prevalence)
print(p_Fig3c_joint_OR)
print(p_Fig3d_joint_forest)
print(p_Fig4a_ASBI_distribution)
print(p_ASBI_categories)
print(p_ASBI_RCS)

# Supplementary ASBI sensitivity figures
print(p_S_ASBI_Full_c0p01)
print(p_S_ASBI_Full_c0p1)
print(p_S_ASBI_Full_c0p5)
print(p_S_ASBI_Exclude_P1S1_c0p01)
print(p_S_ASBI_Exclude_P1S1_c0p1)
print(p_S_ASBI_Exclude_P1S1_c0p5)

# Supplementary PA composition figures
print(p_S_composition_Domain)
print(p_S_composition_Intensity)
print(p_S_composition_DomainIntensity)

# Supplementary mediation figures
print(p_S_mediation_BMI)
print(p_S_mediation_Waist)
print(p_S_mediation_Hypertension)
print(p_S_mediation_Diabetes)
print(p_S_mediation_Hyperlipidemia)
print(p_S_mediation_Uric_acid)
print(p_S_mediation_summary)

# Numerical results behind the figures
Manuscript_results$main$tables$RCS_summary
Manuscript_results$main$tables$Interaction_tests
Manuscript_results$interaction$tables$Test
Manuscript_results$sensitivity$tables$RCS_summary
Manuscript_results$components$tables$Composition
Manuscript_results$mediation$tables$Signed_indirect_total

# Optional: two Excel workbooks, no images or automatic export.
# Excel_files <- export_manuscript_tables(Manuscript_results, directory="tables")
