This folder contains all BigQuery scripts needed to build the cohorts for the ARIVe project.
You will need to adjust some of the table pointers so that they are specific to your bigquery account.

For MIMIC-IV, run:
1) MIMIC_O2deliveryFiO2.sql
2) ARIVe_MIMIC_eligibility.sql
3) ARIVE_MIMIC_timevarying.sql
4) ARIVE_MIMIC_baseline.sql

Then download the resulting tables:
ARIVe_MIMIC_Eligibility
ARIVe_MIMIC_Baseline
ARIVe_MIMIC_Timevarying

For eICU, run:
1) ARIVe_eICU_oxygenFlow.sql
2) ARIVe_eICU_oxygenDevices.sql
3) ARIVe_eICU_eligibility.sql
4) ARIVe_eICU_timevarying.sql
5) ARIVe_eICU_baseline.sql
6) ARIVe_eICU_hospitalCoefficient.sql

Note that eICU oxygen device classifications are available in the spreadsheet "Oxygen device classifications.xlsx"

Then download the resulting tables:
ARIVe_eICU_Eligibility
ARIVe_eICU_Baseline
ARIVe_eICU_Timevarying

Next step is to use the R scripts for final cohort processing and analysis.
The fullEligibility scripts were used for generating the cohort flow charts.