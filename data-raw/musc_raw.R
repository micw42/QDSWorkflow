musc_raw= data.table::fread("featureCounts_musc_regen.txt")
setwd("QDSWorkflow/")
usethis::use_data(musc_raw)