rat_human_hm = data.table::fread("Downloads/rat_human_homology (1).txt")
setwd("QDSWorkflow/")
usethis::use_data(rat_human_hm)