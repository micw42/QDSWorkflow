musc_ann = data.frame(sample_id=colnames(musc_raw)[-1],
                      Condition=rep(c("fixed", "perfused",
                                      "unperfused", "activated"), each=2))
setwd("QDSWorkflow/")
usethis::use_data(musc_ann)