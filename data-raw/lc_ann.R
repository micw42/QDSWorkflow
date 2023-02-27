lc_ann = lc_ann %>% 
  mutate(Development_Stage=case_when(Sample_Origin=="nLung" ~"normal lung",
                                     Sample_Origin=="tLung" ~"lung tumor (early)",
                                     Sample_Origin=="tL/B" ~"lung tumor (late)",
                                     Sample_Origin=="mLN" ~"lymph node (metastasis)",
                                     Sample_Origin=="PE" ~"pleural effusion (metastasis)",
                                     Sample_Origin=="mBrain" ~"brain (metastasis)")) %>%
  mutate(Development_Stage=factor(Development_Stage, ordered = T, levels=c("normal lung", "lung tumor (early)",
                                                                           "lung tumor (late)", "lymph node (metastasis)",
                                                                           "pleural effusion (metastasis)", "brain (metastasis)")))
setwd("QDSWorkflow/")
usethis::use_data(lc_ann, overwrite = T)