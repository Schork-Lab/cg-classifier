library("boot")
library("plyr")
library("dplyr")
library("pROC")
library("doMC")

registerDoMC(10)

fs1 = read.delim("/gpfs/home/ekramer/Projects/cg-classifier/data/NA12877trainedClassifier_predictions_for_NA12878_indels_LR_RF_GBC_unscaled_allfeatures_noCGAXR.tsv")
fs2 = read.delim("/gpfs/home/ekramer/Projects/cg-classifier/data/NA12877trainedClassifier_predictions_for_NA12878_indels_LR_RF_GBC_unscaled_swgrFeatures.tsv")

tmp1 = fs1 %>%
  select(CHROM, 
         POS, 
         REF, 
         LR.FS1=LogReg_C1_l1_probeTP, 
         RF.FS1=RF_500est_gini_split4_leaf2_probeTP,
         GBC.FS1=GBC_150est_learn01_maxdepth5_probeTP,
         TRUTH=truth)

tmp2 = fs2 %>%
  select(CHROM, 
         POS, 
         REF, 
         LR.FS2=LogReg_C1_l1_probeTP, 
         RF.FS2=RF_500est_gini_split4_leaf2_probeTP,
         GBC.FS2=GBC_150est_learn01_maxdepth5_probeTP)

d = inner_join(tmp1, tmp2)

### calculate ROC curves

rocs = llply(d[c("LR.FS1", "RF.FS1", "GBC.FS1", "LR.FS2", "RF.FS2", "GBC.FS2")], 
             function(x) roc(d$TRUTH, as.numeric(as.character(x))),
             .parallel=T)
save(rocs, file="../data/rocs.Rdata")

### calculate bootstrapped CIs for AUCs

ci.aucs = llply(rocs, ci.auc, boot.n=100)
save(ci.aucs, file="../data/ci.aucs.Rdata")

### run pairwise tests 

pairs = combn(names(rocs), 2)
pairs = as.data.frame(t(pairs))
colnames(pairs) = c("ROC1", "ROC2")

roc.tests = mlply(pairs, function(ROC1, ROC2) roc.test(rocs[[ROC1]], rocs[[ROC2]]))
save(roc.tests, "../data/roc.test.Rdata")


