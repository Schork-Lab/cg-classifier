source("vcf_parser.R")

## Parse VCF

training.dir = "/gpfs/home/erscott/Datasets_raw/PlatinumGenomes/NA12877/NA12877_phasing/"

tp = readVcf(paste(training.dir, "tp.vcf.gz", sep=""), genome="hg19")
tp.features = vcfToDf(tp)

fp = readVcf(paste(training.dir, "fp.vcf.gz", sep=""), genome="hg19")
fp.features = vcfToDf(fp)

save(tp.features, fp.features, file="features.Rdata")

## Create a small subset for testing scripts

tp.features$Y = "TP"
fp.features$Y = "FP"

d = rbind(tp.features[sample(nrow(tp.features), 5000), ],
          fp.features[sample(nrow(fp.features), 5000), ])

d.testing = rbind(tp.features[sample(nrow(tp.features), 5000), ],
                  fp.features[sample(nrow(fp.features), 5000), ])

d$Y = factor(d$Y)
d.testing$Y = factor(d$Y)

d = select(d, -VarID, -IID, -FT)
d.testing = select(d.testing, -VarID, -IID, -FT)

save(d, d.testing, file="sample_data.Rdata")
