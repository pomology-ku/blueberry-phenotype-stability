# import libraries
library(updog)
library(VariantAnnotation)
library(future)

# read filtered vcf file
f1 = readVcf("VAC_155006_flt3.vcf")

# make input for updog
refmat = geno(f1)$RO
numind = dim(refmat)[2]
aotemp = geno(f1)$AO
AO = matrix(sapply(aotemp[,],"[[",1),ncol = numind)
sizemat = refmat+AO

# run updog
ploidy = 4
mout <- multidog(refmat = refmat, sizemat = sizemat, ploidy = ploidy, model = "norm", nc = 2)
mout2<-filter_snp(mout,prop_mis<0.05)
res1<-format_multidog(mout2, varname = "geno")
res2<-format_multidog(mout2, varname = "postmean")
write.table(res1,"updog_result_geno.txt",quote=F,sep="\t")
write.table(res2,"updog_result_postmean.txt",quote=F,sep="\t")
