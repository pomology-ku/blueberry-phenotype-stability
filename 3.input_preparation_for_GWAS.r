# import libraries
library(AGHmatrix)

# read genotype and phenotype file
geno <- read.table('updog_result_geno.txt', header=T, sep='\t')
pheno <- read.csv('GWAS_phenotype.csv', header=T)

# individuals of minimum missing genotype of 0.7
miss_ind <- apply(geno, 2, function(x){sum(is.na(x))/length(x)})
geno_flt <- geno[,miss_ind<0.7]

# SNPs of minimum missing genotype of 0.2
miss_snp <- apply(geno_flt, 1, function(x){sum(is.na(x))/length(x)})
geno_flt <- geno_flt[miss_snp < 0.2,]

# minor allele frequency of 0.05
allele_freq <- data.frame(marker=rownames(geno_flt), ref=rep(0, nrow(geno_flt)), alt=rep(0, nrow(geno_flt)))
for(r in 1:nrow(geno_flt)){
  for(c in 1:ncol(geno_flt)){
    if(is.na(geno_flt[r,c])==FALSE){
      allele_freq[r,]$ref <- allele_freq[r,]$ref + geno_flt[r,c]
      allele_freq[r,]$alt <- allele_freq[r,]$alt + (4 - geno_flt[r,c])
    }
  }
}
allele_freq['freq'] <- (allele_freq$ref)/(allele_freq$ref + allele_freq$alt)
geno_flt <- geno_flt[allele_freq$freq>0.05 & allele_freq$freq<0.95,]
kinship <- geno_flt

# create a map file 
tmp_marker <- rownames(geno_flt)
tmp_chr <- c()
tmp_pos <- c()
for(tmp in tmp_marker){
  tmp1 <- strsplit(tmp, 'p0.')[[1]]
  tmp2 <- strsplit(tmp1[2], ':')[[1]]
  tmp3 <- strsplit(tmp2[1], 'Chr')[[1]]
  tmp4 <- strsplit(tmp2[2], '_')[[1]]
  if(length(tmp3)==2){
    tmp_chr <- c(tmp_chr, tmp3[2])
  }else{
    tmp_chr <- c(tmp_chr, NA)
  }
  tmp_pos <- c(tmp_pos, tmp4[1])
}
geno_flt <- cbind(data.frame(marker=tmp_marker, chr=tmp_chr, pos=tmp_pos), geno_flt)
geno_flt <- geno_flt[!is.na(geno_flt$chr),]

# prepare phenotype file
for(c in 4:ncol(pheno)){
  tmp_pheno <- pheno[,c(1, c)]
  tmp_pheno <- na.omit(tmp_pheno)
  tmp_pheno <- tmp_pheno[tmp_pheno$RAPiD_Genomics_Sample_Code%in%colnames(geno_flt),]
  tmp_geno <- geno_flt[,c(rep(TRUE, 3), colnames(geno_flt)[4:ncol(geno_flt)]%in%tmp_pheno$RAPiD_Genomics_Sample_Code)]
  tmp_ave <- apply(tmp_geno[,4:ncol(tmp_geno)], 1, function(x){mean(x, na.rm=T)})
  tmp_ave <- round(tmp_ave)
  for(r in 1:nrow(tmp_geno)){
    tmp_geno[r,][is.na(tmp_geno[r,])] <- tmp_ave[r]
  }
# preparation of kinship matrix
  tmp_kinship <- as.matrix(t(tmp_geno[,4:(ncol(tmp_geno))]))
  tmp_kinship[is.na(tmp_kinship)] <- -9
  tmp_K <- Gmatrix(tmp_kinship, method="VanRaden", ploidy=4)
  file <- paste0('phenofile/phenotype_', colnames(pheno)[c], '.csv')
  write.csv(tmp_pheno, file, quote=F, row.names=F)
  file <- paste0('genofile/geno_', colnames(pheno)[c], '.csv')
  write.csv(tmp_geno, file, quote=F, row.names=F)
  file <- paste0('kinshipfile/Gmatrix_', colnames(pheno)[c], '.csv')
  write.csv(tmp_K, file, quote=F, row.names=T)
}
