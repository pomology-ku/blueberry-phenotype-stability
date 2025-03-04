# minimum alternative allele count of 10, minimum number of genotyped individuals of 100 
vcffilter -f "AC > 9" -f "NS > 99" VAC_155006_RAW_SNPs.vcf > VAC_155006_flt1.vcf

# only biallelic loci
cat VAC_155006_flt1.vcf | vcfbiallelic > VAC_155006_flt2.vcf

# minor allele frequency of 0.01
vcffilter -f "AF > 0.01" -f "AF < 0.99" VAC_155006_flt2.vcf > VAC_155006_flt3.vcf
