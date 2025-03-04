# import libraries
library(GWASpoly)
library(ggplot2)

# read GWASpoly files
for(name in colnames(pheno)[4:ncol(pheno)]){
  phenofile <- paste0('phenofile/phenotype_20240217_', name, '.csv')
  genofile <- paste0('genofile/geno_', name, '.csv')
  kinshipfile <- paste0('kinshipfile/Gmatrix_', name, '.csv')
  data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile, format='numeric', n.traits=1, delim=',')
  tmp_kinship <- as.matrix(read.csv(kinshipfile, header=T, row.names=1))
  data.loco <- set.K(data, K=tmp_kinship, LOCO=F, n.core=2)
  params <- set.params(n.PC=9, MAF=0.05, geno.freq=0.9)
  data.loco.scan <- GWASpoly(data=data.loco,models=c('additive', '1-dom', '2-dom'), params=params, n.core=8)

# plot result
  a = paste0("data.loco.scan@scores$", name)
  res <- cbind(data.loco.scan@map, eval(parse(text=a)))
  res['additive_p'] <- 10^(-res$additive)
  res['1-dom-alt_p'] <- 10^(-res$`1-dom-alt`)
  res['1-dom-ref_p'] <- 10^(-res$`1-dom-ref`)
  res['2-dom-alt_p'] <- 10^(-res$`2-dom-alt`)
  res['2-dom-ref_p'] <- 10^(-res$`2-dom-ref`)
  res['additive_padj'] <- p.adjust(res$additive_p, method="BH")
  res['1-dom-alt_padj'] <- p.adjust(res$`1-dom-alt_p`, method='BH')
  res['1-dom-ref_padj'] <- p.adjust(res$`1-dom-ref_p`, method='BH')
  res['2-dom-alt_padj'] <- p.adjust(res$`2-dom-alt_p`, method='BH')
  res['2-dom-ref_padj'] <- p.adjust(res$`2-dom-ref_p`, method='BH')

  file <- paste0('result_csv/result _', name, '.csv')
  write.csv(res, file, quote=F, row.names=F)

  for(type in c('additive', '1-dom-alt', '1-dom-ref', '2-dom-alt', '2-dom-ref')){
    b <- paste0('tmp_res <- data.frame(res$Marker, res$Chrom, res$Position, res$`', type, '`, res$`', type, '_p`, res$`', type, '_padj`)')
    eval(parse(text=b))
    colnames(tmp_res) <- c('marker', 'chr', 'pos', 'log10p', 'p', 'padj')
    tmp_res <- na.omit(tmp_res)
    
    g <- ggplot()
    g <- g + geom_point(data=tmp_res, mapping=aes(x=pos, y=log10p, color=chr), size=1)
    g <- g + geom_point(data=tmp_res[tmp_res$padj<0.1,], mapping=aes(x=pos, y=log10p), size=1, color='red')
    g <- g + facet_grid(. ~ chr, switch='x', scales='free')
    g <- g + scale_color_manual(values=rep(c("lightblue", "royalblue"), 6))
    g <- g + theme_classic()
    g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=20, color='black'), axis.title=element_text(size=20), strip.placement='outside', strip.background=element_blank(), strip.text=element_text(size=20), panel.spacing=unit(0, "lines"), legend.position='none')
    g <- g + scale_y_continuous(expand=c(0,0), limits=c(0, max(7, max(tmp_res$log10p)*1.1)))
    g <- g + labs(x='Chromosome', y='-log10P')
    plot(g)
    file <- paste0('manhattan _', name, '_', type, '.pdf')
    ggsave(file, g, width=10, height=4)
  }
}

# extract significant results
sigres <- data.frame()
for(name in colnames(pheno)[4:ncol(pheno)]){
  file <- paste0('result_csv/result_', name, '.csv')
  tmp_res <- read.csv(file)
  for(c in 16:ncol(tmp_res)){
    for(r in 1:nrow(tmp_res)){
      if(is.na(tmp_res[r,c]) == F){
        if(tmp_res[r,c] < 0.1){
          sigres <- rbind(sigres,
                          data.frame(trait=strsplit(name, '_')[[1]][1], year=strsplit(name, '_')[[1]][2], model=colnames(tmp_res)[c],
                                     chromosome=tmp_res[r,]$Chrom, position=tmp_res[r,]$Position, p_value=tmp_res[r,c], effect=NA))
        }
      }
    }
  }
}

# calculate marker effect
for(name in colnames(pheno)[4:ncol(pheno)]){
  tmp_name <- strsplit(name, '_')[[1]]
  if(nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2],])>0){
    phenofile <- paste0('phenofile/phenotype_', name, '.csv')
    genofile <- paste0('genofile/geno_', name, '.csv')
    kinshipfile <- paste0('kinshipfile/Gmatrix_', name,  '.csv')
    data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile,
                          format='numeric', n.traits=1, delim=',')
    tmp_kinship <- as.matrix(read.csv(kinshipfile, header=T, row.names=1))
    data.loco <- set.K(data, K=tmp_kinship, LOCO=F, n.core=2)
    params <- set.params(n.PC=9, MAF=0.05, geno.freq=0.9)
    data.loco.scan <- GWASpoly(data=data.loco,models=c('additive', '1-dom', '2-dom'), params=params, n.core=32)
    data2 <- set.threshold(data.loco.scan,method="FDR",level=0.2)
    p <- LD.plot(data2)
    p <- p + xlim(0,30) 
    file <- paste0('ld/ld_GWASpoly_', name, '_.pdf')
    ggsave(file, p, width=5, height=5)
    # estimate SNP effect
    qtl <- get.QTL(data=data2,traits=c(name),models="additive")
    if(nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='additive_padj',])>0){
      for(r in 1:nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='additive_padj',])){
        tmp_chrom <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='additive_padj',][r,]$chromosome
        tmp_pos <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='additive_padj',][r,]$position
        sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='additive_padj',][r,]$effect <- 
          qtl[qtl$Chrom==tmp_chrom & qtl$Position==tmp_pos,]$Effect
      }
    }
    qtl <- get.QTL(data=data2,traits=c(name),models="1-dom")
    if(nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.ref_padj',])>0){
      for(r in 1:nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.ref_padj',])){
        tmp_chrom <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.ref_padj',][r,]$chromosome
        tmp_pos <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.ref_padj',][r,]$position
        sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.ref_padj',][r,]$effect <- 
          qtl[qtl$Chrom==tmp_chrom & qtl$Position==tmp_pos & qtl$Model=='1-dom-ref',]$Effect
      }
    }
    if(nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.alt_padj',])>0){
      for(r in 1:nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.alt_padj',])){
        tmp_chrom <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.alt_padj',][r,]$chromosome
        tmp_pos <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.alt_padj',][r,]$position
        sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X1.dom.alt_padj',][r,]$effect <- 
          qtl[qtl$Chrom==tmp_chrom & qtl$Position==tmp_pos & qtl$Model=='1-dom-alt',]$Effect
      }
    }
    qtl <- get.QTL(data=data2,traits=c(name),models="2-dom")
    if(nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.ref_padj',])>0){
      for(r in 1:nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.ref_padj',])){
        tmp_chrom <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.ref_padj',][r,]$chromosome
        tmp_pos <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.ref_padj',][r,]$position
        sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.ref_padj',][r,]$effect <- 
          qtl[qtl$Chrom==tmp_chrom & qtl$Position==tmp_pos & qtl$Model=='2-dom-ref',]$Effect
      }
    }
    if(nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.alt_padj',])>0){
      for(r in 1:nrow(sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.alt_padj',])){
        tmp_chrom <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.alt_padj',][r,]$chromosome
        tmp_pos <- sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.alt_padj',][r,]$position
        sigres[sigres$trait==tmp_name[1] & sigres$year==tmp_name[2] & sigres$model=='X2.dom.alt_padj',][r,]$effect <- 
          qtl[qtl$Chrom==tmp_chrom & qtl$Position==tmp_pos & qtl$Model=='2-dom-alt',]$Effect
      }
    }
  }
}
