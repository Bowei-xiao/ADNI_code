library('data.table')

ewas_summary=NULL
includeSex=T
if (includeSex){
    name_suffix='withSex'
} else {
    name_suffix='noSex'
}

for (chr_num in 1:22){
    ewas_summary=rbind(ewas_summary,fread(paste0('Univariate_ewas_3phenos_allProbes_',name_suffix,'_chr',chr_num,'.csv'),stringsAsFactors=F))
}

ewas_summary[which.min(ewas_summary$p_a_probe),]


# hand calculate lambda
calculate_lambda=function(ewas_summary,pheno=c('A','M','D')){
    return(median(qchisq(as.vector(unlist(ewas_summary[,paste0('p_',tolower(pheno),'_probe'),with=F])), df=1, lower.tail=FALSE)) / qchisq(0.5, 1))
    
}
lambda_a=calculate_lambda(ewas_summary,'A')
lambda_t=calculate_lambda(ewas_summary,'M')
lambda_n=calculate_lambda(ewas_summary,'D')

# bacon correction
# forgot to change the variable names: a ->A, t->M, n->D
library('bacon')
bc_a <- bacon(teststatistics = NULL,effectsizes =  ewas_summary$beta_a_probe,standarderrors = ewas_summary$se_a_probe,na.exclude = TRUE,verbose = F)
bc_t <- bacon(teststatistics = NULL,effectsizes =  ewas_summary$beta_t_probe,standarderrors = ewas_summary$se_t_probe,na.exclude = TRUE,verbose = F)
bc_n <- bacon(teststatistics = NULL,effectsizes =  ewas_summary$beta_n_probe,standarderrors = ewas_summary$se_n_probe,na.exclude = TRUE,verbose = F)

ewas_a_ca_corrected=data.frame(probe=ewas_summary$probe,chr=ewas_summary$chr,location=ewas_summary$location,
                              beta=bacon::es(bc_a),se=bacon::se(bc_a),pval=pval(bc_a),fdr=p.adjust(pval(bc_a), method = "fdr"),stringsAsFactors = FALSE)

ewas_t_ca_corrected=data.frame(probe=ewas_summary$probe,chr=ewas_summary$chr,location=ewas_summary$location,
                              beta=bacon::es(bc_t),se=bacon::se(bc_t),pval=pval(bc_t),fdr=p.adjust(pval(bc_t), method = "fdr"),stringsAsFactors = FALSE)

ewas_n_ca_corrected=data.frame(probe=ewas_summary$probe,chr=ewas_summary$chr,location=ewas_summary$location,
                              beta=bacon::es(bc_n),se=bacon::se(bc_n),pval=pval(bc_n),fdr=p.adjust(pval(bc_n), method = "fdr"),stringsAsFactors = FALSE)


lambda_a_corrected=median(qchisq(ewas_a_ca_corrected$pval, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
lambda_t_corrected=median(qchisq(ewas_t_ca_corrected$pval, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
lambda_n_corrected=median(qchisq(ewas_n_ca_corrected$pval, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)


# manual qq plot for p-values before and after bacon correction
plot_qq_2pvals=function(pval1,pval2,outName){
    # assume pval1 and pval2 have the same length
    # pval1 is the main plot (colored in black), pval2 is the compared dots (colored in grey)
    o1 = -log10(sort(pval1,decreasing=FALSE))
    o2 = -log10(sort(pval2,decreasing=FALSE))
    e = -log10( ppoints(length(pval1) ))

    png(outName)
    plot(e, o1, pch=20, col='black', cex=1.5,
          xlab=expression(Expected~~-log[10](italic(p))), 
          ylab=expression(Observed~~-log[10](italic(p))))
    points(e,o2,pch=20,col='grey',cex=1.5)
    abline(0,1,col="red")
    dev.off()
}

plot_qq_2pvals(pval2=ewas_summary$p_a_probe,pval1=ewas_a_ca_corrected$pval,outName=paste0('qqPlot_logA_',name_suffix,'_before_after_bacon_correction.png'))
plot_qq_2pvals(pval2=ewas_summary$p_t_probe,pval1=ewas_t_ca_corrected$pval,outName=paste0('qqPlot_logM_',name_suffix,'_before_after_bacon_correction.png'))
plot_qq_2pvals(pval2=ewas_summary$p_n_probe,pval1=ewas_n_ca_corrected$pval,outName=paste0('qqPlot_logD_',name_suffix,'_before_after_bacon_correction.png'))


# plots
png(paste0('ManhattanPlot_logA','_',name_suffix,'_bacon_correction_threshold_1e-5.png'))
manhattan(ewas_a_ca_corrected,chr='chr',bp='location',p='pval',snp='probe',genomewideline=-log10(3.6e-8/3),suggestiveline=-log10(1e-5/3),ylim=c(0,11))
dev.off()
png(paste0('qqPlot_logA_',name_suffix,'_bacon_correction.png'))
qq(ewas_a_ca_corrected$pval,cex=1.5)
dev.off()
png(paste0('hist_pval_logA_',name_suffix,'_bacon_correction.png'))
hist(ewas_a_ca_corrected$pval)
dev.off()

png(paste0('ManhattanPlot_logM','_',name_suffix,'_bacon_correction_threshold_1e-5.png'))
manhattan(ewas_t_ca_corrected,chr='chr',bp='location',p='pval',snp='probe',genomewideline=-log10(3.6e-8/3),suggestiveline=-log10(1e-5/3),ylim=c(0,11))
dev.off()
png(paste0('qqPlot_logM_',name_suffix,'_bacon_correction.png'))
qq(ewas_t_ca_corrected$pval,cex=1.5)
dev.off()
png(paste0('hist_pval_logM_',name_suffix,'_bacon_correction.png'))
hist(ewas_t_ca_corrected$pval)
dev.off()



png(paste0('ManhattanPlot_logD','_',name_suffix,'_bacon_correction_threshold_1e-5.png'))
manhattan(ewas_n_ca_corrected,chr='chr',bp='location',p='pval',snp='probe',genomewideline=-log10(3.6e-8/3),suggestiveline=-log10(1e-5/3),ylim=c(0,11))
dev.off()
png(paste0('qqPlot_logD_',name_suffix,'_bacon_correction.png'))
qq(ewas_n_ca_corrected$pval,cex=1.5)
dev.off()
png(paste0('hist_pval_logD_',name_suffix,'_bacon_correction.png'))
hist(ewas_n_ca_corrected$pval)
dev.off()
