library('data.table')
#check top probes in CLSA data

# load CLSA data first
# covariates loading
load('./CLSA/CLSA.covariate_modified.RData')
# covariates to grab: SEX_ASK_COM, AGE_NMBR_COM, CCC_ALZH_COM (phneotype?), 
# load cell decomposition
load('./CLSA/CLSA.PC_Estimated_Cell_Proportion.RData')

# load methylation
load('./CLSA/CLSA.Beta_funNorm.RData')


# They are already ordered so we just select the probes from emthylation then bind with phenotypes/covariates
misso_probe=fread('Selected_probes_for_great_missoNet_BIC.bed',stringsAsFactors=F)
names(misso_probe)=c('chr','start','end','probe')
# subset methylation, then transpose
methy_sub=t(datMeth[rownames(datMeth) %in% misso_probe$probe,])

# now map to covariate file, match by row name from DATA & methy_sub
# simple solution: check if the rownames are the same order, if so, can directly cbind
sum(rownames(methy_sub)!=rownames(DATA))
# sum is 0, so safe to cbind. Let's just save the useful columns from DATA
replicate_df=cbind(DATA[,c('SEX_ASK_COM','AGE_NMBR_COM','CCC_ALZH_COM','COG_REYI_SCORE_COM','COG_MAT_SCORE_COM','DNAmAge_COM')],methy_sub)
# adding PCA here
PCs=pcaC$x
replicate_df=cbind(replicate_df,PCs)

# calculate EAA
replicate_df$eaa=replicate_df$DNAmAge_COM-replicate_df$AGE_NMBR_COM

# adding year of education. Note, loading this will overwrite the original variable DATA. Make sure to reload if need to use it again
load('./CLSA/CLSA.covariate_modified.SocioEconomic.RData')
# one column called ED_UDR11_COM seems usable. Alternatively, use ED_HSGR_COM to binarize
# are they the same order?
sum(rownames(methy_sub)!=DATA$ADM_EPIGEN2_COM)
replicate_df$yoe=as.factor(DATA$ED_HIGH_COM) 

#alternative pheotypes: COG_REYI_SCORE_COM and CCC_ALZH_COM

pheno='CCC_ALZH_COM'
#pheno='COG_REYI_SCORE_COM'
summary_res=NULL;

for (probe in misso_probe$probe){
    res1=glm(as.formula(paste0(pheno,'~',probe,'+as.factor(SEX_ASK_COM)+eaa+AGE_NMBR_COM+PC1+PC2+PC3+yoe')),data=replicate_df[replicate_df$CCC_ALZH_COM!=8,],family="quasipoisson")
    summary_res=rbind(summary_res,data.frame(probe=probe,beta=coef(summary(res1))[2,1],se=coef(summary(res1))[2,2],pval=coef(summary(res1))[2,4]))
}

write.csv(summary_res,'replication_CLSA_using_AD_phenotype.csv',quote=F,row.names=F)
#write.csv(summary_res,'replication_CLSA_using_REY_phenotype.csv',quote=F,row.names=F)