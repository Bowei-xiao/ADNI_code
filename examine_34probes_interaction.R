library('data.table')

# load the 34 probes identified by missoNet
d1=fread('Selected_probes_for_great_missoNet_BIC.bed',stringsAsFactors=F)
d1$chr_num=as.numeric(gsub('chr','',d1$V1))
#load phenotypes and corresponding methylation residuals
pheno=read.csv('ADNI_covariate_moreCovs_540samples_withEpiage_1905obs.csv')
pheno_sub=pheno[pheno$uniqueID==1 & !is.na(pheno$ABETA),]

pheno_sub$scaleA=scale(log(pheno_sub$ABETA)); pheno_sub$scaleM=(scale(log(pheno_sub$PTAU))+scale(log(pheno_sub$TAU)))/2;
pheno_sub$scaleD=scale(log(pheno_sub$TAU))-scale(log(pheno_sub$PTAU))

# store results in a 4d array. dim1: probe, dim2: variables estimated (3 main+2 interaction+intercept), 
## dim3: columns for summary stats (beta, se, t, p-val), dim4: 3 phenotypes (A/M/D)
res_df=array(0,dim=c(nrow(d1),6,4,3))
for (i in 1:nrow(d1)){
    # load methlation residuals
    pheno_copy=pheno_sub
    probeID=d1$V4[i]
    residual_trans=fread(paste0('Transposed_standardized_residuals_allProbes_withRowNames_noColNames_withSex_Chr',d1$chr_num[i],'.txt'),header=F,stringsAsFactors=F)
    # grep the probe
    methy=residual_trans[residual_trans$V1 == probeID,]
    methy_num=as.vector(t(methy[,-1])); pheno_copy$meth=methy_num
    pheno_name=paste0('scale',c('A','M','D'))
    for (j in 1:3){
     res=lm(paste0(pheno_name[j],'~meth*DX_num+meth*APOE4'),data=pheno_copy)
     res_df[i,,,j]=coef(summary(res))   
    }
}

# res_df is in the same order as the probe
# we are interested in which one has significant interactions
# so we search through the 4d array, find which has a p-value <=0.05, then create a new df
# new df 
interact_pval=NULL
for (i in 1:nrow(d1)){
    new_df=data.frame(probe=d1$V4[i],
                                     pA_ad=res_df[i,5,4,1],pA_apoe=res_df[i,6,4,1],
                                     pM_ad=res_df[i,5,4,2],pM_apoe=res_df[i,6,4,2],
                                     pD_ad=res_df[i,5,4,3],pD_apoe=res_df[i,6,4,3])
    interact_pval=rbind(interact_pval,new_df)
}
alpha=0.05/34/3
interact_pval[interact_pval$pA_ad<=alpha | interact_pval$pM_ad<=alpha| interact_pval$pD_ad<=alpha,]
interact_pval[interact_pval$pA_apoe<=alpha | interact_pval$pM_apoe<=alpha| interact_pval$pD_apoe<=alpha,]

# now output as 3 csv file for each phenotype. For each we output probe per row, beta, se, pval for both interaction terms on column
get_interact_df=function(probe,summary_df){
    # internal function for output_interaction. output a new df for a 2d array (selected by df[i,,,j])
    # row 5 and 6 are interaction needed. row 5 name defaults: probe_ad, row6 name defaults: probe_apoe
    return(data.frame(probe=probe,beta_probe_ad=summary_df[5,1],se_probe_ad=summary_df[5,2],pval_probe_ad=summary_df[5,4],
                        beta_probe_apoe=summary_df[6,1],se_probe_apoe=summary_df[6,2],pval_probe_apoe=summary_df[6,4]))
}
output_interaction=function(df,columnInx,probeID,output_toCSV=F,outfileName=''){
    # the df is the 4d array outputed above (res_df),columnInx \in {1,2,3} for each of the phenotype
    # probeID is a list of probes in the same order (usually taken from d1$V4 )
    intact_df=NULL
    for (i in 1:length(probeID)){
        intact_df=rbind(intact_df,get_interact_df(probeID[i],df[i,,,columnInx]))
    }
    if (output_toCSV){
        write.csv(intact_df,paste0(outfileName,'.csv'),quote=F,row.names=F)
    }
    return(intact_df)
}

intact_df_a=output_interaction(res_df,1,d1$V4,outfileName='Supplement_t5_interaction_summary_forA',output_toCSV=T)
intact_df_m=output_interaction(res_df,2,d1$V4,outfileName='Supplement_t5_interaction_summary_forM',output_toCSV=T)
intact_df_d=output_interaction(res_df,3,d1$V4,outfileName='Supplement_t5_interaction_summary_forD',output_toCSV=T)