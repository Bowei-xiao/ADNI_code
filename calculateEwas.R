library('data.table')
library('car')
library('foreach')
library('doMC')
registerDoMC(22)

args=(commandArgs(TRUE))

chr_num=as.numeric(args[1])

pheno=read.csv('ADNI_covariate_moreCovs_540samples_withEpiage_1905obs.csv',stringsAsFactors=F)
probe_loc=fread('probe_locations.csv')

# load residuals of all individuals (note they have the same order as pheno)
# skip the first row, make the second row as header and then convert all entries into numbers
includeSex=T
if (includeSex){
    name_suffix='withSex'
} else {
    name_suffix='noSex'
}


residual_trans=fread(paste0('Transposed_standardized_residuals_allProbes_withRowNames_noColNames_',name_suffix,'_Chr',chr_num,'.txt'),header=F,stringsAsFactors=F)


runMultivar=function(i, pheno,residual_file,probe_locationFile){
    probe_meth=as.vector(t(residual_file[i,2:ncol(residual_file)]))
    probe_name=as.character(residual_file[i,1])
    pheno$methy=probe_meth
    
    # calculate two new phenotypes
    pheno$meanTN=(scale(log(pheno$PTAU))+scale(log(pheno$TAU)))/2
    pheno$diffTN=scale(log(pheno$TAU))-scale(log(pheno$PTAU))
    # location use start position
    probe_location=probe_locationFile$Start_hg38[probe_locationFile$Name==probe_name]


    # run A/M/D separately, 
    res11=lm(scale(log(ABETA))~methy+DX_num+APOE4, data=pheno)
    res12=lm(meanTN~methy+DX_num+APOE4, data=pheno)
    res13=lm(diffTN~methy+DX_num+APOE4, data=pheno)


    df=data.frame(probe=probe_name,chr=chr_num,location=probe_location,
                    beta_a_probe=coef(summary(res11))['methy',1],se_a_probe=coef(summary(res11))['methy',2],p_a_probe=coef(summary(res11))['methy',4],
                    beta_a_dx=coef(summary(res11))['DX_num',1],se_a_dx=coef(summary(res11))['DX_num',2],p_a_dx=coef(summary(res11))['DX_num',4],
                    beta_a_ap=coef(summary(res11))['APOE4',1],se_a_ap=coef(summary(res11))['APOE4',2],p_a_ap=coef(summary(res11))['APOE4',4],
                    beta_m_probe=coef(summary(res12))['methy',1],se_t_probe=coef(summary(res12))['methy',2],p_m_probe=coef(summary(res12))['methy',4],
                    beta_m_dx=coef(summary(res12))['DX_num',1],se_t_dx=coef(summary(res12))['DX_num',2],p_m_dx=coef(summary(res12))['DX_num',4],
                    beta_m_ap=coef(summary(res12))['APOE4',1],se_t_ap=coef(summary(res12))['APOE4',2],p_m_ap=coef(summary(res12))['APOE4',4],
                    beta_d_probe=coef(summary(res13))['methy',1],se_n_probe=coef(summary(res13))['methy',2],p_d_probe=coef(summary(res13))['methy',4],
                    beta_d_dx=coef(summary(res13))['DX_num',1],se_n_dx=coef(summary(res13))['DX_num',2],p_d_dx=coef(summary(res13))['DX_num',4],
                    beta_d_ap=coef(summary(res13))['APOE4',1],se_n_ap=coef(summary(res13))['APOE4',2],p_d_ap=coef(summary(res13))['APOE4',4], stringsAsFactors=F)
    
    return(df)
}

# run in parallel

res_list=foreach(i=1:nrow(residual_trans), .combine=rbind) %dopar% {
   res=runMultivar(i, pheno=pheno,residual_file=residual_trans,probe_locationFile=probe_loc)
}
write.csv(res_list,paste0('Univariate_ewas_3phenos_allProbes_',name_suffix,'_chr',chr_num,'.csv'),quote = F,row.names = F)
