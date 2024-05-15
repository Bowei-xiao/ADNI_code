library('data.table')
library('cglasso')
library('dplyr')
library('missoNet')

args=(commandArgs(TRUE))
#clusterID=as.numeric(args[1])

chr_num=as.numeric(args[1])
nCluster=30

includeSex=T
if (includeSex){
    name_suffix='withSex'
} else {
    name_suffix='noSex'
}

missoNet_to_df=function(coefB,probe_name){
    # input missonet output beta matrix, and a name vector, return a dataframe
    if (is.null(coefB)){
        coef_df=data.frame(matrix(ncol = 4, nrow = 0))
        names(coef_df)=c('probe','A','T','N')
    } else {
     coef_df=data.frame(probe=c('APOE4','DX_num',probe_name),A=coefB[,1],T=coefB[,2],N=coefB[,3])
    }
    return(coef_df)
}

# load methylation residuals
residual_trans=fread(paste0('Transposed_standardized_residuals_top25percent_withRowNames_noColNames_',name_suffix,'_Chr',chr_num,'.txt'),header=F,stringsAsFactors=F)


# read phenotypes
pheno=read.csv('ADNI_covariate_moreCovs_540samples_withEpiage_1905obs.csv')
pheno_sub=pheno[pheno$uniqueID==1 & !is.na(pheno$ABETA),]
# create the dataset format for cglasso

Y=as.matrix(pheno_sub[,c('ABETA','PTAU','TAU')])
# transformation of Y
#normalize A/T/N separately
Y=as.matrix(cbind(scale(log(pheno_sub$ABETA)),scale(log(pheno_sub$PTAU)),scale(log(pheno_sub$TAU))))

probe_loc=fread('probe_locations.csv')
inx1=which(probe_loc$chr_num == chr_num)
probe_loc_oneChr=probe_loc[inx1,]

# This step is not necessary, the input is sorted. But keep it here just in case.
residual_trans_sorted=residual_trans %>% arrange(factor(V1, levels = probe_loc_oneChr$Name))

nStart=1; nEnd=nrow(residual_trans_sorted)
clusterID=1;
nlam=200

while (nStart<nEnd){
clusterEnd=min(nStart+nCluster-1,nEnd)
residual_trans_inOneCluster=residual_trans_sorted[nStart:clusterEnd,]

print(paste0('This is cluster ',clusterID))
print(paste0('Number of probes used for cgLasso is ', dim(residual_trans_inOneCluster)))


X_oneCluster=t(residual_trans_inOneCluster[,2:ncol(residual_trans_inOneCluster)]) 
X_oneCluster=as.data.frame(X_oneCluster)
names(X_oneCluster)=residual_trans_inOneCluster$V1
# adding two covariates, setting their weights to 0 (not penalized at all)
X_withCov=cbind(pheno_sub$APOE4,pheno_sub$DX_num,X_oneCluster); names(X_withCov)[1:2]=c('APOE4','DX_num')
weightB=matrix(1,nrow=ncol(X_withCov),ncol=ncol(Y))
weightB[1,]=0; weightB[2,]=0
#====== run cglasso
# construct Z
Z=datacggm(Y=Y,X=X_withCov)

#Do not penalize theta. Just use rho=0 and become 1-d search.
print(Sys.time())
res1=cglasso(.~.,data=Z,nlambda=nlam,rho=0,lambda.min.ratio = 1e-5,weights.B=weightB,maxit.em =1e+5)
print(Sys.time())

# best model by BIC
out_res1=summary(res1,GoF=BIC,print.all=F)
# get betas
coef_B1=coef(res1,type='B',lambda.id=out_res1$lambda.id,rho.id=out_res1$rho.id)

write.csv(coef_B1,paste0('cgLasso_chr',chr_num,'_cluster',clusterID,'_clusterSize_',nCluster,'_',name_suffix,'_useBIC_betaOutputs.csv'),quote=F)

#==== run missonet with cv min


res_misso_cv = cv.missoNet(X = X_withCov, Y = Y, kfold=10,
                        lamB.min.ratio = 1e-5, 
                        n.lamB = nlam,
                        lambda.Theta = 0, 
                        lamB.penalty.factor =weightB, lamTh.penalty.factor = matrix(0,3,3),
                        shuffle=FALSE,penalize.diagonal = FALSE,
                        fast = FALSE, verbose = 0)

# extract needed betas
coefB_misso_cvmin = res_misso_cv$est.min$Beta

write.csv(missoNet_to_df(coefB_misso_cvmin,names(X_oneCluster)),paste0('missoNet_chr',chr_num,'_cluster',clusterID,'_clusterSize_',nCluster,'_',name_suffix,'_useCVmin_betaOutputs.csv'),quote=F,row.names=F)

