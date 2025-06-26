library('data.table')
library('cglasso')
library('dplyr')
library('glmnet')
#library('Rcpp')
#library('foreach')
#library('doMC')
library('missoNet')
#registerDoMC(20)
args=(commandArgs(TRUE))
chr_num=22; genePerChunk=100; cvConsensus=3
chunk=as.numeric(args[1])
includeSex=T
if (includeSex){
    name_suffix='withSex'
} else {
    name_suffix='noSex'
}
set.seed(1020)
missoNet_to_df=function(coefB,probe_name){
    # input missonet output beta matrix, and a name vector, return a dataframe
    if (is.null(coefB)){
        coef_df=data.frame(matrix(ncol = 4, nrow = 0))
        names(coef_df)=c('probe','A','T','N')
    } else {
     coef_df=data.frame(probe=c('DX_num','APOE4',probe_name),A=coefB[,1],T=coefB[,2],N=coefB[,3])
    }
    return(coef_df)
}

residual_trans=fread(paste0('Transposed_standardized_residuals_allProbes_withRowNames_noColNames_',name_suffix,'_Chr',chr_num,'.txt'),header=F,stringsAsFactors=F)

summarize_5CVrun=function(array1,consensus){
    # load a 3D array (created in the main loop that stored the betas of 5 CV runs)
    # output a 2d df of the same size as dim(array1)[c(1,2)] and same order
    # for each entry, if there are >=consenus non-zero among all runs then it will become the average values (excluing 0)
    # othersie, this entry will become 0.
    df=matrix(0,dim(array1)[1],dim(array1)[2])
        for (i in 1:dim(df)[1]){
            for (j in 1:dim(df)[2]){
                est=array1[i,j,]
                if (sum(est!=0)>=consensus){
                    # calculate average (removing 0)
                    df[i,j]=mean(est[est!=0])
                }
            }
        }
    return(df)
}

# read phenotypes
pheno=read.csv('ADNI_covariate_moreCovs_540samples_withEpiage_1905obs.csv')
pheno_sub=pheno[pheno$uniqueID==1 & !is.na(pheno$ABETA),]
# create the dataset format for cglasso

# transformation of Y
#normalize A/T/N separately
Y=as.matrix(cbind(scale(log(pheno_sub$ABETA)), (scale(log(pheno_sub$PTAU))+scale(log(pheno_sub$TAU)))/2, scale(log(pheno_sub$TAU))-scale(log(pheno_sub$PTAU))))

probe_loc=fread(paste0('probe_chr',chr_num,'_geneLabel_100k.csv'))
geneList=names(probe_loc)[-1]

# This step is not necessary, the input is sorted. But keep here just in case.
residual_trans_sorted=residual_trans %>% arrange(factor(V1, levels = probe_loc$probe))
nlam=200

# subset the gene list for the chunk
geneList_sub=geneList[(genePerChunk*(chunk-1)+1):min(chunk*genePerChunk,length(geneList))]
for (geneID in geneList_sub){
# here we will do something different. We will cluster probes by genes
# however, based on our little test, it seems safer to cap cluster size at 350.

    probeList=probe_loc$probe[probe_loc[,..geneID] == 1]
    residual_trans_inOneCluster=residual_trans_sorted[residual_trans_sorted$V1 %in% probeList,]

    print(paste0('This is cluster for gene ',geneID))
    print(paste0('Number of probes used for cgLasso is ', dim(residual_trans_inOneCluster)[1]))
    print('Start generating X, Y, weights')
    print(Sys.time())
    X_oneCluster=t(residual_trans_inOneCluster[,2:ncol(residual_trans_inOneCluster)]) 
    X_oneCluster=as.data.frame(X_oneCluster)
    names(X_oneCluster)=residual_trans_inOneCluster$V1
    X_withCov=cbind(pheno_sub$DX_num,pheno_sub$APOE4,X_oneCluster); names(X_withCov)[1]='DX_num'; names(X_withCov)[2]='APOE4'
    # construct Z
    Z=datacggm(Y=Y,X=X_withCov)

    #==== run missonet
    # use BIC, cv min and cv 1se
    print(paste0('missonet-BIC for gene',geneID,' starts at',Sys.time()))
    # change weight back to 0 for missonet
    weightB=matrix(1,nrow=ncol(X_withCov),ncol=ncol(Y))
    # no penalization to DX and apoe4 (1e-9 to avoid error)
    weightB[c(1,2),]=0
    res_misso_bic = missoNet(X = X_withCov, Y = Y, GoF = "eBIC",
                            lamB.min.ratio = 1e-5, 
                            n.lamB = nlam,
                            lambda.Theta = 0, 
                            lamB.penalty.factor = weightB,
                            penalize.diagonal = FALSE,
                            fast = FALSE, verbose = 0)
    print(paste0('missonet-CV for gene',geneID,' starts at',Sys.time()))
    # output lambda series for cglasso
    lambda.Beta = res_misso_bic$lambda.Beta
    # extract needed betas
    coefB_misso_bic = res_misso_bic$est.min$Beta


    # run 5 CVs and for probes >=3 times with non-zero beta on either Y, save the avg of these non-zero betas at the corresponding phenotype
    # If less than 3 times, just assign 0 
    cv1se_5run_array=array(0,dim=c(ncol(X_withCov),ncol(Y),5))
    cvmin_5run_array=array(0,dim=c(ncol(X_withCov),ncol(Y),5))
    for (i in 1:5){
        res_misso_cv = cv.missoNet(X = X_withCov, Y = Y, kfold=10,
                                lamB.min.ratio = 1e-5, 
                                n.lamB = nlam,
                                lambda.Theta = 0, 
                                lamB.penalty.factor =weightB,
                                shuffle=TRUE,penalize.diagonal = FALSE,
                                fast = FALSE, verbose = 0)
        print(paste0('missonet-CV for gene',geneID,' run',i,' ends at',Sys.time()))

        if (!is.null(res_misso_cv$est.min)){
            cvmin_5run_array[,,i] = res_misso_cv$est.min$Beta
        }
        if (!is.null(res_misso_cv$est.1seB)){
            cv1se_5run_array[,,i] = res_misso_cv$est.1seB$Beta
        }
    }
    coefB_misso_cvmin=summarize_5CVrun(cvmin_5run_array,consensus=cvConsensus)
    coefB_misso_cv1se=summarize_5CVrun(cv1se_5run_array,consensus=cvConsensus)

   
    write.csv(missoNet_to_df(coefB_misso_bic,names(X_oneCluster)),paste0('Output_missoNet_chr_',chr_num,'_gene_',geneID,'_useBIC_betaOutputs.csv'),quote=F,row.names=F)
    write.csv(missoNet_to_df(coefB_misso_cvmin,names(X_oneCluster)),paste0('Output_missoNet_chr_',chr_num,'_gene_',geneID,'_useCVmin_betaOutputs.csv'),quote=F,row.names=F)
    write.csv(missoNet_to_df(coefB_misso_cv1se,names(X_oneCluster)),paste0('Output_missoNet_chr_',chr_num,'_gene_',geneID,'_useCV1se_betaOutputs.csv'),quote=F,row.names=F)
    
    #====== run cglasso
    weightB=matrix(1,nrow=ncol(X_withCov),ncol=ncol(Y))
    # no penalization to DX and apoe4 (1e-9 to avoid error)
    weightB[c(1,2),]=1e-9
    #Do not penalize theta. Just use rho=0 and become 1-d search.
    print(paste0('cglasso for gene',geneID,' starts at',Sys.time()))
    res1=cglasso(.~.,data=Z,diagonal = FALSE, rho=0, lambda = lambda.Beta, weights.B=weightB, maxit.em=1e+5)
    print(paste0('cglasso for gene',geneID,' ends at',Sys.time()))

    # best model by BIC
    out_res1=summary(res1,GoF=BIC,print.all=F)
    # get betas
    coef_B1=coef(res1,type='B',lambda.id=out_res1$lambda.id,rho.id=out_res1$rho.id)
    write.csv(coef_B1,paste0('Output_cgLasso_chr_',chr_num,'_gene_',geneID,'_useBIC_betaOutputs.csv'),quote=F)

    
    # ====== run Lasso
    print(paste0('lasso for gene',geneID,' starts at',Sys.time()))
    X_mat=as.matrix(X_withCov,nrow=nrow(X_withCov),ncol=ncol(X_withCov))

    # let's also do 5 run lasso. Note lasso output 1 extra column: intercept.
    cv1se_lasso_5run_array=array(0,dim=c((ncol(X_withCov)+1),ncol(Y),5))
    cvmin_lasso_5run_array=array(0,dim=c((ncol(X_withCov)+1),ncol(Y),5))
    for (i in 1:5){
        coefB_lasso_cvmin=NULL; coefB_lasso_cv1se=NULL
        for (j in 1:3){
            cv_fit=cv.glmnet(x=X_mat,y=Y[,j],lambda.min.ratio=1e-5,nlambda=nlam,alpha=1,penalty.factor=c(1e-9,1e-9,rep(1,ncol(X_mat)-2)))
            best_fitmin=coef(cv_fit,s='lambda.min',exact=T);
            best_fit1se=coef(cv_fit,s='lambda.1se',exact=T);

            coefB_lasso_cvmin=cbind(coefB_lasso_cvmin,as.matrix(best_fitmin))
            coefB_lasso_cv1se=cbind(coefB_lasso_cv1se,as.matrix(best_fit1se))
        }
        cv1se_lasso_5run_array[,,i]=coefB_lasso_cv1se
        cvmin_lasso_5run_array[,,i]=coefB_lasso_cvmin
    }
    coefB_lasso_cvmin_5run=summarize_5CVrun(cvmin_lasso_5run_array,consensus=cvConsensus)
    coefB_lasso_cv1se_5run=summarize_5CVrun(cv1se_lasso_5run_array,consensus=cvConsensus)
    # we still need the rowNames later, it doesn't change in the for loop above, so we can take any iteratant inside.
    output_rowName=rownames(coefB_lasso_cvmin)
    
    print(paste0('lasso for gene',geneID,' ends at',Sys.time()))
    coefB_lasso_cvmin_5run=data.frame(coefB_lasso_cvmin_5run)
    names(coefB_lasso_cvmin_5run)=c('A','T','N'); coefB_lasso_cvmin_5run$probe=output_rowName

    coefB_lasso_cv1se_5run=data.frame(coefB_lasso_cv1se_5run)
    names(coefB_lasso_cv1se_5run)=c('A','T','N'); coefB_lasso_cv1se_5run$probe=output_rowName

   
    write.csv(coefB_lasso_cvmin_5run[,c('probe','A','T','N')],paste0('Output_lasso_chr_',chr_num,'_gene_',geneID,'_useCVmin_betaOutputs.csv'),quote=F,row.names=F)
    write.csv(coefB_lasso_cv1se_5run[,c('probe','A','T','N')],paste0('Output_lasso_chr_',chr_num,'_gene_',geneID,'_useCV1se_betaOutputs.csv'),quote=F,row.names=F)
    print(paste0('All runs for gene ', geneID,' is completed.'))
}

print(paste0('All genes in chr ',chr_num,' chunk ',chunk,' finished!'))

