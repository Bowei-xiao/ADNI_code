# This is the script used to clean the ADNI DNAm data and phenotype matching


#### Part I: basic QC
# iDAT files to input and out put two R objects storing the methylation proportion (beta) and detection p-value (pvals)
# Both are stored in this RData file
load('ADNI_beta-pval_matrices.RData')

# consider detection p-value >=0.01 as fail, check fail rate across samples and probes
pval_threshold=0.01
pvals[pvals>=pval_threshold]=NA
miss_byRow=apply(pvals,1,function(t) sum(is.na(t))/length(t))
miss_byColumn=apply(pvals,2,function(t) sum(is.na(t))/length(t))

# check failed rate
if (F) {
png('hist_missingness_byProbe.png')
hist(miss_byRow,xlab='missing rate')
dev.off()

png('hist_missingness_bySample.png')
hist(miss_byColumn,xlab='missing rate')
dev.off()
}

# This chip is EPIC/850k
# On the methylation data, there are some SNPs and cross-reaction probes that should be removed
library('maxprobes')
xloci <- maxprobes::xreactive_probes(array_type = "EPIC")
xloci_vec=unlist(xloci)
# remove these and then check missing rate again
remove_rowinx=which(rownames(pvals) %in% xloci_vec)
pvals_sub=pvals[-remove_rowinx,]
# row is by probe
miss_byRow=apply(pvals_sub,1,function(t) sum(is.na(t))/length(t))
# column is by sample
miss_byColumn=apply(pvals_sub,2,function(t) sum(is.na(t))/length(t))

# remove probes where missing exists on 20% of all samples
remove_rowinx2=names(miss_byRow[miss_byRow>=0.2])
pvals_sub2=pvals_sub[which(!rownames(pvals_sub) %in% remove_rowinx2),]
miss_byColumn=apply(pvals_sub2,2,function(t) sum(is.na(t))/length(t))
# no need to remove individuals
# output the subsets
keep_row=rownames(pvals)[rownames(pvals) %in% rownames(pvals_sub2)]
write.table(keep_row,'probeLists_afterQC.txt',quote = F,col.names = F,row.names = F)


# match phenotype IDs (and calculate relevant covariates)
# firstly, calculate epi-age by Horvath et al.
#### Part II: calculate epi-age.
library('methylclock')
library('tibble')
library('ExperimentHub')

ExperimentHub::setExperimentHubOption("LOCAL", TRUE)
eh=ExperimentHub()
# subset to the probes passed QC
probeList=keep_row

probeList=read.table('probeLists_afterQC.txt',stringsAsFactors=F)[,1]

subProbe=betas[rownames(betas) %in% probeList,]

# change the rownames to be the first column of the dataset to follow Horvath's format
# release some space
rm(betas); rm(pvals); gc()
# convert subProbe into a tibble format to reduce space/memory useage
subProbe_tbl=as_tibble(subProbe)
# add the first column: probe ID
subProbe_tbl['ProbeID']=rownames(subProbe)
subProbe_tbl=subProbe_tbl %>% relocate('ProbeID')
save(subProbe_tbl,file='ADNI_probes_afterQC_tbl.RData')

# run this by batch for computational feasibility
batch_num=50
# run this parallele
n_job=60;
# so each one should take indices from [n/60*(i-1)+1,n/60*i]+1 (since index starts from 2)
# in case this overflows/underflows, hard cap min and max index
min_inx=max((ncol(subProbe_tbl)/n_job*(batch_num-1)+1),1)
max_inx=min((ncol(subProbe_tbl)/n_job*batch_num),ncol(subProbe_tbl)-1)

subProbe_batch=subProbe_tbl[,c(1,(min_inx:max_inx)+1)]
print(Sys.time())
res1=DNAmAge(subProbe_batch,clocks="Horvath",cell.count=F,normalize=T) 
print(Sys.time())
write.csv(res1,paste0('Horvath_mAge_normalized_batch_',as.character(batch_num),'.txt'),quote=F,row.names = F)


#### Part III: match phenotype data to the DNAm data and obtain important covariates
# plate info from ADNI portal
plate_info=read.csv('ADNI_DNA_Methylation_SampleAnnotation_20170530_sample-annotations_BG.csv',header=T,stringsAsFactors=F)

# load covariates from covariate files
load('20220128_all-vars-env_variables.RData')
# Here we will ignore the fact that people get multiple methylation. Will for now assume each entry is one individual.
# useful columns?
pheno_file=annot.ext[,c('RID.a','barcodes','prop.B','prop.NK','prop.CD4T','prop.CD8T','prop.Mono','prop.Neutro','DX','age.now','PTGENDER','PTEDUCAT','APOE4','ABETA','TAU','PTAU')]

# For immune cells, we calculate its PCs and includes the first 3 PCs as covariates
PCs=prcomp(pheno_file[,c('prop.B','prop.NK','prop.CD4T','prop.CD8T','prop.Mono','prop.Neutro')])
pheno_file=cbind(pheno_file, data.frame(PCs$x[,1:3]))

# change DX to numeric. CN& MCI =0 and Dementia=1
pheno_file$DX_num=0;
pheno_file$DX_num[pheno_file$DX=='Dementia']=1


# add a column to mark which entry should be used as the unique ID 
# if one individual appears multiple time, take the earliest viist (youngest age) will be marked as 1
pheno_file$uniqueID=0
uniqueID=c()
for (i in 1:nrow(pheno_file)){
  phenoID=pheno_file$RID.a[i]
  if (phenoID %in% uniqueID){
    # this is a duplicated individual, choose the youngest age and change the uniqueID entry to 1
    old_index=which(pheno_file$uniqueID==1 & pheno_file$RID.a==phenoID)
    age_old=pheno_file$age.now[old_index]
    # if the current age is younger, this will be selected
    #however, if a/t/n is missing then don't change
    if ((age_old>pheno_file$age.now[i])& (!is.na(pheno_file$ABETA[i])) & (!is.na(pheno_file$TAU[i])) & (!is.na(pheno_file$PTAU[i]))){
      pheno_file$uniqueID[i]=1; pheno_file$uniqueID[old_index]=0
    }
  } else {
    # new individual, set uniqueID to 1
    pheno_file$uniqueID[i]=1
    # add this ID to uniqueID list
    uniqueID=c(uniqueID,phenoID)
  }
}

# bind the plate info then output
pheno_withPlateinfo=merge(pheno,plate_info[,c('barcodes','PlateNumber')])
write.csv(pheno_withPlateinfo,'ADNI_covariate_moreCovs_540samples_withEpiage_1905obs.csv',quote=F,row.names=F)
