#######################################
##
## Script name: PERMUTATION CS SHOCK ANALYSIS FOR LMERSEQ PAPPER
##
## Purpose of script: CREATE AND ANALYZE SIMULATED DATASETS BASED ON PERMUTATIONS OF THE SHOCK DATA
##
## Author: Camille Moore and Elizabeth Wynn
##
## Date Created: 2022-07-15
##
#######################################
##
## Notes:
##   
##
#######################################

## load the packages we will need:
library(lmerSeq)
library(DESeq2)
library(dplyr)
library(rmRNAseq)
library(variancePartition)
library(BiocParallel)
library(openxlsx)
library(edgeR)

############################################################
## 01 - Read in sample data and create 10 permuted datasets
#############################################################

############## Read in Data #################
dat <- readRDS('cs_dataset.RDS')
counts_raw <- dat$counts
sample_data <- dat$sample_data

############## Permute #################
permuted_sample_data <- list()

for(i in 1:10){
  temp <- sample_data
  
  # Remove a subject so we will have a total sample size of 10
  drop_subject <- sample(size=1, x = unique(temp$patient)) #choose 1 to drop 
  temp <- temp[temp$patient != drop_subject,]
  
  # Flip the time points for 5 subjects
  shuffle_subject <- sample(size = 5, unique(temp$patient)) # choose 5 to shuffle
  temp$time_new <- temp$time
  temp[temp$patient %in% shuffle_subject & temp$time == 'T1',]$time_new <- 'T3'
  temp[temp$patient %in% shuffle_subject & temp$time == 'T3',]$time_new <- 'T1'
  temp$time <- temp$time_new
  permuted_sample_data[[i]] <-temp[,c('sample', 'patient', 'time', 'group')]
}

saveRDS(permuted_sample_data, 'permuted_sample_data.RDS')


for(perm in 1:10){
  
  ############## Read in Data #################
  sample_data<-readRDS("permuted_sample_data.RDS")[[perm]]
  counts=counts_raw[,rownames(sample_data)]
  
  ############## Filter #################
  cpm_counts=edgeR::cpm(counts)
  cpm_more1<-apply(cpm_counts, 1, function(x) sum(x>1))
  idx_keep=which(cpm_more1>=11)
  counts=counts[idx_keep,]
  
  
  ########### Adjust the counts to make DE for 2600 genes (1300 up and 1300 down)
  # Genes 1:1302 are up at time 3
  counts[1:1300, colnames(counts) %in% sample_data$sample[sample_data$time=='T3']] <- 2*counts[1:1300, colnames(counts) %in% sample_data$sample[sample_data$time=='T3']]
  counts[10000:11300, colnames(counts) %in% sample_data$sample[sample_data$time=='T3']] <- round(0.5*counts[10000:11300, colnames(counts) %in% sample_data$sample[sample_data$time=='T3']])
  
  
  ############## Methods we're running #################
  methods=c("lmm","dream","rmrnaseq")

  
  ############## Contrasts #################
  original_contrasts=list(
    ## Difference between timepoints in CS group
    c1 = c(0,1,0),
    c2 = c(0,0,1),
    c3 = c(0,-1,1)
  ) 
  
  #######################################
  ## 02 - Run/save analysis
  #######################################
  
  for(method in methods){
    contrasts=original_contrasts
    sim_start=Sys.time()
    print(method)
    print("Fitting Models")
    
    ############## File names for saving results #################
    file_name_fit=paste0("cs_shock_fit_", method,  "_perm_spike_", perm, ".RDS")
    file_name_sum=paste0("cs_shock_summary_", method,"_perm_spike_", perm, ".RDS")
    
    
    ############## Run LMM #################
    if(method=="lmm"){
      
      ############## Transform Data #################
      dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = sample_data,
                                    design = ~ time)
      dds <- DESeq(dds)
      vsd.fixed <- varianceStabilizingTransformation(dds, blind=F)
      counts_trans<- assay(vsd.fixed)
      rm(dds)
      rm(vsd.fixed)
      
      ############## Fit/save models #################
      fit=lmerSeq.fit(form = ~time+(1|patient), expr_mat=counts_trans,
                      sample_data=sample_data,
                      parallel = T, cores = 8)
      saveRDS(fit, file=file_name_fit)
      sim_end_fit=Sys.time()
      time_fit=sim_end_fit-sim_start
      
      ############## Run hypothesis testing on all contrasts #################
      print("Summarizing")
      
      all_contrasts=lapply(contrasts, function(contrast){
        sum=lmerSeq.contrast(fit, contrast=rbind(contrast),sort_results = F)
        sum
      })
      names(all_contrasts)=names(contrasts)
      #all_contrasts[['time_test']] <- lmerSeq.contrast(fit, contrast=do.call(rbind, contrasts[1:2]), sort_results = F)
      sim_end=Sys.time()
      rm(fit)
      
      
      ##############rmRNAseq##################
      
    }else if(method=="rmrnaseq"){
      ############## Define model matrix/reformat contrasts #################
      dmat=model.matrix(~time, data=sample_data)
      contrasts_adj=lapply(lapply(contrasts, rbind), t)
      
      ############## Run models and testing #################
      fit = TC_CAR1(counts = counts,
                    design = dmat,
                    Subject = as.factor(sample_data$patient),
                    Time = as.numeric(factor(sample_data$time)),
                    C.matrix = contrasts_adj,
                    Nboot = 100,
                    ncores = 8,
                    print.progress = T)
      sim_end=Sys.time()
      time_fit=NA
      saveRDS(fit, file=file_name_fit)
      
      ############## Get adjusted p-values #################
      all_contrasts=lapply(paste0("c", 1:3), function(x){
        p_val_adj<-p.adjust(fit$pqvalue$pv[,x], method="BH")
        data.frame(p_val_raw=fit$pqvalue$pv[,x], p_val_adj=p_val_adj, qval=fit$pqvalue$qv[,x])
      })
      names(all_contrasts)=paste0("c", 1:3)
      rm(fit)
      
      
    }else if(method=='dream'){
      
      serialParam <- SerialParam() # this is important, keep this line
      
      # Voom transformation
      expr_dream = voomWithDreamWeights(counts = counts, 
                                        formula = ~ time + (1|patient),
                                        data = sample_data,
                                        BPPARAM = serialParam)
      
      L_mat = cbind(c(0, 1,0), c(0,0,1), c(0,-1,1))
      
      all_contrasts = dream(exprObj = expr_dream, 
                            formula = ~ time + (1|patient),  
                            L = L_mat,
                            BPPARAM = serialParam,
                            data = sample_data)
      
      sim_end=Sys.time()
      time_fit=NA
      
    }
    
    ############## Save summary results  #################
    print("Summary")
    sum=list(contrasts=all_contrasts, time_total=sim_end-sim_start, time_fit=time_fit)
    rm(all_contrasts)
    saveRDS(sum, file_name_sum)
    rm(sum)
  }
  
  
}


#################################################
# Summarize results after all analyses complete
#################################################

res <- NULL

for (i in c(1:10)){
  
  dat_lmm <- readRDS(paste0("cs_shock_summary_lmm_perm_spike_",i,".RDS"))
  
  lmm_n_fp <- table(dat_lmm$contrasts$c2$summary_table[-c(1:1300, 10000:11300),]$p_val_adj < 0.05)['TRUE'] #false positives
  lmm_n_tn <- table(dat_lmm$contrasts$c2$summary_table[-c(1:1300, 10000:11300),]$p_val_adj < 0.05)['FALSE'] # true negatives
  
  lmm_n_tp <- table(dat_lmm$contrasts$c2$summary_table[c(1:1300, 10000:11300),]$p_val_adj < 0.05)['TRUE'] # true positives
  lmm_n_fn <- table(dat_lmm$contrasts$c2$summary_table[c(1:1300, 10000:11300),]$p_val_adj < 0.05)['FALSE'] # false negatives
  
  dat_dream <- readRDS(paste0("cs_shock_summary_dream_perm_spike_",i,".RDS"))
  
  dream_p_adj <- p.adjust(dat_dream$contrasts$p.value[,'L2'], method='BH')
  
  dream_n_fp <- table(dream_p_adj[-c(1:1300, 10000:11300)] < 0.05)['TRUE']
  dream_n_tn <- table(dream_p_adj[-c(1:1300, 10000:11300)] < 0.05)['FALSE']
  
  dream_n_tp <- table(dream_p_adj[c(1:1300, 10000:11300)] < 0.05)['TRUE']
  dream_n_fn <- table(dream_p_adj[c(1:1300, 10000:11300)] < 0.05)['FALSE']
  
  
  dat_rm <- readRDS(paste0("cs_shock_summary_rmrnaseq_perm_spike_",i,".RDS"))
  
  rm_n_fp <- table(dat_rm$contrasts$c2[-c(1:1300, 10000:11300),]$p_val_adj < 0.05)['TRUE']
  rm_n_tn <- table(dat_rm$contrasts$c2[-c(1:1300, 10000:11300),]$p_val_adj < 0.05)['FALSE']
  
  rm_n_tp <- table(dat_rm$contrasts$c2[c(1:1300, 10000:11300),]$p_val_adj < 0.05)['TRUE']
  rm_n_fn <- table(dat_rm$contrasts$c2[c(1:1300, 10000:11300),]$p_val_adj < 0.05)['FALSE']
  
  
  
  res <- rbind(res, data.frame(dream_n_fp = dream_n_fp, 
                               dream_n_tn = dream_n_tn, 
                               dream_n_tp = dream_n_tp, 
                               dream_n_fn = dream_n_fn,
                               lmm_n_fp = lmm_n_fp, 
                               lmm_n_tn = lmm_n_tn, 
                               lmm_n_tp = lmm_n_tp, 
                               lmm_n_fn = lmm_n_fn, 
                               rm_n_fp = rm_n_fp, 
                               rm_n_tn = rm_n_tn, 
                               rm_n_tp = rm_n_tp, 
                               rm_n_fn = rm_n_fn
  ))
  
}

res$dream_fdr <- res$dream_n_fp / (res$dream_n_fp + res$dream_n_tp)
res$dream_power <- res$dream_n_tp / (res$dream_n_fn + res$dream_n_tp)

res$lmm_fdr <- res$lmm_n_fp / (res$lmm_n_fp + res$lmm_n_tp)
res$lmm_power <- res$lmm_n_tp / (res$lmm_n_fn + res$lmm_n_tp)

res$rm_fdr <- res$rm_n_fp / (res$rm_n_fp + res$rm_n_tp)
res$rm_power <- res$rm_n_tp / (res$rm_n_fn + res$rm_n_tp)

colMeans(res)

write.csv(res,'shock_perm_spike_result_summary.csv')



