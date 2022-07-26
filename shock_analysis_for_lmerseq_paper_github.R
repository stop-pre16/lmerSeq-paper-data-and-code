#######################################
##
## Script name: CS SHOCK ANALYSIS FOR LMERSEQ PAPPER
##
## Purpose of script: Analysis of shock dataset using LMERSEQ, DREAM AND RMRNASEQ
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

#######################################
## 01 - Read, format, filter data
#######################################

############## Read in Data #################
dat <- readRDS('cs_dataset.RDS')
counts <- dat$counts
sample_data <- dat$sample_data

############## Filter #################
cpm_counts=edgeR::cpm(counts)
cpm_more1<-apply(cpm_counts, 1, function(x) sum(x>1))
idx_keep=which(cpm_more1>=11)
counts=counts[idx_keep,]

############## Methods we're running #################
methods=c("lmm","dream","rmrnaseq")

############## Contrasts #################
original_contrasts=list(
  c1 = c(0,1,0), # diff between 48 hours and baseline
  c2 = c(0,0,1), # diff between 1 week and baseline
  c3 = c(0,-1,1) # diff between 1 week and 48 hours
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
  file_name_fit=paste0("cs_shock_fit_", method,  ".RDS")
  file_name_sum=paste0("cs_shock_summary_", method,".RDS")
  
  
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
    
    serialParam <- SerialParam() # parallelization was not working on our HPC
    
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




###################################################################
# Summarize results
###################################################################

rmrnaseq <- readRDS('cs_shock_summary_rmrnaseq.RDS')
dream <- readRDS('cs_shock_summary_dream.RDS')
lmm <- readRDS('cs_shock_summary_lmm.RDS')

# Number of significant genes by contrast for each method
lmm_nc1 <- table(lmm$contrasts$c1$summary_table$p_val_adj < 0.05)['TRUE'] # 1
lmm_nc2 <- table(lmm$contrasts$c2$summary_table$p_val_adj < 0.05)['TRUE'] # 1452
lmm_nc3 <- table(lmm$contrasts$c3$summary_table$p_val_adj < 0.05)['TRUE'] # 362
lmm_n_sing <- length(lmm$contrasts$c3$singular_fits)

rmrnaseq_nc1 <- table(rmrnaseq$contrasts$c1$p_val_adj < 0.05)['TRUE']
rmrnaseq_nc2 <- table(rmrnaseq$contrasts$c2$p_val_adj < 0.05)['TRUE']
rmrnaseq_nc3 <- table(rmrnaseq$contrasts$c3$p_val_adj < 0.05)['TRUE']

dream_nc1 <- table(p.adjust(dream$contrasts$p.value[,"L1"], method='BH')< 0.05)['TRUE']
dream_nc2 <- table(p.adjust(dream$contrasts$p.value[,"L2"], method='BH')< 0.05)['TRUE']
dream_nc3 <- table(p.adjust(dream$contrasts$p.value[,"L3"], method='BH')< 0.05)['TRUE']


# Test for normality and heteroskedasticity (DREAM and lmerSeq)
library(car)
lmm_fit <- readRDS('cs_shock_fit_lmm.RDS')

# Levine's test for heteroskedasticity
lt <- lapply(lmm_fit, function(x) leveneTest(residuals(x$fit, scale=T)~factor(sample_data$time))$`Pr(>F)`[1])
table(unlist(lt) < 0.05) #206 
table(p.adjust(unlist(lt), method='BH')<0.05) #0

lt_dream <- apply(dream$residuals, 1,function(x) leveneTest(x~factor(sample_data$time))$`Pr(>F)`[1])
table(unlist(lt_dream) < 0.05) #149 
table(p.adjust(unlist(lt_dream), method='BH')<0.05) #0

# KS test for normality
ks <- lapply(lmm_fit, function(x) ks.test(residuals(x$fit, scale=T), 'pnorm')$p.value)
table(unlist(ks) < 0.05) #24 
table(p.adjust(unlist(ks), method='BH')<0.05) #0

# need to scale the dream residuals
dream_scaled_res <- dream_$residuals

for(i in 1:nrow(dream_scaled_res)){
  dream_scaled_res[i,] <- dream_scaled_res[i,]/dream$sigma[i]
}

ks_dream <- apply(dream_scaled_res,1, function(x) ks.test(x, 'pnorm')$p.value)
table(unlist(ks_dream) < 0.05) #17
table(p.adjust(unlist(ks_dream), method='BH')<0.05) #0

# To get number of singular fits for DREAM, refit with lmerSeq
# Voom transformation
expr_dream = voomWithDreamWeights(counts = counts, 
                                  formula = ~ time + (1|patient),
                                  data = sample_data,
                                  BPPARAM = serialParam)


fit=lmerSeq.fit(form = ~time+(1|patient), expr_mat=expr_dream$E,
                weights = expr_dream$weights,
                sample_data=sample_data,
                parallel = T, cores = 8)

all_contrasts=lapply(contrasts, function(contrast){
  sum=lmerSeq.contrast(fit, contrast=rbind(contrast),sort_results = F)
  sum
})
names(all_contrasts)=names(contrasts)
dream_n_sing <- length(all_contrasts$contrasts$c3$singular_fits)


##############################################
# Enrichment for lmerSeq C2 DEGs
#############################################
library(enrichR)
library(biomaRt)

# Get gene symbols for ensembl IDs
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(counts)

G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

# Get lmerSeq up and down DEGs
lmm_up <- lmm$contrasts$c2$summary_table[ lmm$contrasts$c2$summary_table$p_val_adj < 0.05 & lmm$contrasts$c2$summary_table$Estimate > 0,]$gene
lmm_down <- lmm$contrasts$c2$summary_table[ lmm$contrasts$c2$summary_table$p_val_adj < 0.05 & lmm$contrasts$c2$summary_table$Estimate < 0,]$gene

symbols_up <- G_list[G_list$ensembl_gene_id %in% lmm_up,]$hgnc_symbol
symbols_down <- G_list[G_list$ensembl_gene_id %in% lmm_down,]$hgnc_symbol

# Perform enrichment analysis
dbs <- c("BioPlanet_2019")

tmp <- enrichr(symbols_up, dbs)
enrichr_list_up <- do.call(rbind, tmp)
enrichr_list_up <- enrichr_list_up[enrichr_list_up$Adjusted.P.value < 0.05,]
enrichr_list_up <- enrichr_list_up[order(enrichr_list_up$Adjusted.P.value),]

tmp <- enrichr(symbols_down, dbs)
enrichr_list_down <- do.call(rbind, tmp)
enrichr_list_down <- enrichr_list_down[enrichr_list_down$Adjusted.P.value < 0.05,]
enrichr_list_down <- enrichr_list_down[order(enrichr_list_down$Adjusted.P.value),]

write.xlsx(list(up = enrichr_list_up, 
                down = enrichr_list_down),
           'lmerseq_shock_enrichment_bioplanet.xlsx')











