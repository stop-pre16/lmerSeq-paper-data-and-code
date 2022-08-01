library(lmerSeq)
library(DESeq2)
library(limma)
library(dplyr)
library(pbapply)
library(variancePartition)

# data_path <- "~/Documents/longitudinal_rna_seq/sims_based_real_data/sim_data_cj/"
data_path = "/path/to/two_timepoint/simulated/datasets/"
results_path <- "/path/to/save/results/"


sim_type_name <- "ROC_MIX" 
param_mat <- expand.grid(n_sub_group = c(3, 5, 10, 20))

options(MulticoreParam=quote(MulticoreParam(workers=1)))
options(mc.cores=1)

serialParam <- SerialParam()

for(row_id in 1:4){
  # row_id <- 1
  
  n_sample_per_group <- param_mat$n_sub_group[row_id]
  n_sim <- 10
  
  sim_results <- pblapply(1:n_sim, function(j){
    # j = 1
    sim_data <- readRDS(paste0(data_path, "sim_data_", sim_type_name, "_n_",
                               n_sample_per_group,"_",j, ".RDS"))
    
    cnts = sim_data$counts
    smd = data.frame(group = sim_data$groups, 
                     time = sim_data$time, 
                     id = sim_data$ids)
    
    expr_dream = voomWithDreamWeights(counts = cnts, 
                                      formula = ~ group * time + (1|id),
                                      data = smd,
                                      BPPARAM = serialParam)
    
    L_mat = cbind(c(0, 0, 0, 1),
                  c(0, 1, 0, 1),
                  c(0, 0, 1, 1))
    
    res_dream_cont = dream(exprObj = expr_dream, 
                           formula = ~ group * time + (1|id),  
                           L = L_mat,
                           BPPARAM = serialParam,
                           data = smd)
    
    tt_int_dream = topTable(fit = res_dream_cont, 
                            coef = "L1", 
                            number = Inf, 
                            sort.by = 'none')
    head(tt_int_dream)
    
    tt_btw_dream = topTable(fit = res_dream_cont, 
                            coef = "L2", 
                            number = Inf, 
                            sort.by = 'none')
    head(tt_btw_dream)
    
    tt_wtn_dream = topTable(fit = res_dream_cont, 
                            coef = "L3", 
                            number = Inf, 
                            sort.by = 'none')
    head(tt_wtn_dream)
    
    dream_tab_out_int = cbind(sim_data$param, tt_int_dream) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      # left_join(tt_int_dream$summary_table) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Interaction')
    
    dream_tab_out_btw = cbind(sim_data$param, tt_btw_dream) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      # left_join(tt_btw_dream$summary_table) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Between')
    
    dream_tab_out_wtn = cbind(sim_data$param, tt_wtn_dream) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      # left_join(tt_wtn_dream$summary_table) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Within')
    
    out_tab_dream = rbind(dream_tab_out_int, dream_tab_out_btw, dream_tab_out_wtn)
    
    ret = list(dream_tab_out = out_tab_dream)
    return(ret)
  })
  
  res_all_dream <- do.call(rbind, lapply(sim_results, function(x){return(x$dream_tab_out)}))
  
  res_both = list(res_all_dream = res_all_dream)
  
  file_out <- paste0("just_dream_ROC_MIX_n_", n_sample_per_group, ".RDS")
  saveRDS(object = res_both,
          file = paste0(results_path,
                        file_out))
  
}
