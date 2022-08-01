library(lmerSeq)
library(DESeq2)
library(limma)
library(dplyr)
library(pbapply)

data_path = "/path/to/two_timepoint/simulated/datasets/"
results_path <- "/path/to/save/results/"


sim_type_name <- "ROC_MIX" 
param_mat <- expand.grid(n_sub_group = c(3, 5, 10, 20))


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
    
    ###   getting vst data
    dds <- DESeqDataSetFromMatrix(countData = cnts,
                                  colData = smd,
                                  design = ~ group * time)
    dds_vst = varianceStabilizingTransformation(dds, blind = F, fitType = "parametric")
    expr_vst = assay(dds_vst)
    res_vst = lmerSeq.fit(form = ~ group * time + (1|id), 
                          expr_mat = expr_vst, 
                          sample_data = smd, 
                          parallel = F, 
                          cores = 6)
    
    tt_int_vst = lmerSeq.contrast(lmerSeq_results = res_vst, 
                                  contrast_mat = rbind(c(0, 0, 0, 1)),
                                  p_adj_method = "BH",
                                  sort_results = F, 
                                  include_singular = T)
    
    tt_btw_vst = lmerSeq.contrast(lmerSeq_results = res_vst, 
                                  contrast_mat = rbind(c(0, 1, 0, 1)),
                                  p_adj_method = "BH",
                                  sort_results = F, 
                                  include_singular = T)
    
    tt_wtn_vst = lmerSeq.contrast(lmerSeq_results = res_vst, 
                                  contrast_mat = rbind(c(0, 0, 1, 1)),
                                  p_adj_method = "BH",
                                  sort_results = F, 
                                  include_singular = T)
    
    vst_tab_out_int = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      left_join(tt_int_vst$summary_table) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Interaction')
    
    vst_tab_out_btw = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      left_join(tt_btw_vst$summary_table) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Between')
    
    vst_tab_out_wtn = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      left_join(tt_wtn_vst$summary_table) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Within')
    
    out_tab_vst = rbind(vst_tab_out_int, vst_tab_out_btw, vst_tab_out_wtn)
    
    
    ret = list(vst_tab_out = out_tab_vst)
    return(ret)
  })
  
  
  res_all_vst <- do.call(rbind, lapply(sim_results, function(x){return(x$vst_tab_out)}))
  
  res_both = list(res_all_vst = res_all_vst)
  
  file_out <- paste0("just_vst_w_singular_ROC_MIX_n_", n_sample_per_group, ".RDS")
  saveRDS(object = res_both,
          file = paste0(results_path,
                        file_out))
}