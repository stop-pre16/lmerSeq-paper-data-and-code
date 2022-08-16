library(dplyr)
library(pbapply)
library(variancePartition)

data_path <- "/path/to/simulated/data/"
save_path = "/path/to/results/folder/"

# Simualtion type
sim_type_name <- "ROC" # "FDR" "ROC" "T1E" "POW"

param_mat <- expand.grid(n_sub_group = c(3, 5, 10, 20))
options(MulticoreParam=quote(MulticoreParam(workers=1)))
options(mc.cores=1)

serialParam <- SerialParam()

for(row_id in 1:4){
  # row_id = 1
  n_sample_per_group <- param_mat$n_sub_group[row_id]
  n_sim <- 10
  
  sim_results <- pblapply(1:n_sim, function(j){
    # j = 1
    sim_data <- readRDS(paste0(data_path, "sim_data_4tp_", sim_type_name, "_n_",
                               n_sample_per_group,"_",j, ".RDS"))
    
    cnts = sim_data$counts
    smd = data.frame(group = sim_data$groups, 
                     time = sim_data$time, 
                     subject = sim_data$ids)
    smd$time_factor <- factor(smd$time)
    
    expr_dream = voomWithDreamWeights(counts = cnts, 
                                      formula = ~ group * time_factor + (1|subject),
                                      data = smd,
                                      BPPARAM = serialParam)
    fit.lmerSeq.continuous_ri <- lmerSeq.fit(form = ~ group*time_factor + (1|subject),
                                             expr_mat = expr_dream$E,
                                             weights = expr_dream$weights,
                                                sample_data = smd,
                                                REML = T, cores = 3,
                                                parallel=T)
    # Difference in change over time
    interaction.ttest.cont <- rbind(c(0,0,0,0,0,0,0,1))
    
    cat_ri_ineraction.ttest <- lmerSeq.contrast(lmerSeq_results = fit.lmerSeq.continuous_ri,
                                                    contrast_mat = interaction.ttest.cont,
                                                    p_adj_method = 'BH',
                                                    ddf = 'Satterthwaite',
                                                    sort_results = F, include_singular=T)
    
    # Change over time in treatment group
    time.ttest.cont <- group.ttest.cat <- rbind(c(0,1,0,0,0,0,0,1))
    
    cat_ri_time.ttest <- lmerSeq.contrast(lmerSeq_results = fit.lmerSeq.continuous_ri,
                                              contrast_mat = time.ttest.cont,
                                              p_adj_method = 'BH',
                                              ddf = 'Satterthwaite',
                                              sort_results = F, include_singular=T)
    
    # Difference between groups at t = 3
    group.ttest.cont <- rbind(c(0,0,0,0,1,0,0,1))
    
    cat_ri_group.ttest <- lmerSeq.contrast(lmerSeq_results = fit.lmerSeq.continuous_ri,
                                               contrast_mat = group.ttest.cont,
                                               p_adj_method = 'BH',
                                               ddf = 'Satterthwaite',
                                               sort_results = F, include_singular=T)
    rm(fit.lmerSeq.continuous_ri)
    gc()
    gc()
    
    cat_ri_tab_out_int = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      left_join(cat_ri_ineraction.ttest$summary_table %>% dplyr::select(-lower, -upper)) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Interaction', fit = "cat_ri")
    
    cat_ri_tab_out_btw = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      left_join(cat_ri_group.ttest$summary_table %>% dplyr::select(-lower, -upper)) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Between', fit = "cat_ri")
    
    cat_ri_tab_out_wtn = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      left_join(cat_ri_time.ttest$summary_table %>% dplyr::select(-lower, -upper)) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Within', fit = "cat_ri")
    
    out_tab_vst = rbind(cat_ri_tab_out_int, cat_ri_tab_out_btw, cat_ri_tab_out_wtn)
    ret = list(vst_tab_out = out_tab_vst)
    return(ret)
  })
  
  
  res_all_vst <- do.call(rbind, lapply(sim_results, function(x){return(x$vst_tab_out)}))
  
  
  res_both = list(res_all_vst = res_all_vst)
  
  file_out <- paste0("DREAM_4tp_cat_ri_n_", n_sample_per_group, ".RDS")
  saveRDS(object = res_both,
          file = paste0(save_path,
                        file_out))
}