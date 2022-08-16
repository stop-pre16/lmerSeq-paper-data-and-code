library(dplyr)
library(pbapply)
library(rmRNAseq)

data_path <- "/path/to/simulated/data/"
save_path = "/path/to/results/folder/"

# Simualtion type
sim_type_name <- "ROC" # "FDR" "ROC" "T1E" "POW"

param_mat <- expand.grid(n_sub_group = c(3, 5, 10, 20))

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
                     id = sim_data$ids)
    smd$time_factor <- factor(smd$time)
    dmat = model.matrix(~ smd$group * smd$time)
    
    Cmat = list()
    Cmat[[1]] = cbind(c(0,0,0,1))
    Cmat[[2]] = cbind(c(0,1,0,3))
    Cmat[[3]] = cbind(c(0,0,1,1))
    Cmat[[4]] = cbind(c(0,1,0,0),
                      c(0,0,1,0),
                      c(0,0,0,1))
    names(Cmat) = c("Interaction", "Between", "Within", "Overall")
    
    res_rmrnaseq = TC_CAR1(counts = cnts[1:100, ], 
                           design = dmat, 
                           Subject = as.factor(smd$id), 
                           Time = smd$time, 
                           C.matrix = Cmat,
                           Nboot = 100,
                           ncores = 8, 
                           print.progress = T)
    
    # res_rmrnaseq$pqvalue$qv
    
    tab_out_int = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      mutate(p_val_raw = as.numeric(res_rmrnaseq$pqvalue$pv$Interaction)) %>% 
      mutate(q_val =  as.numeric(res_rmrnaseq$pqvalue$qv$Interaction)) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Interaction')
    
    tab_out_btw = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      mutate(p_val_raw = as.numeric(res_rmrnaseq$pqvalue$pv$Between)) %>% 
      mutate(q_val =  as.numeric(res_rmrnaseq$pqvalue$qv$Between)) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Between')
    
    tab_out_wtn = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      mutate(p_val_raw = as.numeric(res_rmrnaseq$pqvalue$pv$Within)) %>% 
      mutate(q_val =  as.numeric(res_rmrnaseq$pqvalue$qv$Within)) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Within')
    
    tab_out_ovr = data.frame(sim_data$param) %>% 
      mutate(gene = as.factor(1:nrow(sim_data$param))) %>% 
      mutate(true_de = beta_ints != 0) %>% 
      mutate(p_val_raw = as.numeric(res_rmrnaseq$pqvalue$pv$Overall)) %>% 
      mutate(q_val =  as.numeric(res_rmrnaseq$pqvalue$qv$Overall)) %>% 
      mutate(sim_num = j, n_per_grp = n_sample_per_group, contrast = 'Overall')
    
    out_tab_all = rbind(tab_out_int, tab_out_btw, tab_out_wtn, tab_out_ovr)
    
    ret = list(all_tab_out = out_tab_all)
    return(ret)
  })
  
  res_all <- do.call(rbind, lapply(sim_results, function(x){return(x$all_tab_out)}))
  
  res_both = list(res_all = res_all)
  
  file_out <- paste0("just_rmrnaseq_4tp_cont_n_", n_sample_per_group, ".RDS")
  saveRDS(object = res_both,
          file = paste0(save_path,
                        file_out))
  
}