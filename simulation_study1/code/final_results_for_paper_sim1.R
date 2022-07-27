library(dplyr)
library(ggplot2)
library(reshape2)
library(pbapply)


n_vec = c(3, 5, 10, 20)
res_path = "/path/to/two_timepoint_sim/results/files/"

####   Looking at type 1 error rates    ####

all_t1e = bind_rows(pblapply(n_vec, function(n_sub){
  # n_sub = 3
  
  file_path2 = paste0(res_path, "just_dream_ROC_MIX_n_", n_sub, ".RDS")
  res_dream2 = readRDS(file_path2)$res_all_dream
  
  file_path3 = paste0(res_path, "just_rmrnaseq_ROC_MIX_n_", n_sub, ".RDS")
  res_rmrnaseq = readRDS(file_path3)$res_all
  
  file_path4 = paste0(res_path, "just_vst_w_singular_ROC_MIX_n_", n_sub, ".RDS")
  res_vst_singular = readRDS(file_path4)$res_all_vst
  
  file_path5 = paste0(res_path, "just_vst_w_gls_ROC_MIX_n_", n_sub, ".RDS")
  res_gls = readRDS(file_path5)$res_all_gls
  p_val_levels = c(0.05, 0.01, 0.001, 0.0001)
  
  sim_ids = unique(res_vst_singular$sim_num)
  
  t1e_vst_singular = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_vst_singular %>% filter(sim_num == sim_sub)
    t1e_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      t1e_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val_raw < p_sub, levels = c(F, T)))
        t1e_est = tab_sub[1, 2] / sum(tab_sub[1, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         t1e = t1e_est, rel_t1e = t1e_est / p_sub, transform = 'lmerSeq (WS)', 
                         contrast = c_sub)
      }))
      return(t1e_sub)
    }))
    
  }))
  t1e_vst_singular
  
  t1e_dream2 =bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_dream2 %>% filter(sim_num == sim_sub)
    t1e_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      t1e_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$P.Value < p_sub, levels = c(F, T)))
        t1e_est = tab_sub[1, 2] / sum(tab_sub[1, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         t1e = t1e_est, rel_t1e = t1e_est / p_sub, transform = 'DREAM (WS)', 
                         contrast = c_sub)
      }))
      return(t1e_sub)
    }))
    
  }))
  t1e_dream2
  
  t1e_rmrnaseq = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_rmrnaseq %>% filter(sim_num == sim_sub)
    t1e_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      t1e_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val_raw < p_sub, levels = c(F, T)))
        t1e_est = tab_sub[1, 2] / sum(tab_sub[1, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         t1e = t1e_est, rel_t1e = t1e_est / p_sub, transform = 'rmRNAseq', 
                         contrast = c_sub)
      }))
      return(t1e_sub)
    }))
    
  }))
  t1e_rmrnaseq
  
  t1e_gls =bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_gls %>% filter(sim_num == sim_sub)
    res_sub_true = res_vst_singular %>% filter(sim_num == sim_sub)
    t1e_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      res_sub_true2 = res_sub_true %>% filter(contrast == c_sub)
      t1e_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub_true2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val < p_sub, levels = c(F, T)))
        t1e_est = tab_sub[1, 2] / sum(tab_sub[1, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         t1e = t1e_est, rel_t1e = t1e_est / p_sub, transform = 'GLS', 
                         contrast = c_sub)
      }))
      return(t1e_sub)
    }))
    
  }))
  t1e_gls
  
  t1e_all = rbind(t1e_dream2, t1e_rmrnaseq, t1e_vst_singular, t1e_gls)
  return(t1e_all)
}))


t1e_summary = all_t1e %>% 
  group_by(transform, contrast, level, n_per_group) %>% 
  summarize(med_t1e = median(t1e), med_rel_t1e = median(rel_t1e),
            mean_t1e = mean(t1e), mean_rel_t1e = mean(rel_t1e))
t1e_summary$level = factor(t1e_summary$level, levels = c(0.0001, 0.001, 0.01, 0.05))

ggplot(t1e_summary) +
  aes(x = level, y = log2(mean_rel_t1e), color = transform, group = transform) +
  geom_point(size = 2) +
  # ylim(c(0, NA)) +
  facet_grid(rows = vars(contrast), cols = vars(n_per_group)) +
  geom_path() +
  geom_hline(yintercept = 0)


####   Looking at false discovery rates    ####

all_fdr = bind_rows(pblapply(n_vec, function(n_sub){
  # n_sub = 3
  
  file_path2 = paste0(res_path, "just_dream_ROC_MIX_n_", n_sub, ".RDS")
  res_dream2 = readRDS(file_path2)$res_all_dream
  
  file_path3 = paste0(res_path, "just_rmrnaseq_ROC_MIX_n_", n_sub, ".RDS")
  res_rmrnaseq = readRDS(file_path3)$res_all
  
  file_path4 = paste0(res_path, "just_vst_w_singular_ROC_MIX_n_", n_sub, ".RDS")
  res_vst_singular = readRDS(file_path4)$res_all_vst
  
  file_path5 = paste0(res_path, "just_vst_w_gls_ROC_MIX_n_", n_sub, ".RDS")
  res_gls = readRDS(file_path5)$res_all_gls
  
  p_val_levels = c(0.10, 0.05, 0.01)
  
  sim_ids = unique(res_vst_singular$sim_num)
  
  fdr_vst_singular = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_vst_singular %>% filter(sim_num == sim_sub)
    fdr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      fdr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val_adj < p_sub, levels = c(F, T)))
        fdr_est = tab_sub[1, 2] / sum(tab_sub[, 2])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         fdr = fdr_est, rel_fdr = fdr_est / p_sub, transform = 'lmerSeq (WS)', 
                         contrast = c_sub)
      }))
      return(fdr_sub)
    }))
    
  }))
  fdr_vst_singular
  
  fdr_gls = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_gls %>% filter(sim_num == sim_sub)
    res_sub_true = res_vst_singular %>% filter(sim_num == sim_sub)
    fdr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      res_sub_true2 = res_sub_true %>% filter(contrast == c_sub)
      fdr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub_true2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val_adj < p_sub, levels = c(F, T)))
        fdr_est = tab_sub[1, 2] / sum(tab_sub[, 2])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         fdr = fdr_est, rel_fdr = fdr_est / p_sub, transform = 'GLS', 
                         contrast = c_sub)
      }))
      return(fdr_sub)
    }))
    
  }))
  fdr_gls
  
  fdr_dream2 = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_dream2 %>% filter(sim_num == sim_sub)
    fdr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      fdr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$adj.P.Val < p_sub, levels = c(F, T)))
        fdr_est = tab_sub[1, 2] / sum(tab_sub[, 2])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         fdr = fdr_est, rel_fdr = fdr_est / p_sub, transform = 'DREAM (WS)', 
                         contrast = c_sub)
      }))
      return(fdr_sub)
    }))
    
  }))
  
  fdr_rmrnaseq = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_rmrnaseq %>% filter(sim_num == sim_sub)
    fdr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      fdr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$q_val < p_sub, levels = c(F, T)))
        fdr_est = tab_sub[1, 2] / sum(tab_sub[, 2])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         fdr = fdr_est, rel_fdr = fdr_est / p_sub, transform = 'rmRNAseq', 
                         contrast = c_sub)
      }))
      return(fdr_sub)
    }))
    
  }))
  
  fdr_rmrnaseq_bh = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_rmrnaseq %>% filter(sim_num == sim_sub) %>% mutate(p_val_bh = p.adjust(p_val_raw, method = "BH"))
    fdr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      fdr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val_bh < p_sub, levels = c(F, T)))
        fdr_est = tab_sub[1, 2] / sum(tab_sub[, 2])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         fdr = fdr_est, rel_fdr = fdr_est / p_sub, transform = 'rmRNAseq_BH', 
                         contrast = c_sub)
      }))
      return(fdr_sub)
    }))
    
  }))
  
  
  fdr_all = rbind(fdr_dream2, fdr_rmrnaseq, fdr_rmrnaseq_bh, fdr_vst_singular, fdr_gls)
  return(fdr_all)
}))

fdr_summary = all_fdr %>% 
  group_by(transform, contrast, level, n_per_group) %>% 
  summarize(med_fdr = median(fdr, na.rm = T), med_rel_fdr = median(rel_fdr, na.rm = T),
            mean_fdr = mean(fdr, na.rm = T), mean_rel_fdr = mean(rel_fdr, na.rm = T))
fdr_summary$level = factor(fdr_summary$level, levels = c( 0.01, 0.05, 0.10))

fdr_fig = fdr_summary %>% 
  # filter(transform %in% c("VST (WS)", "DREAM (WS)", "rmRNAseq")) %>% 
  ggplot() +
  aes(x = level, y = log2(mean_rel_fdr), color = transform, group = transform) +
  geom_point(size = 2) +
  # ylim(c(0, NA)) +
  facet_grid(rows = vars(contrast), cols = vars(n_per_group)) +
  geom_path() +
  scale_color_discrete(name = "Method") +
  ylab("log2 of relative FDR") +
  xlab("FDR level") +
  geom_hline(yintercept = 0)

print(fdr_fig)

####    Looking at statistical power    ####
all_pwr = bind_rows(pblapply(n_vec, function(n_sub){
  # n_sub = 3
  
  file_path2 = paste0(res_path, "just_dream_ROC_MIX_n_", n_sub, ".RDS")
  res_dream2 = readRDS(file_path2)$res_all_dream
  
  file_path3 = paste0(res_path, "just_rmrnaseq_ROC_MIX_n_", n_sub, ".RDS")
  res_rmrnaseq = readRDS(file_path3)$res_all
  
  file_path4 = paste0(res_path, "just_vst_w_singular_ROC_MIX_n_", n_sub, ".RDS")
  res_vst_singular = readRDS(file_path4)$res_all_vst
  
  file_path5 = paste0(res_path, "just_vst_w_gls_ROC_MIX_n_", n_sub, ".RDS")
  res_gls = readRDS(file_path5)$res_all_gls
  
  p_val_levels = c(0.10, 0.05, 0.01)
  
  sim_ids = unique(res_vst_singular$sim_num)
  
  
  pwr_vst_singular = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_vst_singular %>% filter(sim_num == sim_sub)
    pwr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      pwr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val_adj < p_sub, levels = c(F, T)))
        pwr_est = tab_sub[2, 2] / sum(tab_sub[2, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         pwr = pwr_est, transform = 'lmerSeq (WS)', 
                         contrast = c_sub)
      }))
      return(pwr_sub)
    }))
    
  }))
  
  pwr_gls = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_gls %>% filter(sim_num == sim_sub)
    res_sub_true = res_vst_singular %>% filter(sim_num == sim_sub)
    pwr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      res_sub_true2 = res_sub_true %>% filter(contrast == c_sub)
      pwr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub_true2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val_adj < p_sub, levels = c(F, T)))
        pwr_est = tab_sub[2, 2] / sum(tab_sub[2, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         pwr = pwr_est, transform = 'GLS', 
                         contrast = c_sub)
      }))
      return(pwr_sub)
    }))
    
  }))
  
  pwr_dream2 = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_dream2 %>% filter(sim_num == sim_sub)
    pwr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      pwr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$adj.P.Val < p_sub, levels = c(F, T)))
        pwr_est = tab_sub[2, 2] / sum(tab_sub[2, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         pwr = pwr_est, transform = 'DREAM (WS)', 
                         contrast = c_sub)
      }))
      return(pwr_sub)
    }))
    
  }))
  
  pwr_rmrnaseq = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_rmrnaseq %>% filter(sim_num == sim_sub)
    pwr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      pwr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$q_val < p_sub, levels = c(F, T)))
        pwr_est = tab_sub[2, 2] / sum(tab_sub[2, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         pwr = pwr_est, transform = 'rmRNAseq', 
                         contrast = c_sub)
      }))
      return(pwr_sub)
    }))
    
  }))
  
  pwr_rmrnaseq_bh = bind_rows(lapply(sim_ids, function(sim_sub){
    # sim_sub = 1
    res_sub = res_rmrnaseq %>% filter(sim_num == sim_sub) %>% mutate(p_val_bh = p.adjust(p_val_raw, method = "BH"))
    pwr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
      # c_sub = 'Interaction'
      res_sub2 = res_sub %>% filter(contrast == c_sub)
      pwr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
        # p_sub = .05
        tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                        est_de = factor(res_sub2$p_val_bh < p_sub, levels = c(F, T)))
        pwr_est = tab_sub[2, 2] / sum(tab_sub[2, ])
        ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                         pwr = pwr_est, transform = 'rmRNAseq_BH', 
                         contrast = c_sub)
      }))
      return(pwr_sub)
    }))
    
  }))
  
  
  pwr_all = rbind(pwr_dream2, pwr_rmrnaseq, pwr_rmrnaseq_bh, pwr_vst_singular, pwr_gls)
  return(pwr_all)
}))

pwr_summary = all_pwr %>% 
  group_by(transform, contrast, level, n_per_group) %>% 
  summarize(med_pwr = median(pwr, na.rm = T),
            mean_pwr = mean(pwr, na.rm = T))
pwr_summary$level = factor(pwr_summary$level, levels = c( 0.01, 0.05, 0.10))

df_plot_combo = fdr_summary %>% left_join(pwr_summary)
head(df_plot_combo)

df_plot_combo$contrast = factor(df_plot_combo$contrast, levels = c("Between", "Within", "Interaction"))

####   Generating power by FDR plots as seen in the paper and supplement    ####

###   0.05 FDR threshold
plt_combo_05 = df_plot_combo %>% 
  filter(transform %in% c("lmerSeq (WS)", "DREAM (WS)", "rmRNAseq_BH", "GLS")) %>%
  filter(level == "0.05") %>% 
  ggplot() +
  aes(x = log2(mean_rel_fdr), y = mean_pwr, color = transform, shape = as.factor(n_per_group)) +
  # geom_point(size = 3) +
  geom_path(aes(group = transform)) +
  geom_jitter(size = 3, width = .01, height = 0) +
  ylim(c(0, 1)) +
  theme_bw() +
  facet_grid(cols = vars(contrast)) +
  scale_shape_discrete(name = "N per Group") +
  scale_color_discrete(name = "Method", labels = c("DREAM", "lmerSeq (CS)", "lmerSeq (RI)", "rmRNAseq")) +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_vline(xintercept = log2(0.8), linetype='dotted') +
  ggtitle("Simulation 1: 0.05 FDR Level") +
  theme(plot.title = element_text(hjust = .5)) +
  # theme(text = element_text(size = 14)) +
  ylab("Sensitivity") +
  xlab("Log2 of Relative FDR (Observed / Nominal)") 

print(plt_combo_05)

###   0.01 FDR threshold
plt_combo_01 = df_plot_combo %>% 
  filter(transform %in% c("lmerSeq (WS)", "DREAM (WS)", "rmRNAseq_BH", "GLS")) %>%
  filter(level == "0.01") %>% 
  ggplot() +
  aes(x = log2(mean_rel_fdr), y = mean_pwr, color = transform, shape = as.factor(n_per_group)) +
  # geom_point(size = 3) +
  geom_path(aes(group = transform)) +
  geom_jitter(size = 3, width = .01, height = 0) +
  ylim(c(0, 1)) +
  theme_bw() +
  facet_grid(cols = vars(contrast)) +
  scale_shape_discrete(name = "N per Group") +
  scale_color_discrete(name = "Method", labels = c("DREAM", "lmerSeq (CS)", "lmerSeq (RI)", "rmRNAseq")) +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_vline(xintercept = log2(0.8), linetype='dotted') +
  ggtitle("Simulation 1: 0.01 FDR Level") +
  theme(plot.title = element_text(hjust = .5)) +
  # theme(text = element_text(size = 14)) +
  ylab("Sensitivity") +
  xlab("Log2 of Relative FDR (Observed / Nominal)") 

print(plt_combo_01)

###   0.10 FDR threshold
plt_combo_10 = df_plot_combo %>% 
  filter(transform %in% c("lmerSeq (WS)", "DREAM (WS)", "rmRNAseq_BH", "GLS")) %>%
  filter(level == "0.1") %>% 
  ggplot() +
  aes(x = log2(mean_rel_fdr), y = mean_pwr, color = transform, shape = as.factor(n_per_group)) +
  # geom_point(size = 3) +
  geom_path(aes(group = transform)) +
  geom_jitter(size = 3, width = .01, height = 0) +
  ylim(c(0, 1)) +
  theme_bw() +
  facet_grid(cols = vars(contrast)) +
  scale_shape_discrete(name = "N per Group") +
  scale_color_discrete(name = "Method", labels = c("DREAM", "lmerSeq (CS)", "lmerSeq (RI)", "rmRNAseq")) +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_vline(xintercept = log2(0.8), linetype='dotted') +
  ggtitle("Simulation 1: 0.10 FDR Level") +
  theme(plot.title = element_text(hjust = .5)) +
  # theme(text = element_text(size = 14)) +
  ylab("Sensitivity") +
  xlab("Log2 of Relative FDR (Observed / Nominal)") 

print(plt_combo_10)

####    Getting p-values for null features    ####
all_null_ps = bind_rows(pblapply(n_vec, function(n_sub){
  # n_sub = 3
  
  file_path2 = paste0(res_path, "just_dream_ROC_MIX_n_", n_sub, ".RDS")
  res_dream2 = readRDS(file_path2)$res_all_dream
  
  file_path3 = paste0(res_path, "just_rmrnaseq_ROC_MIX_n_", n_sub, ".RDS")
  res_rmrnaseq = readRDS(file_path3)$res_all
  
  file_path4 = paste0(res_path, "just_vst_w_singular_ROC_MIX_n_", n_sub, ".RDS")
  res_vst_singular = readRDS(file_path4)$res_all_vst
  
  file_path5 = paste0(res_path, "just_vst_w_gls_ROC_MIX_n_", n_sub, ".RDS")
  res_gls = readRDS(file_path5)$res_all_gls
  
  df_ps = rbind(res_dream2 %>% 
                  filter(beta_ints == 0) %>% 
                  dplyr::select(P.Value, n_per_grp, contrast) %>% 
                  mutate(Method = "DREAM"),
                res_rmrnaseq %>% 
                  filter(beta_ints == 0) %>% 
                  mutate(P.Value = p_val_raw) %>% 
                  dplyr::select(P.Value, n_per_grp, contrast) %>% 
                  mutate(Method = "rmRNAseq"),
                res_gls[res_vst_singular$beta_ints == 0, ] %>% 
                  # filter(beta_ints == 0) %>% 
                  mutate(P.Value = p_val) %>% 
                  dplyr::select(P.Value, n_per_grp, contrast) %>% 
                  mutate(Method = "lmerSeq (CS)"),
                res_vst_singular %>% 
                  filter(beta_ints == 0) %>% 
                  mutate(P.Value = p_val_raw) %>% 
                  dplyr::select(P.Value, n_per_grp, contrast) %>% 
                  mutate(Method = "lmerSeq (RI)"))
  return(df_ps)
}))

all_null_ps$contrast = factor(all_null_ps$contrast, levels = c("Between", "Within", "Interaction"))

####    Generating histograms of the null p-values   ####
plt_null_ps3 = all_null_ps %>% 
  filter(n_per_grp == 3) %>% 
  ggplot() +
  aes(x = P.Value, fill = Method) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(Method), cols = vars(contrast), scales = "free") +
  ggtitle("Simulation 1: N = 3 per Group") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 12),
        plot.title = element_text(hjust = .5)) +
  xlab("p-value") +
  ylab("Count")
print(plt_null_ps3)

plt_null_ps5 = all_null_ps %>% 
  filter(n_per_grp == 5) %>% 
  ggplot() +
  aes(x = P.Value, fill = Method) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(Method), cols = vars(contrast), scales = "free") +
  ggtitle("Simulation 1: N = 5 per Group") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 12),
        plot.title = element_text(hjust = .5)) +
  xlab("p-value") +
  ylab("Count")
print(plt_null_ps5)

plt_null_ps10 = all_null_ps %>% 
  filter(n_per_grp == 10) %>% 
  ggplot() +
  aes(x = P.Value, fill = Method) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(Method), cols = vars(contrast), scales = "free") +
  ggtitle("Simulation 1: N = 10 per Group") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 12),
        plot.title = element_text(hjust = .5)) +
  xlab("p-value") +
  ylab("Count")
print(plt_null_ps10)

plt_null_ps20 = all_null_ps %>% 
  filter(n_per_grp == 20) %>% 
  ggplot() +
  aes(x = P.Value, fill = Method) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(Method), cols = vars(contrast), scales = "free") +
  ggtitle("Simulation 1: N = 20 per Group") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 12),
        plot.title = element_text(hjust = .5)) +
  xlab("p-value") +
  ylab("Count")
print(plt_null_ps20)






