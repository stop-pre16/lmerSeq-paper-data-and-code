library(dplyr)
library(ggplot2)
library(reshape2)
library(pbapply)


n_vec = c(3, 5, 10, 20)
res_path = "/path/to/four_timepoint_sim/results/files/"

mtab = data.frame(method = c("lmerSeq_cat_ri", "lmerSeq_cont_ri", "lmerSeq_cont_ri_rs", "lmerSeq_cat_un",
                             "DREAM_cat_ri", "DREAM_cont_ri", "DREAM_cont_ri_rs",
                             "rmRNASeq_cat", "rmRNASeq_cont"),
                  file_prefix = c("lmerSeq_4tp_cat_ri_n_", "lmerSeq_4tp_cont_ri_n_", 
                                  "lmerSeq_4tp_cont_ri_rs_n_", "lmerSeq_4tp_cat_un_n_",
                                  "DREAM_4tp_cat_ri_n_", "DREAM_4tp_cont_ri_n_", 
                                  "DREAM_4tp_cont_ri_rs_n_",
                                  "just_rmrnaseq_4tp_cat_n_", "just_rmrnaseq_4tp_cont_n_"))
n_meth = nrow(mtab)

all_t1e = bind_rows(pblapply(n_vec, function(n_sub){
  # n_sub = 3
  
  res_list = lapply(1:n_meth, function(ii){
    # ii = 1
    ret = readRDS(paste0(res_path, mtab$file_prefix[ii], n_sub, ".RDS"))[[1]]
    return(ret)
  })
  # lapply(res_list, head)
  p_val_levels = c(0.05, 0.01, 0.001, 0.0001)
  
  sim_ids = 1:10
  
  t1e_all = bind_rows(lapply(1:n_meth, function(ii){
    res_meth = bind_rows(lapply(sim_ids, function(sim_sub){
      # sim_sub = 1
      res_sub = res_list[[ii]] %>% filter(sim_num == sim_sub)
      t1e_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
        # c_sub = 'Interaction'
        res_sub2 = res_sub %>% filter(contrast == c_sub)
        t1e_sub = bind_rows(lapply(p_val_levels, function(p_sub){
          # p_sub = .05
          tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                          est_de = factor(res_sub2$p_val_raw < p_sub, levels = c(F, T)))
          t1e_est = tab_sub[1, 2] / sum(tab_sub[1, ])
          ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                           t1e = t1e_est, rel_t1e = t1e_est / p_sub, method = mtab$method[ii], 
                           contrast = c_sub)
        }))
        return(t1e_sub)
      }))
      return(t1e_res)
    }))
    return(res_meth)
  }))
  return(t1e_all)
}))


t1e_summary = all_t1e %>% 
  group_by(method, contrast, level, n_per_group) %>% 
  summarize(med_t1e = median(t1e), med_rel_t1e = median(rel_t1e),
            mean_t1e = mean(t1e), mean_rel_t1e = mean(rel_t1e))
t1e_summary$level = factor(t1e_summary$level, levels = c(0.0001, 0.001, 0.01, 0.05))

ggplot(t1e_summary) +
  aes(x = level, y = log2(mean_rel_t1e), color = method, group = method) +
  geom_point(size = 2) +
  # ylim(c(0, NA)) +
  facet_grid(rows = vars(contrast), cols = vars(n_per_group)) +
  geom_path() +
  geom_hline(yintercept = 0)

all_fdr = bind_rows(pblapply(n_vec, function(n_sub){
  # n_sub = 3
  
  res_list = lapply(1:n_meth, function(ii){
    # ii = 1
    ret = readRDS(paste0(res_path, mtab$file_prefix[ii], n_sub, ".RDS"))[[1]]
    return(ret)
  })
  # lapply(res_list, head)
  # names(res_list[[8]])[10] = 'p_val_adj'
  # names(res_list[[9]])[10] = 'p_val_adj'
  res_list[[8]] = res_list[[8]] %>% mutate(p_val_adj = p.adjust(p_val_raw, method = "BH"))
  res_list[[9]] = res_list[[9]] %>% mutate(p_val_adj = p.adjust(p_val_raw, method = "BH"))
  p_val_levels = c(0.10, 0.05, 0.01)
  
  sim_ids = 1:10
  
  fdr_all = bind_rows(lapply(1:n_meth, function(ii){
    res_meth = bind_rows(lapply(sim_ids, function(sim_sub){
      # sim_sub = 1
      res_sub = res_list[[ii]] %>% filter(sim_num == sim_sub)
      fdr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
        # c_sub = 'Interaction'
        res_sub2 = res_sub %>% filter(contrast == c_sub)
        fdr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
          # p_sub = .05
          tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                          est_de = factor(res_sub2$p_val_adj < p_sub, levels = c(F, T)))
          fdr_est = tab_sub[1, 2] / sum(tab_sub[, 2])
          ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                           fdr = fdr_est, rel_fdr = fdr_est / p_sub, method = mtab$method[ii], 
                           contrast = c_sub)
        }))
        return(fdr_sub)
      }))
      return(fdr_res)
    }))
    return(res_meth)
  }))
  return(fdr_all)
}))


fdr_summary = all_fdr %>% 
  group_by(method, contrast, level, n_per_group) %>% 
  summarize(med_fdr = median(fdr, na.rm = T), med_rel_fdr = median(rel_fdr, na.rm = T),
            mean_fdr = mean(fdr, na.rm = T), mean_rel_fdr = mean(rel_fdr, na.rm = T))

fdr_summary$level = factor(fdr_summary$level, levels = c( 0.01, 0.05, 0.10))
ggplot(fdr_summary) +
  aes(x = level, y = log2(mean_rel_fdr), color = method, group = method) +
  geom_point(size = 2) +
  # ylim(c(0, NA)) +
  facet_grid(rows = vars(contrast), cols = vars(n_per_group)) +
  geom_path() +
  geom_hline(yintercept = 0)

all_pwr = bind_rows(pblapply(n_vec, function(n_sub){
  # n_sub = 3
  
  res_list = lapply(1:n_meth, function(ii){
    # ii = 1
    ret = readRDS(paste0(res_path, mtab$file_prefix[ii], n_sub, ".RDS"))[[1]]
    return(ret)
  })
  # lapply(res_list, head)
  # names(res_list[[8]])[10] = 'p_val_adj'
  # names(res_list[[9]])[10] = 'p_val_adj'
  res_list[[8]] = res_list[[8]] %>% mutate(p_val_adj = p.adjust(p_val_raw, method = "BH"))
  res_list[[9]] = res_list[[9]] %>% mutate(p_val_adj = p.adjust(p_val_raw, method = "BH"))
  p_val_levels = c(0.10, 0.05, 0.01)
  
  sim_ids = 1:10
  
  pwr_all = bind_rows(lapply(1:n_meth, function(ii){
    res_meth = bind_rows(lapply(sim_ids, function(sim_sub){
      # sim_sub = 1
      res_sub = res_list[[ii]] %>% filter(sim_num == sim_sub)
      pwr_res = bind_rows(lapply(unique(res_sub$contrast), function(c_sub){
        # c_sub = 'Interaction'
        res_sub2 = res_sub %>% filter(contrast == c_sub)
        pwr_sub = bind_rows(lapply(p_val_levels, function(p_sub){
          # p_sub = .05
          tab_sub = table(true_de = factor(res_sub2$true_de, levels = c(F, T)),
                          est_de = factor(res_sub2$p_val_adj < p_sub, levels = c(F, T)))
          pwr_est = tab_sub[2, 2] / sum(tab_sub[2, ])
          ret = data.frame(n_per_group = n_sub, sim_num = sim_sub, level = p_sub, 
                           pwr = pwr_est, method = mtab$method[ii], 
                           contrast = c_sub)
        }))
        return(pwr_sub)
      }))
      return(pwr_res)
    }))
    return(res_meth)
  }))
  return(pwr_all)
}))

pwr_summary = all_pwr %>% 
  group_by(method, contrast, level, n_per_group) %>% 
  summarize(med_pwr = median(pwr, na.rm = T),
            mean_pwr = mean(pwr, na.rm = T))
pwr_summary$level = factor(pwr_summary$level, levels = c( 0.01, 0.05, 0.10))

df_plot_combo = fdr_summary %>% left_join(pwr_summary)
head(df_plot_combo)
df_plot_combo = df_plot_combo %>% 
  rowwise() %>% 
  mutate(ovr_method = strsplit(method, split = "_")[[1]][1]) %>% 
  ungroup() %>% 
  data.frame
df_plot_combo$mean_rel_fdr[df_plot_combo$mean_rel_fdr == 0 | is.na(df_plot_combo$mean_rel_fdr)] = min(df_plot_combo$mean_rel_fdr[!(df_plot_combo$mean_rel_fdr == 0 | is.na(df_plot_combo$mean_rel_fdr))])

summary(df_plot_combo)

mtab$display = c("lmerSeq (Cat, RI)",
                 "lmerSeq (Cont, RI)",
                 "lmerSeq (Cont, RI + RS)",
                 "lmerSeq (Cat, Un)",
                 "DREAM (Cat, RI)",
                 "DREAM (Cont, RI)",
                 "DREAM (Cont, RI + RS)",
                 "rmRNAseq (Cat, CAR)",
                 "rmRNAseq (Cont, CAR)")

df_plot_combo$contrast = factor(df_plot_combo$contrast, levels = c("Between", "Within", "Interaction"))
plt_combo5 = df_plot_combo %>% 
  rowwise() %>% 
  mutate(mean_rel_fdr2 = max(mean_rel_fdr, 2^-4)) %>%
  filter(contrast != "Overall") %>% 
  # filter(transform %in% c("lmerSeq (WS)", "DREAM (WS)", "rmRNAseq", "GLS")) %>%
  filter(level == "0.05") %>% 
  left_join(mtab) %>% 
  filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  ggplot() +
  # aes(x = log2(med_rel_fdr), y = med_pwr, color = method, shape = as.factor(n_per_group)) +
  aes(x = log2(mean_rel_fdr2), y = mean_pwr, color = display, shape = as.factor(n_per_group)) +
  # aes(x = log2(mean_rel_fdr), y = mean_pwr, color = ovr_method, shape = as.factor(n_per_group)) +
  # geom_point(size = 3) +
  geom_path(aes(group = ovr_method)) +
  geom_jitter(size = 3, width = 0, height = 0) +
  ylim(c(0, 1)) +
  # facet_grid(rows = vars(contrast), cols = vars(ovr_method)) +
  facet_grid(cols = vars(contrast)) +
  scale_shape_discrete(name = "N per Group") +
  scale_color_discrete(name = "Method") +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_vline(xintercept = log2(0.8), linetype='dotted') +
  ylab("Sensitivity") +
  theme_bw() +
  ggtitle("Simulation 2: \"Correct\" Model Specification at 0.05 FDR Level") +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = .5)) +
  xlab("Log2 of Relative FDR (Observed / Nominal)") 
print(plt_combo5)

plt_combo6 = df_plot_combo %>% 
  rowwise() %>% 
  mutate(mean_rel_fdr2 = max(mean_rel_fdr, 2^-4)) %>% 
  filter(contrast != "Overall") %>% 
  # filter(transform %in% c("lmerSeq (WS)", "DREAM (WS)", "rmRNAseq", "GLS")) %>%
  filter(level == "0.05") %>% 
  left_join(mtab) %>% 
  filter(!(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont"))) %>% 
  ggplot() +
  # aes(x = log2(med_rel_fdr), y = med_pwr, color = method, shape = as.factor(n_per_group)) +
  aes(x = log2(mean_rel_fdr2), y = mean_pwr, color = display, shape = as.factor(n_per_group)) +
  # geom_point(size = 3) +
  geom_path(aes(group = method)) +
  geom_jitter(size = 3, width = 0, height = 0) +
  ylim(c(0, 1)) +
  # facet_grid(rows = vars(contrast), cols = vars(ovr_method)) +
  facet_grid(cols = vars(contrast)) +
  theme_bw() +
  ggtitle("Simulation 2: Models with Misspecification at 0.05 FDR Level") +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = .5)) +
  scale_shape_discrete(name = "N per Group") +
  scale_color_discrete(name = "Method") +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_vline(xintercept = log2(0.8), linetype='dotted') +
  ylab("Sensitivity") +
  xlab("Log2 of Relative FDR (Observed / Nominal)") 
print(plt_combo6)



plt_combo7 = df_plot_combo %>% 
  rowwise() %>% 
  mutate(mean_rel_fdr2 = max(mean_rel_fdr, 2^-4)) %>% 
  filter(contrast != "Overall") %>% 
  # filter(transform %in% c("lmerSeq (WS)", "DREAM (WS)", "rmRNAseq", "GLS")) %>%
  filter(level == "0.05") %>% 
  left_join(mtab) %>% 
  # filter(!(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont"))) %>% 
  ggplot() +
  # aes(x = log2(med_rel_fdr), y = med_pwr, color = method, shape = as.factor(n_per_group)) +
  aes(x = log2(mean_rel_fdr2), y = mean_pwr, color = display, shape = as.factor(n_per_group)) +
  # geom_point(size = 3) +
  geom_path(aes(group = method)) +
  geom_jitter(size = 3, width = 0, height = 0) +
  ylim(c(0, 1)) +
  # facet_grid(rows = vars(contrast), cols = vars(ovr_method)) +
  facet_grid(cols = vars(contrast), rows = vars(ovr_method)) +
  theme_bw(base_size = 12) +
  ggtitle("Simulation 2: All Models at 0.05 FDR Level") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_shape_discrete(name = "N per Group") +
  scale_color_discrete(name = "Method") +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_vline(xintercept = log2(0.8), linetype='dotted') +
  ylab("Sensitivity") +
  xlab("Log2 of Relative FDR (Observed / Nominal)") 
print(plt_combo7)


plt_combo8 = df_plot_combo %>% 
  rowwise() %>% 
  mutate(mean_rel_fdr2 = max(mean_rel_fdr, 2^-4)) %>% 
  filter(contrast != "Overall") %>% 
  # filter(transform %in% c("lmerSeq (WS)", "DREAM (WS)", "rmRNAseq", "GLS")) %>%
  filter(level == "0.1") %>% 
  left_join(mtab) %>% 
  # filter(!(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont"))) %>% 
  ggplot() +
  # aes(x = log2(med_rel_fdr), y = med_pwr, color = method, shape = as.factor(n_per_group)) +
  aes(x = log2(mean_rel_fdr2), y = mean_pwr, color = display, shape = as.factor(n_per_group)) +
  # geom_point(size = 3) +
  geom_path(aes(group = method)) +
  geom_jitter(size = 3, width = 0, height = 0) +
  ylim(c(0, 1)) +
  # facet_grid(rows = vars(contrast), cols = vars(ovr_method)) +
  facet_grid(cols = vars(contrast), rows = vars(ovr_method)) +
  theme_bw(base_size = 12) +
  ggtitle("Simulation 2: All Models at 0.10 FDR Level") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_shape_discrete(name = "N per Group") +
  scale_color_discrete(name = "Method") +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_vline(xintercept = log2(0.8), linetype='dotted') +
  ylab("Sensitivity") +
  xlab("Log2 of Relative FDR (Observed / Nominal)") 
print(plt_combo8)


plt_combo9 = df_plot_combo %>% 
  rowwise() %>% 
  mutate(mean_rel_fdr2 = max(mean_rel_fdr, 2^-4)) %>% 
  filter(contrast != "Overall") %>% 
  # filter(transform %in% c("lmerSeq (WS)", "DREAM (WS)", "rmRNAseq", "GLS")) %>%
  filter(level == "0.01") %>% 
  left_join(mtab) %>% 
  # filter(!(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont"))) %>% 
  ggplot() +
  # aes(x = log2(med_rel_fdr), y = med_pwr, color = method, shape = as.factor(n_per_group)) +
  aes(x = log2(mean_rel_fdr2), y = mean_pwr, color = display, shape = as.factor(n_per_group)) +
  # geom_point(size = 3) +
  geom_path(aes(group = method)) +
  geom_jitter(size = 3, width = 0, height = 0) +
  ylim(c(0, 1)) +
  # facet_grid(rows = vars(contrast), cols = vars(ovr_method)) +
  facet_grid(cols = vars(contrast), rows = vars(ovr_method)) +
  theme_bw(base_size = 12) +
  ggtitle("Simulation 2: All Models at 0.01 FDR Level") +
  theme(plot.title = element_text(hjust = .5)) +
  scale_shape_discrete(name = "N per Group") +
  scale_color_discrete(name = "Method") +
  geom_vline(xintercept = 0, linetype='dashed') +
  geom_vline(xintercept = log2(0.8), linetype='dotted') +
  ylab("Sensitivity") +
  xlab("Log2 of Relative FDR (Observed / Nominal)") 
print(plt_combo9)


all_null_ps = bind_rows(pblapply(n_vec, function(n_sub){
  # n_sub = 3
  
  res_list = bind_rows(lapply(1:n_meth, function(ii){
    # ii = 1
    ret = readRDS(paste0(res_path, mtab$file_prefix[ii], n_sub, ".RDS"))[[1]]
    ret = ret %>% filter(beta_ints == 0) %>% dplyr::select(p_val_raw, contrast) %>% mutate(method = mtab$method[ii])
    return(ret)
  }))
  # lapply(res_list, head)
  res_list = res_list %>% mutate(n_per_group = n_sub)
  return(res_list)
}))

###  Next line will take a while to run unfortunately
all_null_ps = all_null_ps %>% 
  rowwise() %>% 
  mutate(ovr_method = strsplit(method, split = "_")[[1]][1]) %>% 
  ungroup() %>% 
  data.frame

all_null_ps$contrast = factor(all_null_ps$contrast, levels = c("Between", "Within", "Interaction"))

plt_null_ps3 = all_null_ps %>% 
  filter(n_per_group == 3) %>% 
  filter(contrast != "Overall") %>% 
  filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  left_join(mtab) %>% 
  ggplot() +
  aes(x = p_val_raw, fill = ovr_method) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(display), cols = vars(contrast), scales = "free") +
  ggtitle("Simulation 2: \"Correct\" Model Specification at N = 3 per Group") +
  ylab("Count") +
  xlab("p-value") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 12),
        plot.title = element_text(hjust = .5))
print(plt_null_ps3)

plt_null_ps5 = all_null_ps %>% 
  filter(n_per_group == 5) %>% 
  filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  filter(contrast != "Overall") %>% 
  left_join(mtab) %>% 
  ggplot() +
  aes(x = p_val_raw, fill = ovr_method) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(display), cols = vars(contrast), scales = "free") +
  theme_bw() +
  ggtitle("Simulation 2: \"Correct\" Model Specification at N = 5 per Group") +
  ylab("Count") +
  xlab("p-value") +
  theme(legend.position = "none",
        text = element_text(size = 12),
        plot.title = element_text(hjust = .5))
print(plt_null_ps5)

plt_null_ps5_all = all_null_ps %>% 
  filter(n_per_group == 5) %>% 
  # filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  filter(contrast != "Overall") %>% 
  left_join(mtab) %>% 
  ggplot() +
  aes(x = p_val_raw, fill = display) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(display), cols = vars(contrast), scales = "free") +
  theme_bw() +
  ggtitle("Simulation 2: All Models at N = 5 per Group") +
  ylab("Count") +
  xlab("p-value") +
  theme(legend.position = "none",
        text = element_text(size = 7),
        plot.title = element_text(hjust = .5))
print(plt_null_ps5_all)

plt_null_ps3_all = all_null_ps %>% 
  filter(n_per_group == 3) %>% 
  # filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  filter(contrast != "Overall") %>% 
  left_join(mtab) %>% 
  ggplot() +
  aes(x = p_val_raw, fill = display) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(display), cols = vars(contrast), scales = "free") +
  theme_bw() +
  ggtitle("Simulation 2: All Models at N = 3 per Group") +
  ylab("Count") +
  xlab("p-value") +
  theme(legend.position = "none",
        text = element_text(size = 7),
        plot.title = element_text(hjust = .5))
print(plt_null_ps3_all)

plt_null_ps10_all = all_null_ps %>% 
  filter(n_per_group == 10) %>% 
  # filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  filter(contrast != "Overall") %>% 
  left_join(mtab) %>% 
  ggplot() +
  aes(x = p_val_raw, fill = display) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(display), cols = vars(contrast), scales = "free") +
  theme_bw() +
  ggtitle("Simulation 2: All Models at N = 10 per Group") +
  ylab("Count") +
  xlab("p-value") +
  theme(legend.position = "none",
        text = element_text(size = 7),
        plot.title = element_text(hjust = .5))
print(plt_null_ps10_all)

plt_null_ps20_all = all_null_ps %>% 
  filter(n_per_group == 20) %>% 
  # filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  filter(contrast != "Overall") %>% 
  left_join(mtab) %>% 
  ggplot() +
  aes(x = p_val_raw, fill = display) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(display), cols = vars(contrast), scales = "free") +
  theme_bw() +
  ggtitle("Simulation 2: All Models at N = 20 per Group") +
  ylab("Count") +
  xlab("p-value") +
  theme(legend.position = "none",
        text = element_text(size = 7),
        plot.title = element_text(hjust = .5))
print(plt_null_ps20_all)

plt_null_ps10 = all_null_ps %>% 
  filter(n_per_group == 10) %>% 
  filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  filter(contrast != "Overall") %>% 
  left_join(mtab) %>% 
  ggplot() +
  aes(x = p_val_raw, fill = ovr_method) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(display), cols = vars(contrast), scales = "free") +
  theme_bw() +
  ggtitle("Simulation 2: \"Correct\" Model Specification at N = 10 per Group") +
  ylab("Count") +
  xlab("p-value") +
  theme(legend.position = "none",
        text = element_text(size = 12),
        plot.title = element_text(hjust = .5))
print(plt_null_ps10)

plt_null_ps20 = all_null_ps %>% 
  filter(n_per_group == 20) %>% 
  filter(method %in% c("lmerSeq_cont_ri_rs", "DREAM_cont_ri_rs", "rmRNASeq_cont")) %>% 
  filter(contrast != "Overall") %>% 
  ggplot() +
  aes(x = p_val_raw, fill = ovr_method) +
  geom_histogram(alpha = .7, breaks = seq(0, 1, by = .02)) +
  facet_grid(rows = vars(ovr_method), cols = vars(contrast), scales = "free") +
  ylab("Count") +
  xlab("p-value") +
  ggtitle("N = 20 per Group") +
  theme(legend.position = "none",
        text = element_text(size = 14),
        plot.title = element_text(hjust = .5))
print(plt_null_ps20)
