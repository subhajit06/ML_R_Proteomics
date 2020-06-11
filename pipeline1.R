
case_crtl_pair_t_test <-function(type){
  
  load("./IGM_IGG_med.Rdata")
  load('./sample_F_by_B_signals.RData')
  
  n_total_probes = length(sample_F_by_B_signals$row_names)
  n_samples =  length(sample_F_by_B_signals$col_names)
  
  #iso_ID_name_arr = map_isoform_ID_gene_name(IGM_IGG_med) 
  load("./iso_ID_name_arr.RData")
  
  #### IGM signal only 
  #### IGM_signals: sample in row and signals/features in column
  
  if(type == 'dup1'){
    dup1_signals = gen_max_signals(sample_F_by_B_signals,'dup1')
    IGM_signals = t(dup1_signals$mat_IGM_F_by_B)
    probe_names = dup1_signals$row_names
    sample_names = dup1_signals$col_names
    
  }
  if(type == 'dup2'){
    dup2_signals = gen_max_signals(sample_F_by_B_signals,'dup2')
    IGM_signals = t(dup2_signals$mat_IGM_F_by_B)
    probe_names = dup2_signals$row_names
    sample_names = dup2_signals$col_names
    
  }
  
  all_subjects_sample_cnt = as.vector(table(unlist(lapply(sample_names,function(x) 
    sprintf('%s',strsplit(x,'_')[[1]][1])))))
  
  subj_time_info = subject_and_time(sample_names)
  
  probe_set = 1:length(probe_names)
 
  #################################### 
  #### pre_post paired t-test analysis
  
  case_pre_post_samples_names = c("S1_m_12","S1_p_6",
                                  "S2_m_6","S2_p_5",
                                  "S3_m_11","S3_p_3")
  
  ctrl_pre_post_samples_names = c("S4_m_8","S4_p_6",
                                  "S5_m_11","S5_p_6",
                                  "S6_m_14","S6_p_4")
  
  both_pre_post_samples_names = c(case_pre_post_samples_names,ctrl_pre_post_samples_names)
  
  n_both = length(both_pre_post_samples_names)/2
  n_case = length(case_pre_post_samples_names)/2
  n_ctrl = length(ctrl_pre_post_samples_names)/2
  
  id_set_both = find_pre_post_case_ctrl_id(both_pre_post_samples_names,sample_names)
  id_set_case = find_pre_post_case_ctrl_id(case_pre_post_samples_names,sample_names)
  id_set_ctrl = find_pre_post_case_ctrl_id(ctrl_pre_post_samples_names,sample_names)
  
  all_sample_id_case = c(id_set_case$grp_pre_id,id_set_case$grp_post_id)
  all_sample_id_ctrl = c(id_set_ctrl$grp_pre_id,id_set_ctrl$grp_post_id)
  
  analysis_name_both = 'pre_vs_post'
  analysis_name_case = 'pre_vs_post_case'
  analysis_name_ctrl = 'pre_vs_post_ctrl'
  
  id_pre_post_thresholded_both = paired_t_test_pre_post(analysis_name_both,IGM_signals,iso_ID_name_arr,id_set_both,n_both,type='IGM',threshold = 0.05,test_type=0)
  id_pre_post_thresholded_case = paired_t_test_pre_post(analysis_name_case,IGM_signals,iso_ID_name_arr,id_set_case,n_case,type='IGM',threshold = 0.05,test_type=0)
  id_pre_post_thresholded_ctrl = paired_t_test_pre_post(analysis_name_ctrl,IGM_signals,iso_ID_name_arr,id_set_ctrl,n_ctrl,type='IGM',threshold = 0.05,test_type=0)
  
  #######################################
  ### case_control paired t_test analysis
  pre_case_ctrl_samples_names = c("S1_m_12","S2_m_6","S3_m_11",
                                  "S4_m_8","S5_m_11","S6_m_14")
                                  

  post_case_ctrl_samples_names = c("S1_p_6","S2_p_5","S3_p_3",
                                   "S4_p_6","S5_p_6","S6_p_4");
  
  both_case_ctrl_samples_names = c(pre_case_ctrl_samples_names,post_case_ctrl_samples_names)
  
  n_both = length(both_case_ctrl_samples_names)/2
  n_pre  = length(pre_case_ctrl_samples_names)/2
  n_post = length(post_case_ctrl_samples_names)/2
  
  id_set_both = find_pre_post_case_ctrl_id(both_pre_post_samples_names,sample_names)
  id_set_pre = find_pre_post_case_ctrl_id(pre_case_ctrl_samples_names,sample_names)
  id_set_post = find_pre_post_case_ctrl_id(post_case_ctrl_samples_names,sample_names)
  
  all_sample_id_pre = c(id_set_case$grp_case_id,id_set_case$grp_ctrl_id)
  all_sample_id_post = c(id_set_ctrl$grp_case_id,id_set_ctrl$grp_ctrl_id)
  
  analysis_name_both = 'case_vs_ctrl'
  analysis_name_pre  = 'case_vs_ctrl_pre'
  analysis_name_post = 'case_vs_ctrl_post'
  
  id_case_ctrl_thresholded_both = unpaired_t_test_pre_post(analysis_name_both,IGM_signals,iso_ID_name_arr,id_set_both,n_both,type='IGM',threshold = 0.05,test_type=0)
  id_case_ctrl_thresholded_pre  = unpaired_t_test_pre_post(analysis_name_pre,IGM_signals,iso_ID_name_arr,id_set_pre,n_pre,type='IGM',threshold = 0.05,test_type=0)
  id_case_ctrl_thresholded_post = unpaired_t_test_pre_post(analysis_name_post,IGM_signals,iso_ID_name_arr,id_set_post,n_post,type='IGM',threshold = 0.05,test_type=0)
  

  
  ##########################
  ##########################
  signal_and_isoform = NULL
  signal_and_isoform$id_pre_post_case = id_pre_post_thresholded_case
  signal_and_isoform$id_pre_post_ctrl = id_pre_post_thresholded_ctrl
  signal_and_isoform$id_pre_post_both = id_pre_post_thresholded_both
  
  signal_and_isoform$id_case_ctrl_pre = id_pre_post_thresholded_case
  signal_and_isoform$id_case_ctrl_post = id_pre_post_thresholded_ctrl
  signal_and_isoform$id_case_ctrl_both = id_pre_post_thresholded_both
  
  
  
  signal_and_isoform$signal = IGM_signals
  signal_and_isoform$iso_ID_name_arr = iso_ID_name_arr
  
  return(signal_and_isoform)
}

unpaired_t_test_case_ctrl <- function(analysis_name,IGM_signals_without_noise,
                                   iso_ID_name_arr,id_set,n,type,threshold,test_type){
  
  #### first n are in the group 1 and next n are in group 2 
  n = as.integer(n)
  
  grp_id_1 = c(1:n)
  grp_id_2 = c((n+1):(2*n))
  
  F_by_B_case_ctrl_n_paired_sample = t(IGM_signals_without_noise[c(id_set$grp_case_id,id_set$grp_ctrl_id),])
  
  if(test_type == 0){
    raw_p_val_n_sample_case_ctrl_unpaired_t_test = apply(
      as.data.frame(F_by_B_case_ctrl_n_paired_sample),1,ttest_unpaired, 
      grp1 = grp_id_1, grp2 = grp_id_2)
  }
  
  #### isoform_id name 
  id_case_ctrl_thresholded = which(raw_p_val_n_sample_case_ctrl_unpaired_t_test < threshold)
  
  return(id_case_ctrl_thresholded)
}


paired_t_test_pre_post <- function(analysis_name,IGM_signals_without_noise,
                                   iso_ID_name_arr,id_set,n,type,threshold,test_type){
  
  #### first n are in the group 1 and next n are in group 2 
  n = as.integer(n)
  
  grp_id_1 = c(1:n)
  grp_id_2 = c((n+1):(2*n))
  
  F_by_B_pre_post_n_paired_sample = t(IGM_signals_without_noise[c(id_set$grp_pre_id,id_set$grp_post_id),])
  
  if(test_type == 0){
    raw_p_val_n_sample_pre_post_paired_t_test = apply(
      as.data.frame(F_by_B_pre_post_n_paired_sample),1,ttest_paired, 
      grp1 = grp_id_1, grp2 = grp_id_2)
  }
  if(test_type == 1){
    raw_p_val_n_sample_pre_post_paired_t_test = apply(
      as.data.frame(F_by_B_pre_post_n_paired_sample),1,ttest_paired_greater, 
      grp1 = grp_id_1, grp2 = grp_id_2)
    
  }
  if(test_type == 2){
    raw_p_val_n_sample_pre_post_paired_t_test = apply(
      as.data.frame(F_by_B_pre_post_n_paired_sample),1,ttest_paired_less, 
      grp1 = grp_id_1,grp2 = grp_id_2)
  }
  
  #### isoform_id name 
  id_pre_post_thresholded = which(raw_p_val_n_sample_pre_post_paired_t_test < threshold)
  
  return(id_pre_post_thresholded)
}

remove_F_by_B_less_one <- function(IGM_signals,all_sample_id){
  
  set_rem = NULL
  n_sample = length(all_sample_id)
  
  for(i in 1:n_sample){
    set_rem_IGM = intersect(set_rem, which(IGM_signals[all_sample_id[i],] < 1))
  }
  
  return(set_rem)
}


gene_name_from_iso_id <-function(id_list,iso_ID_name_arr){
  
  gene_list = as.vector(unlist(lapply(id_list, function(x){ 
    iso_ID_name_arr$gene_name[which(iso_ID_name_arr$iso_id==x)]})))
  
  return(gene_list)
}

### type = 1 result with t-test with 0.05 significance
### type = 2 result with high-variable set of isoforms
main_run <-function(type){
  
  setwd("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4")
  
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/GPR_data_input.R")
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/glmnet_train.R")
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/use_psygenet.R")

  load_libs()
  
  signal_and_isoform_dup1 = case_crtl_pair_t_test(type='dup1')
  signal_and_isoform_dup2 = case_crtl_pair_t_test(type='dup2')  
  
  id_case_ctrl_intersection_dup1_dup2_union = intersect(signal_and_isoform_dup1$id_case_ctrl_union,
                                                        signal_and_isoform_dup2$id_case_ctrl_union)
  
  id_case_ctrl_intersection_dup1_dup2_intrsct = intersect(signal_and_isoform_dup1$id_case_ctrl_intrsct,
                                                          signal_and_isoform_dup2$id_case_ctrl_intrsct)
  
  id_both_intersection_dup1_dup2 = intersect(signal_and_isoform_dup1$id_both,
                                             signal_and_isoform_dup2$id_both)
  
  iso_ID_name_arr = signal_and_isoform_dup1$iso_ID_name_arr
  
  gn_case_ctrl_intersection_dup1_dup2_union = unique(iso_ID_name_arr$gene_name[id_case_ctrl_intersection_dup1_dup2_union])
  #gn_case_ctrl_intersection_dup1_dup2_intrsct = unique(iso_ID_name_arr$gene_name[id_case_ctrl_intersection_dup1_dup2_intrsct])
  gn_case_ctrl_intersection_dup1 = unique(iso_ID_name_arr$gene_name[signal_and_isoform_dup1$id_case_ctrl_intrsct])
  gn_case_ctrl_intersection_dup2 = unique(iso_ID_name_arr$gene_name[signal_and_isoform_dup2$id_case_ctrl_intrsct])
  gn_both_intersection_dup1_dup2 = unique(iso_ID_name_arr$gene_name[id_both_intersection_dup1_dup2])
  ###################
  X = signal_and_isoform_dup1$signal
  list_genes_dup1 = run_EN(X,iso_ID_name_arr)
  
  X = signal_and_isoform_dup2$signal
  list_genes_dup2 = run_EN(X,iso_ID_name_arr)
  
  library( psygenet2r )
  
  kk=2

  #for(kk in 1:length(list_genes)){
  genesOfInterest = unique(union(intersect(gn_case_ctrl_intersection_dup1_dup2_union,list_genes_dup1[[2]]),
                                 intersect(gn_case_ctrl_intersection_dup1_dup2_union,list_genes_dup2[[2]])))
  
  genesOfInterest = unique(union(intersect(gn_both_intersection_dup1_dup2,list_genes_dup1[[2]]),
                                 intersect(gn_both_intersection_dup1_dup2,list_genes_dup2[[2]])))
  
  m1 = psygenetGene(gene = genesOfInterest, database = "ALL", verbose  = FALSE, warnings = FALSE)
  fName = sprintf("./gene_disease_ratio_%0.2f_dup1_dup2.pdf",kk*0.1)
  pdf(fName,width=10,height=12)
  plot( m1 )
  dev.off()
  
  fName = sprintf("./gene_disease_hmap_ratio_%0.2f_dup1_dup2.pdf",kk*0.1)
  pdf(fName,width=11,height=9)
  plot( m1, type="GDA heatmap")
  dev.off()
  
  genesOfInterest = unique(list_genes_dup2[[kk]])
  m2 = psygenetGene(gene = genesOfInterest, database = "ALL", verbose  = FALSE, warnings = FALSE)
  fName = sprintf("./gene_disease_ratio_%0.2f_dup2.pdf",kk*0.1)
  pdf(fName,width=10,height=12)
  plot( m2 )
  dev.off()
  
  fName = sprintf("./gene_disease_hmap_ratio_%0.2f_dup2.pdf",kk*0.1)
  pdf(fName,width=11,height=9)
  plot( m2, type="GDA heatmap")
  dev.off()
  
}

