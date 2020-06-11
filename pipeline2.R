case_crtl_up_down_pair_t_test <-function(type){
  
  load("./IGM_IGG_med.Rdata")
  load('./sample_F_by_B_signals.RData')
  
  n_total_probes = length(sample_F_by_B_signals$row_names)
  n_samples =  length(sample_F_by_B_signals$col_names)
  
  #iso_ID_name_arr = map_isoform_ID_gene_name(IGM_IGG_med) 
  load("./iso_ID_name_arr.RData")
  
  #### We use IGM signal only 
  #### IGM_signals: sample in row and signals/features in column
  if(type == 'avg'){
    avg_signals = gen_max_signals(sample_F_by_B_signals,'dup1')
    IGM_signals = t(avg_signals$mat_IGM_F_by_B)
    probe_names = avg_signals$row_names
    sample_names = avg_signals$col_names
  }
  
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
  
  case_pre_post_samples_names = c("S1_m_6","S1_p_6",
                                  "S2_m_6","S2_p_5",
                                  "S3_m_5","S3_p_3")
  
  ctrl_pre_post_samples_names = c("S4_m_8","S4_p_6",
                                  "S5_m_11","S5_p_6",
                                  "S6_m_14","S6_p_4")
  
  n_case = length(case_pre_post_samples_names)/2
  n_ctrl = length(ctrl_pre_post_samples_names)/2
  
  id_set_case = find_pre_post_case_ctrl_id(case_pre_post_samples_names,sample_names)
  id_set_ctrl = find_pre_post_case_ctrl_id(ctrl_pre_post_samples_names,sample_names)
  
  analysis_name_case = 'pre_vs_post_case'
  analysis_name_ctrl = 'pre_vs_post_ctrl'
  
  ######## case
  
  down_reg_case_pre_greater_post = paired_t_test_pre_post(analysis_name_case,IGM_signals,iso_ID_name_arr,id_set_case,n_case,type='IGM',threshold = 0.025,test_type=1)
  up_reg_case_pre_less_post = paired_t_test_pre_post(analysis_name_case,IGM_signals,iso_ID_name_arr,id_set_case,n_case,type='IGM',threshold = 0.025,test_type=2)
  
  flat_reg_case_pre_equiv_post = setdiff(probe_set,union(down_reg_case_pre_greater_post,up_reg_case_pre_less_post))
  
  ######## control
  
  down_reg_ctrl_pre_greater_post = paired_t_test_pre_post(analysis_name_ctrl,IGM_signals,iso_ID_name_arr,id_set_ctrl,n_ctrl,type='IGM',threshold = 0.025,test_type=1)
  up_reg_ctrl_pre_less_post = paired_t_test_pre_post(analysis_name_ctrl,IGM_signals,iso_ID_name_arr,id_set_ctrl,n_ctrl,type='IGM',threshold = 0.025,test_type=2)
  
  flat_reg_ctrl_pre_equiv_post = setdiff(probe_set,union(down_reg_ctrl_pre_greater_post,up_reg_ctrl_pre_less_post))

  ##########################
  ##########################
  signal_and_isoform = NULL
  signal_and_isoform$id_D_case = down_reg_case_pre_greater_post
  signal_and_isoform$id_U_case = up_reg_case_pre_less_post
  signal_and_isoform$id_F_case = flat_reg_case_pre_equiv_post
  
  signal_and_isoform$id_D_ctrl = down_reg_ctrl_pre_greater_post
  signal_and_isoform$id_U_ctrl = up_reg_ctrl_pre_less_post
  signal_and_isoform$id_F_ctrl = flat_reg_ctrl_pre_equiv_post
  
  signal_and_isoform$signal = IGM_signals
  signal_and_isoform$iso_ID_name_arr = iso_ID_name_arr
  
  return(signal_and_isoform)
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
  raw_p_val_n_sample_pre_post_paired_t_test_threshold_sorted = 
    sort(raw_p_val_n_sample_pre_post_paired_t_test[id_pre_post_thresholded],index.return=T)
  
  id_pre_post_thresholded_sorted = id_pre_post_thresholded[raw_p_val_n_sample_pre_post_paired_t_test_threshold_sorted$ix]
    
  return(id_pre_post_thresholded_sorted)
}

gene_name_from_iso_id <-function(id_list,iso_ID_name_arr){
  
  gene_list = as.vector(unlist(lapply(id_list, function(x){ 
    iso_ID_name_arr$gene_name[which(iso_ID_name_arr$iso_id==x)]})))
  
  return(gene_list)
}

main_run <-function(){
  
  setwd("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4")
  
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/GPR_data_input.R")
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/glmnet_GA_train.R")
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/use_psygenet.R")
  
  load_libs()
  
  #signal_and_isoform_dup1 = case_crtl_up_down_pair_t_test(type='dup1')
  #signal_and_isoform_dup2 = case_crtl_up_down_pair_t_test(type='dup2')  

  signal_and_isoform_avg = case_crtl_up_down_pair_t_test(type='avg')
  
  ##########
  
  iso_ID_name_arr = signal_and_isoform_avg$iso_ID_name_arr
  ##############
  gn_D_case_avg = unique(iso_ID_name_arr$gene_name[signal_and_isoform_avg$id_D_case])
  gn_U_case_avg = unique(iso_ID_name_arr$gene_name[signal_and_isoform_avg$id_U_case])
  gn_F_case_avg = unique(iso_ID_name_arr$gene_name[signal_and_isoform_avg$id_F_case])
  
  gn_U_ctrl_avg = unique(iso_ID_name_arr$gene_name[signal_and_isoform_avg$id_U_ctrl])
  gn_D_ctrl_avg = unique(iso_ID_name_arr$gene_name[signal_and_isoform_avg$id_D_ctrl])
  gn_F_ctrl_avg = unique(iso_ID_name_arr$gene_name[signal_and_isoform_avg$id_F_ctrl])
  
  #######################
  
  #### interesting cases
  gn_D_case_U_ctrl_avg = intersect(gn_D_case_avg,gn_U_ctrl_avg)
  gn_D_case_F_ctrl_avg = intersect(gn_D_case_avg,gn_F_ctrl_avg)
  gn_D_case_D_ctrl_avg = intersect(gn_D_case_avg,gn_D_ctrl_avg)
  
  gn_U_case_D_ctrl_avg = intersect(gn_U_case_avg,gn_D_ctrl_avg)
  gn_U_case_F_ctrl_avg = intersect(gn_U_case_avg,gn_F_ctrl_avg)
  gn_U_case_U_ctrl_avg = intersect(gn_U_case_avg,gn_U_ctrl_avg)
  
  gn_F_case_D_ctrl_avg = intersect(gn_F_case_avg,gn_D_ctrl_avg)
  gn_F_case_U_ctrl_avg = intersect(gn_F_case_avg,gn_U_ctrl_avg)
  gn_F_case_F_ctrl_avg = intersect(gn_F_case_avg,gn_F_ctrl_avg)
  
  ##############################
  gn_all_9_types = list()
  gn_all_9_types[[1]] = gn_D_case_U_ctrl_avg
  gn_all_9_types[[2]] = gn_D_case_F_ctrl_avg
  #gn_all_9_types[[3]] = gn_D_case_D_ctrl_avg
  gn_all_9_types[[3]] = gn_U_case_D_ctrl_avg
  gn_all_9_types[[4]] = gn_U_case_F_ctrl_avg
  #gn_all_9_types[[6]] = gn_U_case_U_ctrl_avg
  gn_all_9_types[[5]] = gn_F_case_D_ctrl_avg
  gn_all_9_types[[6]] = gn_F_case_U_ctrl_avg
  #gn_all_9_types[[9]] = gn_F_case_F_ctrl_avg
  
  t_test_gn_list = Reduce(union,gn_all_9_types)
  
  ################################
  #### run EN model ####
  
  X = signal_and_isoform_avg$signal
  list_genes_all_avg_1 = run_EN(X,iso_ID_name_arr)
  list_genes_all_avg_2 = run_EN_1(X,iso_ID_name_arr,type=1)
  list_genes_all_avg_3 = run_EN_1(X,iso_ID_name_arr,type=2)
  
  
  ##############

  #unlist(lapply(list_genes_all_avg_1,function(x){length(intersect(list_genes_all_avg_2[[2]],x))}))
  unlist(lapply(gn_all_9_types,function(x){length(intersect(list_genes_all_avg_2[[2]],x))}))
  unlist(lapply(gn_all_9_types,function(x){length(intersect(list_genes_all_avg_1[[2]],x))}))
  
  ##############
  sample_names = rownames(X)
  
  V1 = as.vector(unlist(lapply(sample_names, function(x) {ifelse(strsplit(x,"_")[[1]][2] == 'm',-1,1)}))) 
  V2 = as.vector(as.numeric(unlist(lapply(sample_names, function(x) {strsplit(x,"_")[[1]][3]}))))
  GA = V1*V2

  Z = c(rep(1,16),rep(0,9))
  Y = GA 
  rf=randomForest(X, Y, proximity=TRUE, importance=FALSE, norm.votes=TRUE) 
  id_set_rf = data.frame(importance(rf))
  rf_sel_features = row.names(id_set_rf)[id_set_rf$IncNodePurity>0]
  rf_genes = unique(unlist(lapply(rf_sel_features,function(x)
    {iso_ID_name_arr$gene_name[iso_ID_name_arr$iso_id == x]})))
  
  glist1 = union(intersect(rf_genes,gn_D_case_U_ctrl_avg),intersect(rf_genes,gn_U_case_D_ctrl_avg))
  genesOfInterest = glist1
  #list_genes_case_avg = run_EN(X,iso_ID_name_arr,id_sel=c(1:16))
  #list_genes_ctrl_avg = run_EN(X,iso_ID_name_arr,id_sel=c(17:25))
  
  ################################
  
  ################################
  library( psygenet2r )
  
  kk=2
  
  #for(kk in 1:length(list_genes)){
  #genesOfInterest = unique(union(intersect(gn_case_ctrl_intersection_dup1_dup2_union,list_genes_dup1[[2]]),
  #                               intersect(gn_case_ctrl_intersection_dup1_dup2_union,list_genes_dup2[[2]])))
  
  #genesOfInterest = unique(union(intersect(gn_both_intersection_dup1_dup2,list_genes_dup1[[2]]),
  #                               intersect(gn_both_intersection_dup1_dup2,list_genes_dup2[[2]])))
  
  genesOfInterest = intersect(list_genes_all_avg_2[[2]],t_test_gn_list)
  
  m1 = psygenetGene(gene = genesOfInterest, database = "ALL", verbose  = FALSE, warnings = FALSE)
  fName = sprintf("./gene_disease_ratio_avg_2.pdf")
  pdf(fName,width=10,height=12)
  plot( m1 )
  dev.off()
  
  #fName = sprintf("./gene_disease_hmap_ratio_avg.pdf")
  #pdf(fName,width=11,height=9)
  #plot( m1, type="GDA heatmap")
  #dev.off()
  
  ### GWAS catalog search
  gwas_catalog_tbl = read.table('./gwas_catalog_v1.0-associations_e93_r2018-12-21.tsv',sep='\t',header=TRUE,fill=TRUE)
  
  id_list = vector()
  for(i in 1:length(genesOfInterest)){
    ret = which(grepl(genesOfInterest[i],gwas_catalog_tbl$REPORTED.GENE.S.) == TRUE)
    id_list = c(id_list,ret)
  }
  id_list = unique(id_list)
  
  d_list = list()
  for(i in 1:length(id_list)){
    id_list[i]
    d_list[[i]] = gwas_catalog_tbl$DISEASE.TRAIT[i] 
  }
  write.table(d_list,)
  
}

