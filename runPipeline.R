
find_most_changing_isoform_id <-function(type = 'dup1', maxNum = 5){
  
  load("./IGM_IGG_med.Rdata")
  load('./sample_F_by_B_signals.RData')
  
  n_total_probes = length(sample_F_by_B_signals$row_names)
  n_samples =  length(sample_F_by_B_signals$col_names)
  
  iso_ID_name_arr = map_isoform_ID_gene_name(IGM_IGG_med) 
  
  #### IGM signal only 
  
  dup1_signals = gen_max_signals(sample_F_by_B_signals,'dup1')
  dup2_signals = gen_max_signals(sample_F_by_B_signals,'dup2')
  
  probe_names = dup1_signals$row_names
  sample_names = dup1_signals$col_names
  
  all_subjects_sample_cnt = as.vector(table(unlist(lapply(sample_names,function(x) 
    sprintf('%s',strsplit(x,'_')[[1]][1])))))
                                                                 
  subj_time_info = subject_and_time(sample_names)
  
  if(type == 'dup1'){
    IGM_signals = t(dup1_signals$mat_IGM_F_by_B)
  }
  if(type == 'dup2'){
    IGM_signals = t(dup2_signals$mat_IGM_F_by_B)
  }
  
  pre_post_samples_names = c("S1_m_12","S1_p_6",
                             "S2_m_6","S2_p_5",
                             "S3_m_11","S3_p_3",
                             "S4_m_8","S4_p_6",
                             "S5_m_11","S5_p_6",
                             "S6_m_14","S6_p_4")
  
  n = length(pre_post_samples_names)/2
  
  id_set = find_pre_post_id(pre_post_samples_names,sample_names)
  all_sample_id = c(id_set$grp_pre_id,id_set$grp_post_id)
  
  #### IGM_signals: sample in row and signals/features in column
  
  set_rem = remove_F_by_B_less_one(IGM_signals,all_sample_id)
  if(length(set_rem) > 0){  
    IGM_signals_without_noise = IGM_signals[,-set_rem]
  }else{
    IGM_signals_without_noise =  IGM_signals
  }
  ##################################
  ##################################
  ### pre and post changing signal by values
  analysis_name = 'pre_vs_post'
  pre_post_samples_names = c("S1_m_12","S1_p_6",
                             "S2_m_6","S2_p_5",
                             "S3_m_11","S3_p_3",
                             "S4_m_8","S4_p_6",
                             "S5_m_11","S5_p_6",
                             "S6_m_14","S6_p_4")
  
  n = length(pre_post_samples_names)/2
  
  id_set = find_pre_post_id(pre_post_samples_names,sample_names)
  
  pre_post_signals = NULL
  pre_post_signals = matrix(rep(0,length(probe_names)*n),c(n,length(probe_names)))
  
  max_change = NULL
  
  v1 = v2 = v3 = v4 = NULL
  
  max_sig_pre_post_change_id_set = NULL
  
  for (i in 1:n){
    pre_post_signals[i,] = IGM_signals_without_noise[id_set$grp_pre_id[i],]-IGM_signals_without_noise[id_set$grp_post_id[i],]
    
    val_id_1 = sort(pre_post_signals[i,],index.return=T,decreasing = F) ### + to -
    val_id_2 = sort(pre_post_signals[i,],index.return=T,decreasing = T) ### - to +
    
    l1_1 = IGM_signals_without_noise[id_set$grp_pre_id[i],val_id_1$ix[1:maxNum]]
    l2_1 = IGM_signals_without_noise[id_set$grp_post_id[i],val_id_1$ix[1:maxNum]]
    
    l1_2 = IGM_signals_without_noise[id_set$grp_pre_id[i],val_id_2$ix[1:maxNum]]
    l2_2 = IGM_signals_without_noise[id_set$grp_post_id[i],val_id_2$ix[1:maxNum]]
    
    v1 = c(v1,c(l1_1,l1_2))
    v2 = c(v2,c(l2_1,l2_2))
   
    s_name = sprintf('S%d',i)
    l3 = rep(s_name,2*maxNum) 
    v3 = c(v3,l3)
    
    l4 = c(rep('up',maxNum),rep('dn',maxNum))
    v4 = c(v4,l4)
    
    val_id_set = union(val_id_1$ix[1:maxNum],val_id_2$ix[1:maxNum])
    max_sig_pre_post_change_id_set = union(max_sig_pre_post_change_id_set, val_id_set) 
  }
  
  max_change$pre_sig_y     = v1
  max_change$post_sig_y    = v2
  max_change$pre_sig_x     = rep(0,length(v1))
  max_change$post_sig_x    = rep(1,length(v2))
  
  max_change$sample = v3
  max_change$up_dn       = v4
  
  max_change = as.data.frame(max_change)
  
  p = ggplot(max_change) + 
    geom_segment(aes(x=pre_sig_x,y=pre_sig_y,xend = post_sig_x,yend = post_sig_y, color = up_dn)) + 
    facet_wrap(~ sample, nrow = 2)
  
  fName = sprintf("./plot_all_sample_max_change_sigs_%s.pdf",type)
  pdf(file=fName,width=9,height=9)
  p=grid.arrange(p,nrow = 1)
  dev.off()
  
  ##################################
  ### paired t-test for pre and post
  probe_set = 1:length(probe_names)
  
  id_pre_post_thresholded = paired_t_test_pre_post(analysis_name,IGM_signals_without_noise,iso_ID_name_arr,id_set,n,type='IGM',threshold = 0.05,test_type=0)

  isoform_id = NULL
  isoform_id$set1 = max_sig_pre_post_change_id_set
  isoform_id$set2 = id_pre_post_thresholded
  isoform_id$signal = IGM_signals_without_noise
  isoform_id$iso_ID_name_arr = iso_ID_name_arr
  
  return(isoform_id)
}

paired_t_test_pre_post <- function(analysis_name,IGM_signals_without_noise,
                                       iso_ID_name_arr,id_set,n,type,threshold,test_type){
  
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

### type = 1 result with t-test with 0.05 significance
### type = 2 result with high-variable set of isoforms
main_run <-function(type){
  
  setwd("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4")
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/GPR_data_input.R")
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/glmnet_train.R")
  source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/use_psygenet.R")
  load_libs()
  
  maxNum = 1000
  type = 'dup1'
  isoform_dup1 = find_most_changing_isoform_id(type, maxNum)
  
  type = 'dup2'
  isoform_dup2 = find_most_changing_isoform_id(type, maxNum)  

  set1_intersection = intersect(isoform_dup1$set1,isoform_dup2$set1)
  set2_intersection = intersect(isoform_dup1$set2,isoform_dup2$set2)
  
  type = 2 ### wih t-test
  #### either through set1_intersection or set2_intersection
  if(type == 1){
    X = isoform_dup1$signal[,set1_intersection]
  }
  if(type == 2){
    X = isoform_dup1$signal[,set2_intersection]
  }
  
  iso_ID_name_arr = isoform_dup1$iso_ID_name_arr
  
  list_genes = run_EN(X,iso_ID_name_arr)

  
  library( psygenet2r )
  
  kk=2
  #for(kk in 1:length(list_genes)){
    genesOfInterest = unique(list_genes[[i]])
    m1 = psygenetGene(gene = genesOfInterest, database = "ALL", verbose  = FALSE, warnings = FALSE)
    fName = sprintf("./gene_disease_%d_type_%d.pdf",kk,type)
    pdf(fName,width=10,height=12)
    plot( m1 )
    dev.off()
    
    fName = sprintf("./gene_disease_hmap_%d_type_%d.pdf",kk,type)
    pdf(fName,width=11,height=9)
    plot( m1, type="GDA heatmap")
    dev.off()
    
  #}
    
    #################
    
    isoID = (lapply(list_genes[[kk]],function(x){(iso_ID_name_arr$iso_id[(iso_ID_name_arr$gene_name == x)])}))
    
    isoID_filtered = unlist(lapply(isoID,function(x){x[length(x) < 10]}))
    
    Z = X
    indx_set = unlist(lapply(isoID_filtered,function(x){which(colnames(Z)==x)}))
    
    #################
    
    source("~/Desktop/ToDo/Pablo/plasma_sample_data/newCode_Dec4/BGM_EN.R")
    init_library()
    
    Y1 = Y_case = X[1:16,indx_set]
    Y2 = Y_comvert = X[17:25,indx_set]
    
    ##### run my BGM_EN network
    n_burn_in = 200
    n_mcmc_samples = 500
    
    tic("total MCMC")
    data_with_init_with_MCMC_samples_Y1 = MCMC_sampling_BGN_EN(Y1,n_burn_in,n_mcmc_samples)
    data_with_init_with_MCMC_samples_Y2 = MCMC_sampling_BGN_EN(Y2,n_burn_in,n_mcmc_samples)
    
    toc()
    beta_bar_sym_Y1 = data_with_init_with_MCMC_samples_Y1$beta_bar_sym
    fName = sprintf("./data_with_init_with_MCMC_samples_Y1_type_%d.RData",type)
    save(data_with_init_with_MCMC_samples_Y1,file=fName)
    
    beta_bar_sym_Y2 = data_with_init_with_MCMC_samples_Y2$beta_bar_sym
    fName = sprintf("./data_with_init_with_MCMC_samples_Y2_type_%d.RData",type)
    save(data_with_init_with_MCMC_samples_Y2,file=fName)
    
    ### 148 x 148
    generate_and_plot_network(beta_bar_sym_Y1,type,control_convert=1,th=0.125)
    generate_and_plot_network(beta_bar_sym_Y2,type,control_convert=0,th=0.125)
}

