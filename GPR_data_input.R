####
## reading the files
####

load_libs <- function(){
  
  library(limma)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(RColorBrewer)
  library(corrplot)
  library(multiClust)
  library(reshape2)
  library(preprocessCore)
  library(MASS)
  library(car)
  library(gplots)
  library(VennDiagram)
  library(org.Hs.eg.db)
  library(scales)
  library(tidyverse)
  library(glmnet)
  library(randomForest)
}

sym_diff <- function(a,b) setdiff(union(a,b), intersect(a,b))

read_GPR_files <- function(){
    
  #### IGM (cy5 and 635nm) - IGM channel  - F635 Median/B635 Medain
  #### IGG (cy3 and 532nm) - IGG channel  - F532 Median/B532 Median
  #### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1810453/
  
  files_full_path = list.files("~/Desktop/ToDo/Pablo/plasma_sample_data/GPR_files",full.names = TRUE)
  ### read the red and green channel median value both foreground and background
  IGM_IGG_med = read.maimages(files_full_path,"genepix.median")
 
  save(IGM_IGG_med,file='./IGM_IGG_med.Rdata')
  
  return(IGM_IGG_med) 
}
  
signal_extract_from_GPR_files <-function(IGM_IGG_med){

  fNames = dir("~/Desktop/ToDo/Pablo/plasma_sample_data/GPR_files")
  
  gpr_id = which(grepl(".gpr",fNames) == TRUE)
  fNames = fNames[gpr_id]
  
  n_sample = length(IGM_IGG_med$targets$FileName)
  ### remove the isoform ID which are empty/Control/buffer
  ### take only JHU-prefix isoform
  id_JHU = which(grepl("JHU",IGM_IGG_med$genes$ID)==TRUE)
  
  isoform_id = unique(IGM_IGG_med$genes$ID[id_JHU])
  L1 = length(isoform_id)
  
  
  #### each column is a sample
  mat_IGM_F_by_B_dup_1 = matrix(0,nrow = L1, ncol = n_sample)
  mat_IGM_F_by_B_dup_2 = matrix(0,nrow = L1, ncol = n_sample)
  
  mat_IGG_F_by_B_dup_1 = matrix(0,nrow = L1, ncol = n_sample)
  mat_IGG_F_by_B_dup_2 = matrix(0,nrow = L1, ncol = n_sample)
  
  IGM_IGG_med$targets$FileName = fNames
  col_names_tmp = IGM_IGG_med$targets$FileName[1:n_sample]
  
  map_id = read.table("~/Desktop/ToDo/Pablo/plasma_sample_data/samples_map.txt")
  

  R_1 = vector()
  R_2 = vector()
  Rb_1 = vector()
  Rb_2 = vector()
  
  G_1 = vector()
  G_2 = vector()
  Gb_1 = vector()
  Gb_2 = vector()
  
  col_names = NULL
  for (i in 1:n_sample){
    id1 = which(map_id[,1] == col_names_tmp[i])
    col_names[i] = as.character(map_id[id1,2])
  }
  
  row_names = isoform_id[1:L1]
  
  dimnames(mat_IGM_F_by_B_dup_1) = list(row_names,col_names)
  dimnames(mat_IGM_F_by_B_dup_2) = list(row_names,col_names)
  
  dimnames(mat_IGG_F_by_B_dup_1) = list(row_names,col_names)
  dimnames(mat_IGG_F_by_B_dup_2) = list(row_names,col_names)
  
  
  ##### get the corresponding two indices which are same
  v = lapply(isoform_id,function(x) which(IGM_IGG_med$genes$ID == x))
  
  ####
  ## mat_IGM_F_dup_1 is 21154 x 25
  ## mat_IGM_B_dup_2 is 21154 x 25
  
  ## mat_IGG_F_dup_1 is 21154 x 25
  ## mat_IGG_B_dup_2 is 21154 x 25
  
  sample_all_signals = list()
  
  for (i in 1:n_sample){
    
    mat_IGM_F_by_B_dup_1[,i] = as.vector(unlist(lapply(v,function(x) {IGM_IGG_med$R[x[1],i]/IGM_IGG_med$Rb[x[1],i]} )))
    mat_IGM_F_by_B_dup_2[,i] = as.vector(unlist(lapply(v,function(x) {IGM_IGG_med$R[x[2],i]/IGM_IGG_med$Rb[x[2],i]} )))
    
    mat_IGG_F_by_B_dup_1[,i] = as.vector(unlist(lapply(v,function(x) {IGM_IGG_med$G[x[1],i]/IGM_IGG_med$Gb[x[1],i]} )))
    mat_IGG_F_by_B_dup_2[,i] = as.vector(unlist(lapply(v,function(x) {IGM_IGG_med$G[x[2],i]/IGM_IGG_med$Gb[x[2],i]} )))
    
    R_1 = as.vector(unlist(lapply(v,function(x) IGM_IGG_med$R[x[1],i])))
    R_2 = as.vector(unlist(lapply(v,function(x) IGM_IGG_med$R[x[2],i])))
    
    Rb_1 = as.vector(unlist(lapply(v,function(x) IGM_IGG_med$Rb[x[1],i])))
    Rb_2 = as.vector(unlist(lapply(v,function(x) IGM_IGG_med$Rb[x[2],i])))
    
    G_1 = as.vector(unlist(lapply(v,function(x) IGM_IGG_med$G[x[1],i])))
    G_2 = as.vector(unlist(lapply(v,function(x) IGM_IGG_med$G[x[2],i])))
    
    Gb_1 = as.vector(unlist(lapply(v,function(x) IGM_IGG_med$Gb[x[1],i])))
    Gb_2 = as.vector(unlist(lapply(v,function(x) IGM_IGG_med$Gb[x[2],i])))
    
    df_R_G_Rb_Gb = data.frame(R_1,R_2,Rb_1,Rb_2,G_1,G_2,Gb_1,Gb_2)
                              
    sample_all_signals[[i]] = df_R_G_Rb_Gb
    
    cat("sample ", i, " is done!!\n")
  }

  save(sample_all_signals,file="./sample_all_signals.RData")
  
  sample_F_by_B_signals = NULL
  
  sample_F_by_B_signals$mat_IGM_F_by_B_dup_1 = mat_IGM_F_by_B_dup_1
  sample_F_by_B_signals$mat_IGM_F_by_B_dup_2 = mat_IGM_F_by_B_dup_2
  
  sample_F_by_B_signals$mat_IGG_F_by_B_dup_1 = mat_IGG_F_by_B_dup_1
  sample_F_by_B_signals$mat_IGG_F_by_B_dup_2 = mat_IGG_F_by_B_dup_2
    
  sample_F_by_B_signals$row_names = row_names
  sample_F_by_B_signals$col_names = col_names
  
  save(sample_F_by_B_signals,file="./sample_F_by_B_signals.RData")
  
  return(sample_F_by_B_signals)
}


####################
map_isoform_ID_gene_name <- function(IGM_IGG_med){

  id_JHU = which(grepl("JHU",IGM_IGG_med$genes$ID)==TRUE)
  
  uniq_isoform_ID = unique(IGM_IGG_med$gene$ID[id_JHU])
  L2 = length(uniq_isoform_ID)
  
  iso_ID_name_arr = NULL
  iso_id_list = vector()
  gene_name_list = vector()
  
  for (i in 1:L2){
    iso_id = uniq_isoform_ID[i]
    indx = which(IGM_IGG_med$gene$ID == iso_id)
    all_gene_names = IGM_IGG_med$gene$Name[indx]
    if(length(unique(all_gene_names))>1){
      cat("Error !! \n")
    }else{
      gene_name = unique(all_gene_names)
    }
    
    iso_id_list[i] = iso_id
    gene_name_list[i] = gene_name
    
  }
  
  iso_ID_name_arr$iso_id = iso_id_list
  iso_ID_name_arr$gene_name = gene_name_list
  
  
  save(iso_ID_name_arr,file="iso_ID_name_arr.RData")

  return(iso_ID_name_arr)
}


gen_QN_signals <-function(sample_F_by_B_signals){
  
  row_names = sample_F_by_B_signals$row_names 
  col_names = sample_F_by_B_signals$col_names
  
  mat_IGM_F_by_B_dup_1 = sample_F_by_B_signals$mat_IGM_F_by_B_dup_1 
  mat_IGM_F_by_B_dup_2 = sample_F_by_B_signals$mat_IGM_F_by_B_dup_2
  
  mat_IGG_F_by_B_dup_1 = sample_F_by_B_signals$mat_IGG_F_by_B_dup_1
  mat_IGG_F_by_B_dup_2 = sample_F_by_B_signals$mat_IGG_F_by_B_dup_2
  
  ##### Quantile Normalization
  
  mat_IGM_F_by_B_dup_1_QN = normalize.quantiles((mat_IGM_F_by_B_dup_1))
  mat_IGM_F_by_B_dup_2_QN = normalize.quantiles((mat_IGM_F_by_B_dup_2))
  mat_IGG_F_by_B_dup_1_QN = normalize.quantiles((mat_IGG_F_by_B_dup_1))
  mat_IGG_F_by_B_dup_2_QN = normalize.quantiles((mat_IGG_F_by_B_dup_2))
  
  dimnames(mat_IGM_F_by_B_dup_1_QN) = list(row_names, col_names)
  dimnames(mat_IGM_F_by_B_dup_2_QN) = list(row_names, col_names)
  dimnames(mat_IGG_F_by_B_dup_1_QN) = list(row_names, col_names)
  dimnames(mat_IGG_F_by_B_dup_2_QN) = list(row_names, col_names)
  
  
  QN_signals = NULL
  
  QN_signals$mat_IGM_F_by_B_dup_1 = mat_IGM_F_by_B_dup_1_QN
  QN_signals$mat_IGM_F_by_B_dup_2 = mat_IGM_F_by_B_dup_2_QN
  QN_signals$mat_IGG_F_by_B_dup_1 = mat_IGG_F_by_B_dup_1_QN
  QN_signals$mat_IGG_F_by_B_dup_2 = mat_IGG_F_by_B_dup_2_QN
  
  QN_signals$row_names = row_names
  QN_signals$col_names = col_names
  
  
  return(QN_signals)

}
  
gen_max_signals <- function(signals,sel_type){

  row_names = signals$row_names 
  col_names = signals$col_names
  
  L1 = length(row_names)
  n_sample = length(col_names)
  
  mat_IGM_F_by_B_dup_1 = signals$mat_IGM_F_by_B_dup_1 
  mat_IGM_F_by_B_dup_2 = signals$mat_IGM_F_by_B_dup_2
  
  mat_IGG_F_by_B_dup_1 = signals$mat_IGG_F_by_B_dup_1
  mat_IGG_F_by_B_dup_2 = signals$mat_IGG_F_by_B_dup_2
  
  if((sel_type == 'max')|(sel_type == 'avg')) {
    mat_IGM_F_by_B_tmp = array(c(mat_IGM_F_by_B_dup_1,mat_IGM_F_by_B_dup_2),
                                dim = c(L1,n_sample,2))
    mat_IGG_F_by_B_tmp = array(c(mat_IGG_F_by_B_dup_1,mat_IGG_F_by_B_dup_2),
                               dim = c(L1,n_sample,2))
  }
  
  if(sel_type == 'dup1'){
    mat_IGM_F_by_B_tmp = array(c(mat_IGM_F_by_B_dup_1,mat_IGM_F_by_B_dup_1),
                               dim = c(L1,n_sample,2))
    mat_IGG_F_by_B_tmp = array(c(mat_IGG_F_by_B_dup_1,mat_IGG_F_by_B_dup_1),
                               dim = c(L1,n_sample,2))
  }
  
  if(sel_type == 'dup2'){
    mat_IGM_F_by_B_tmp = array(c(mat_IGM_F_by_B_dup_2,mat_IGM_F_by_B_dup_2),
                               dim = c(L1,n_sample,2))
    mat_IGG_F_by_B_tmp = array(c(mat_IGG_F_by_B_dup_2,mat_IGG_F_by_B_dup_2),
                               dim = c(L1,n_sample,2))
  }
  
  if((sel_type == 'max')|(sel_type == 'dup1')|(sel_type == 'dup2')){
    mat_IGM_F_by_B_max = apply(mat_IGM_F_by_B_tmp,c(1,2),'max')
    mat_IGG_F_by_B_max = apply(mat_IGG_F_by_B_tmp,c(1,2),'max')
  }
  if(sel_type == 'avg'){
    mat_IGM_F_by_B_max = apply(mat_IGM_F_by_B_tmp,c(1,2),'mean')
    mat_IGG_F_by_B_max = apply(mat_IGG_F_by_B_tmp,c(1,2),'mean')
  }
  
  
  dimnames(mat_IGM_F_by_B_max) = list(row_names, col_names)
  dimnames(mat_IGG_F_by_B_max) = list(row_names, col_names)
  
  max_signals = NULL
  max_signals$mat_IGM_F_by_B = mat_IGM_F_by_B_max
  max_signals$mat_IGG_F_by_B = mat_IGG_F_by_B_max
  max_signals$row_names = row_names
  max_signals$col_names = col_names
  
  return(max_signals)
  
}

plot_duplicate_signals <-function(sample_all_signals){  
  
  sample = sample_all_signals
  
  n_sample = length(samples)
  
  for (i in 1:n_sample){
    s = sprintf("./dup/IGM_dup_F_by_B_%d.png",i)
    png(s,width=2000,height=2000,res=200)
    IGM_tmp = ggplot(data=sample[[i]], aes(x=sample[[i]]$R_1/sample[[i]]$Rb_1, y=sample[[i]]$R_2/sample[[i]]$Rb_2, group=1)) +geom_point()+ 
      xlab("IGM-Isoform 1 F/B")+ylab("IGM-Isoform 2 F/B")+ggtitle("S1")+
      theme(plot.title = element_text(hjust = 0.5))+
      geom_abline(slope=1, intercept=0,colour='#E41A1C')
    
    yy = grid.arrange(IGM_tmp,nrow=1)
    op = par(no.readonly=TRUE)
    par(op,pty="s")
    dev.off()
  }
  
  for (i in 1:n_sample){
    s = sprintf("./dup/IGG_dup_F_by_B_%d.png",i)
    png(s,width=2000,height=2000,res=200)
    IGG_tmp = ggplot(data=sample[[i]], aes(x=sample[[i]]$G_1/sample[[i]]$Gb_1, y=sample[[i]]$G_2/sample[[i]]$Gb_2, group=1)) +geom_point()+ 
      xlab("IGG-Isoform 1 F/B")+ylab("IGG-Isoform 2 F/B")+ggtitle("S1")+
      theme(plot.title = element_text(hjust = 0.5))+
      geom_abline(slope=1, intercept=0, colour='#E41A1C')
    
    yy = grid.arrange(IGG_tmp,nrow=1)
    op = par(no.readonly=TRUE)
    par(op,pty="s")
    dev.off()
  }

}  

remove_noisy_probes <- function(signals){
  
  set_rem_IGM = NULL
  set_rem_IGG = NULL
  
  n_sample = length(signals$col_names)
  
  for(i in 1:n_sample){
    set_rem_IGM = union(set_rem_IGM, which(signals$mat_IGM_F_by_B[,i] < 1))
  }
  
  for(i in 1:n_sample){
    set_rem_IGG = union(set_rem_IGG, which(signals$mat_IGG_F_by_B[,i] < 1))
  }
  set_rem = NULL
  set_rem$IGM = sort(set_rem_IGM)
  set_rem$IGG = sort(set_rem_IGG)
  
  return(set_rem)
}
  
get_threshold_probes <- function(signals,sd_cutoff_factor){

  mat_IGM_F_by_B = signals$mat_IGM_F_by_B
  mat_IGG_F_by_B = signals$mat_IGG_F_by_B
  row_names = signals$row_names 
  col_names = signals$col_names 
  n_sample = length(col_names)
  
  ### IGM channel
  IGM_pos_prob_len = vector()
  IGM_probe_set = list()
  
  png(file="./all_sample_IGM_F_by_B_density_pos_probes.png",width=2000,height=2000,res=200)
  par(mfrow=c(5,5),mar=c(5,6,4,2)+0.1)
  
  
  for (i in 1:n_sample){
    
    indx_less_one = which(mat_IGM_F_by_B[,i] <= 1.0)
    tmp_data1 = mat_IGM_F_by_B[indx_less_one,i]-1
    tmp_data2 = -tmp_data1
    tmp_data = c(tmp_data1,tmp_data2)+1
    f_distr = fitdistr(tmp_data,'normal')
    mean_est = f_distr$estimate['mean']
    sd_est = f_distr$estimate['sd']
    
    cutoff = mean_est + sd_cutoff_factor*sd_est 
    indx_greater_cutoff = which(mat_IGM_F_by_B[,i] > cutoff)
    IGM_pos_prob_len[i] = length(indx_greater_cutoff)
    
    densityPlot(log(mat_IGM_F_by_B[indx_greater_cutoff,i]),
                main=col_names[i],xlab = 'IGM F/B (log)') 
    
    IGM_probe_set[[i]] = indx_greater_cutoff
    
  }
  dev.off()
  
  ### IGG channel
  IGG_pos_prob_len = vector()
  IGG_probe_set = list()
  
  png(file="./all_sample_IGG_F_by_B_density_pos_probes.png",width=2000,height=2000,res=200)
  par(mfrow=c(5,5),mar=c(5,6,4,2)+0.1)
  
  for (i in 1:n_sample){
    
    indx_less_one = which(mat_IGG_F_by_B[,i] <= 1.0)
    tmp_data1 = mat_IGG_F_by_B[indx_less_one,i]-1
    tmp_data2 = -tmp_data1
    tmp_data = c(tmp_data1,tmp_data2)+1
    f_distr = fitdistr(tmp_data,'normal')
    mean_est = f_distr$estimate['mean']
    sd_est = f_distr$estimate['sd']
    
    cutoff = mean_est + sd_cutoff_factor*sd_est
    indx_greater_cutoff = which(mat_IGG_F_by_B[,i] > cutoff)
    IGG_pos_prob_len[i] = length(indx_greater_cutoff)
    
    densityPlot(log(mat_IGG_F_by_B[indx_greater_cutoff,i]),
                main=col_names[i],xlab = 'IGG F/B (log)') 
    
    IGG_probe_set[[i]] = indx_greater_cutoff
    
  }
  dev.off()
  
  
  IGM_probe_set_union = sort(Reduce(union,IGM_probe_set))
  IGG_probe_set_union = sort(Reduce(union,IGG_probe_set))
  
  IGM_probe_set_intersect = sort(Reduce(intersect,IGM_probe_set))
  IGG_probe_set_intersect = sort(Reduce(intersect,IGG_probe_set))
  
  
  probe_set = NULL
  probe_set$IGM_probe_set = IGM_probe_set
  probe_set$IGG_probe_set = IGG_probe_set
  
  probe_set$IGM_probe_set_union = IGM_probe_set_union
  probe_set$IGG_probe_set_union = IGG_probe_set_union
  
  probe_set$IGM_probe_set_intersect = IGM_probe_set_intersect
  probe_set$IGG_probe_set_intersect = IGG_probe_set_intersect
  
  return(probe_set)
  
}

subject_and_time <- function(sample_names){
  #sample_names = names(average_sig_IGM_sel)
  #n_sample = dim(IGM_signals_sel)[2]
  n_sample = length(sample_names)
  all_tokens = unlist(lapply(sample_names,function(x) strsplit(x,'_')[[1]]))
  all_subject_rep = all_tokens[seq(1,n_sample*3,by=3)]
  all_subj = unique(all_subject_rep)
  n_subj = length(all_subj)
  all_times_rep = unlist(lapply(sample_names,function(x) sprintf('%s_%s',strsplit(x,'_')[[1]][2],
                                                                 strsplit(x,'_')[[1]][3])))
  all_times = unique(all_times_rep)
  n_times = length(all_times)
  
  subj_time_info = NULL
  subj_time_info$all_times = all_times
  subj_time_info$all_subj = all_subj
  
  return(subj_time_info)

}

overall_immune_signal <- function(IGM_signals_sel,IGG_signals_sel,sample_names){
  average_sig_IGM_sel = apply(IGM_signals_sel,2,mean)
  average_sig_IGG_sel = apply(IGG_signals_sel,2,mean)
  
  n_sample = length(sample_names)
  all_tokens = unlist(lapply(sample_names,function(x) strsplit(x,'_')[[1]]))
  all_subject_rep = all_tokens[seq(1,n_sample*3,by=3)]
  all_subj = unique(all_subject_rep)
  n_subj = length(all_subj)
  all_times_rep = unlist(lapply(sample_names,function(x) sprintf('%s_%s',strsplit(x,'_')[[1]][2],
                                                                 strsplit(x,'_')[[1]][3])))
  
  
  overall_sig = NULL
  overall_sig$IGM = average_sig_IGM_sel
  overall_sig$IGG = average_sig_IGG_sel
  overall_sig$subj = all_subject_rep
  overall_sig$times = all_times_rep

  overall_sig = as.data.frame(overall_sig)
  return(overall_sig)  
}

plot_overall_signal <- function(overall_sig_tmp,s){
  
  L1 = dim(overall_sig_tmp)[1]
  
  overall_sig = NULL 
  overall_sig$sig = c(overall_sig_tmp$IGM,overall_sig_tmp$IGG)
  overall_sig$times = c(overall_sig_tmp$times,overall_sig_tmp$times)
  overall_sig$subj = c(overall_sig_tmp$subj,overall_sig_tmp$subj)
  overall_sig$type = c(rep('IGM',L1),rep('IGG',L1))
  
  overall_sig$times = factor(overall_sig_tmp$times,levels=c('m_14','m_12','m_11','m_10',
                                                        'm_8','m_7','m_6','m_5', 
                                                        'm_3','m_2','d_0','p_3',
                                                        'p_4','p_5','p_6','p_14'))
  
  overall_sig$type = factor(overall_sig$type,levels=c('IGM','IGG'))
  
  overall_sig = as.data.frame(overall_sig)
  
  fName = sprintf("./all_sample_IGM_IGG_overall_mean_sig_%s.png",s)
  t1 = sprintf('Mean signal across all selected probes (%s)',s)
  
  png(file=fName,width=3500,height=1500,res=270)
  p1 = ggplot(overall_sig,aes(x=times,y=sig,color = factor(subj),group = subj)) + 
    geom_line()+
    geom_point()+
    facet_grid(.~type)+
    scale_x_discrete(name='Timepoint')+
    scale_y_continuous(name=t1)+
    scale_color_discrete("Subject")+
    geom_vline(xintercept=11, color = "black",linetype='dashed',size=0.5)
  
  p=grid.arrange(p1,nrow = 1)
  dev.off()
  
  
}


plot_one_gene_signal <- function(overall_sig_tmp,s){
  
  L1 = dim(overall_sig_tmp)[1]
  
  overall_sig = NULL 
  overall_sig$sig = overall_sig_tmp$IGM
  overall_sig$times = overall_sig_tmp$times
  overall_sig$subj = overall_sig_tmp$subj
  overall_sig$type = rep('IGM',L1)
  
  overall_sig$times = factor(overall_sig_tmp$times,levels=c('m_14','m_12','m_11','m_10',
                                                            'm_8','m_7','m_6','m_5', 
                                                            'm_3','m_2','d_0','p_3',
                                                            'p_4','p_5','p_6','p_14'))
  
  overall_sig = as.data.frame(overall_sig)
  
  fName = sprintf("./all_sample_IGM_sig_%s.pdf",s)
  t1 = sprintf('Signal across subjects for the selected probe (%s)',s)
  
  pdf(file=fName,width=10,height=7)
  p1 = ggplot(overall_sig,aes(x=times,y=sig,color = factor(subj), group = subj)) + 
    geom_line()+
    geom_point()+
    scale_x_discrete(name='Timepoints')+
    scale_y_continuous(name=t1)+
    scale_color_discrete("Subject")+
    geom_vline(xintercept=11, color = "black",linetype='dashed',size=0.5)+
    theme(axis.text = element_text(size = 12),legend.title=element_text(size=14),
        legend.text=element_text(size=14),legend.position="top",
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13))
  
  p=grid.arrange(p1,nrow = 1)
  dev.off()
  
  
}


#### BRAIN
plot_brain_related_signal <- function(max_signals,brain_GN_list,
                                      iso_ID_name_arr,subj_time_info){
  
  all_subj_brain_probe = list()
  kk = 0
  for(subj in subj_time_info$all_subj){
  
    ts_id = which(grepl(subj,max_signals$col_names) == TRUE)
    ts = max_signals$col_names[ts_id]
    
    IGM_sig = max_signals$mat_IGM_F_by_B[,c(ts_id)]
    IGG_sig = max_signals$mat_IGG_F_by_B[,c(ts_id)]
    
    L1 = length(brain_GN_list)
    
    probe_names = NULL
    
    for(i in 1:L1){
      id_tmp = which(iso_ID_name_arr$gene_name == brain_GN_list[i])
      probe_names = c(probe_names,iso_ID_name_arr$iso_id[id_tmp])
    }
    
    L2 = length(probe_names)
    indx_set = vector()
    brain_probe = NULL
    indx_set = vector()
    
    tmp_L = length(ts_id)
    k = 1
    
    for(i in 1:L2){
      indx_set[i] = which(max_signals$row_names == probe_names[i])
      
      brain_probe$sig[k:(k+tmp_L-1)] = IGM_sig[indx_set[i],]
      brain_probe$type[k:(k+tmp_L-1)] = rep('IGM',tmp_L)
      brain_probe$ts[k:(k+tmp_L-1)] = ts
      brain_probe$probe[k:(k+tmp_L-1)] = probe_names[i]   
      brain_probe$sd_igm[i] = sd(brain_probe$sig[k:(k+tmp_L-1)])
      
      k = k+tmp_L
      brain_probe$sig[k:(k+tmp_L-1)] = IGG_sig[indx_set[i],]
      brain_probe$type[k:(k+tmp_L-1)] = rep('IGG',tmp_L)
      brain_probe$ts[k:(k+tmp_L-1)] = ts
      brain_probe$probe[k:(k+tmp_L-1)] = probe_names[i]   
      brain_probe$sd_igg[i] = sd(brain_probe$sig[k:(k+tmp_L-1)])

      brain_probe$pn[i] = probe_names[i]
      id_tmp = which(iso_ID_name_arr$iso_id == brain_probe$pn[i])
      if(length(id_tmp) > 1){
        id_tmp = id_tmp[1]
      }
      brain_probe$gn[i] = iso_ID_name_arr$gene_name[id_tmp]
      
      k = k+tmp_L
      
    }
    
    brain_probe_plot = as.data.frame(brain_probe)
    brain_probe_plot$ts = factor(brain_probe_plot$ts,levels=c(ts))
    brain_probe_plot$type = factor(brain_probe_plot$type,levels=c('IGM','IGG'))
    
    fname1 = sprintf("./%s_IGM_IGG_brain_sig.png", subj)
    title = sprintf("Brain probe signal across all time point for subject %s",subj)
    
    png(file=fname1,width=3500,height=1500,res=270)
    p1 = ggplot(brain_probe_plot,aes(x=ts,y=sig,group=probe)) + 
      geom_line()+
      geom_point()+
      facet_grid(.~type)+
      scale_x_discrete(name='Timepoints')+
      scale_y_continuous(name='Signal value (in log 10)',
                         trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))+
      #scale_y_continuous(name='Signal value')+
      #coord_trans(y="log10")+
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
      #scale_color_discrete("Samples")+
      #geom_vline(xintercept=11, color = "black",linetype='dashed',size=0.5)
    
    p=grid.arrange(p1,nrow = 1)
    dev.off()
    
    kk = kk+1
    all_subj_brain_probe$sig[[kk]] = brain_probe$sig
    all_subj_brain_probe$type[[kk]] = brain_probe$type
    all_subj_brain_probe$ts[[kk]] = brain_probe$ts
    all_subj_brain_probe$probe[[kk]] = brain_probe$probe
    all_subj_brain_probe$std_dev_igm[[kk]] = brain_probe$sd_igm
    all_subj_brain_probe$std_dev_igg[[kk]] = brain_probe$sd_igg
    all_subj_brain_probe$probe_name[[kk]] = brain_probe$pn
    all_subj_brain_probe$gene_name[[kk]] = brain_probe$gn
    
    cat("subject ", subj, " is done!!\n")
  }
  return(all_subj_brain_probe)
}

plot_top_varying_brain_related_signal <- function(max_signals,all_subj_brain_probe,
                                      iso_ID_name_arr,subj_time_info,topK){
  
  
  top_varying_pn_gn = NULL
  i = 0
  for (subj in subj_time_info$all_subj){
    
    i = i+1
    
    id_sorted_igm = sort(all_subj_brain_probe$std_dev_igm[[i]],
                     decreasing = T, index.return=T)
    
    id_sorted_igg = sort(all_subj_brain_probe$std_dev_igg[[i]],
                     decreasing = T, index.return=T)
    
    
    id_tmp_igm = id_sorted_igm$ix[1:topK]
    id_tmp_igg = id_sorted_igg$ix[1:topK]
    
    #id_tmp_igm = which(all_subj_brain_probe$std_dev_igm[[i]] > 4.0) ##
    #id_tmp_igg = which(all_subj_brain_probe$std_dev_igg[[i]] > 4.0) ## 
    
    cat('number of highly varying sig: igm: ', length(id_tmp_igm),
        ' igg: ',length(id_tmp_igg),'\n')
    
    top_varying_pn_gn$pn_igm[[i]] = all_subj_brain_probe$probe_name[[i]][id_tmp_igm]
    top_varying_pn_gn$pn_igg[[i]] = all_subj_brain_probe$probe_name[[i]][id_tmp_igg]
    
    top_varying_pn_gn$gn_igm[[i]] = all_subj_brain_probe$gene_name[[i]][id_tmp_igm]
    top_varying_pn_gn$gn_igg[[i]] = all_subj_brain_probe$gene_name[[i]][id_tmp_igg]
    
  }
  
  kk = 0
  for (subj in subj_time_info$all_subj){
    
    ts_id = which(grepl(subj,max_signals$col_names) == TRUE)
    ts = max_signals$col_names[ts_id]
    
    IGM_sig = max_signals$mat_IGM_F_by_B[,c(ts_id)]
    IGG_sig = max_signals$mat_IGG_F_by_B[,c(ts_id)]
    
    kk = kk+1
    
    L2_igm = length(top_varying_pn_gn$pn_igm[[kk]])
    L2_igg = length(top_varying_pn_gn$pn_igg[[kk]])
    
    brain_probe_igm = NULL
    brain_probe_igg = NULL
    indx_set_igm = vector()
    indx_set_igg = vector()
    
    tmp_L = length(ts_id)
    k = 1
    
    for(j in 1:L2_igm){
      indx_set_igm[j] = which(max_signals$row_names ==  top_varying_pn_gn$pn_igm[[kk]][j])
      brain_probe_igm$sig[k:(k+tmp_L-1)] = IGM_sig[indx_set_igm[j],]
      brain_probe_igm$ts[k:(k+tmp_L-1)] = ts
      brain_probe_igm$probe[k:(k+tmp_L-1)] = top_varying_pn_gn$pn_igm[[kk]][j]   
      k = k+tmp_L
    }
    
    k = 1
    for(j in 1:L2_igg){
      indx_set_igg[j] = which(max_signals$row_names ==  top_varying_pn_gn$pn_igg[[kk]][j])
      brain_probe_igg$sig[k:(k+tmp_L-1)] = IGG_sig[indx_set_igg[j],]
      brain_probe_igg$ts[k:(k+tmp_L-1)] = ts
      brain_probe_igg$probe[k:(k+tmp_L-1)] = top_varying_pn_gn$pn_igg[[kk]][j]  
      k = k+tmp_L
    }
    
    brain_probe_plot_igm = as.data.frame(brain_probe_igm)
    brain_probe_plot_igm$ts = factor(brain_probe_plot_igm$ts,levels=c(ts))

    brain_probe_plot_igg = as.data.frame(brain_probe_igg)
    brain_probe_plot_igg$ts = factor(brain_probe_plot_igg$ts,levels=c(ts))
    
    fname12 = sprintf("./%s_IGM_IGG_brain_top_varying_sig_topK_%d.png", subj,topK)
    
    title1 = sprintf("Top Varying Brain probe signal (IGM) across all time point for subject %s",subj)
    title2 = sprintf("Top Varying Brain probe signal (IGG) across all time point for subject %s",subj)
    
    png(file=fname12,width=3500,height=1500,res=270)
    
    p1 = ggplot(brain_probe_plot_igm,aes(x=ts,y=sig,group=probe)) + 
      geom_line()+
      geom_point()+
      scale_x_discrete(name='Timepoints')+
      scale_y_continuous(name='Signal value (in log 10)',
                         trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))+
      #coord_trans(y="log10")+
      ggtitle(title1) + theme(plot.title = element_text(hjust = 0.5))
    
    p2 = ggplot(brain_probe_plot_igg,aes(x=ts,y=sig,group=probe)) + 
      geom_line()+
      geom_point()+
      scale_x_discrete(name='Timepoints')+
      scale_y_continuous(name='Signal value (in log 10)',
                         trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))+
      #coord_trans(y="log10")+
      ggtitle(title2) + theme(plot.title = element_text(hjust = 0.5))
    
    p=grid.arrange(p1,p2,nrow = 1)
    dev.off()
    
  }
    
  return(top_varying_pn_gn)   
}


############ PLACENTA ######
plot_placenta_related_signal <- function(max_signals,placenta_GN_list,
                                      iso_ID_name_arr,subj_time_info){
  
  all_subj_placenta_probe = list()
  kk = 0
  for(subj in subj_time_info$all_subj){
    
    ts_id = which(grepl(subj,max_signals$col_names) == TRUE)
    ts = max_signals$col_names[ts_id]
    
    IGM_sig = max_signals$mat_IGM_F_by_B[,c(ts_id)]
    IGG_sig = max_signals$mat_IGG_F_by_B[,c(ts_id)]
    
    L1 = length(placenta_GN_list)
    
    probe_names = NULL
    
    for(i in 1:L1){
      id_tmp = which(iso_ID_name_arr$gene_name == placenta_GN_list[i])
      probe_names = c(probe_names,iso_ID_name_arr$iso_id[id_tmp])
    }
    
    L2 = length(probe_names)
    indx_set = vector()
    
    placenta_probe = NULL
    indx_set = vector()
    
    tmp_L = length(ts_id)
    k = 1
    
    for(i in 1:L2){
      indx_set[i] = which(max_signals$row_names == probe_names[i])
      
      placenta_probe$sig[k:(k+tmp_L-1)] = IGM_sig[indx_set[i],]
      placenta_probe$type[k:(k+tmp_L-1)] = rep('IGM',tmp_L)
      placenta_probe$ts[k:(k+tmp_L-1)] = ts
      placenta_probe$probe[k:(k+tmp_L-1)] = probe_names[i]   
      placenta_probe$sd_igm[i] = sd(placenta_probe$sig[k:(k+tmp_L-1)])
      
      k = k+tmp_L
      placenta_probe$sig[k:(k+tmp_L-1)] = IGG_sig[indx_set[i],]
      placenta_probe$type[k:(k+tmp_L-1)] = rep('IGG',tmp_L)
      placenta_probe$ts[k:(k+tmp_L-1)] = ts
      placenta_probe$probe[k:(k+tmp_L-1)] = probe_names[i]   
      placenta_probe$sd_igg[i] = sd(placenta_probe$sig[k:(k+tmp_L-1)])
      
      placenta_probe$pn[i] = probe_names[i]
      id_tmp = which(iso_ID_name_arr$iso_id == placenta_probe$pn[i])
      if(length(id_tmp) > 1){
        id_tmp = id_tmp[1]
      }
      placenta_probe$gn[i] = iso_ID_name_arr$gene_name[id_tmp]
      
      k = k+tmp_L
      
    }
    
    placenta_probe_plot = as.data.frame(placenta_probe)
    placenta_probe_plot$ts = factor(placenta_probe_plot$ts,levels=c(ts))
    placenta_probe_plot$type = factor(placenta_probe_plot$type,levels=c('IGM','IGG'))
    
    fname1 = sprintf("./%s_IGM_IGG_placenta_sig.png", subj)
    title = sprintf("Placenta probe signal across all time point for subject %s",subj)
    
    png(file=fname1,width=3500,height=1500,res=270)
    p1 = ggplot(placenta_probe_plot,aes(x=ts,y=sig,group=probe)) + 
      geom_line()+
      geom_point()+
      facet_grid(.~type)+
      scale_x_discrete(name='Timepoints')+
      scale_y_continuous(name='Signal value (in log 10)',
                         trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))+
      #scale_y_continuous(name='Signal value')+
      #coord_trans(y="log10")+
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    #scale_color_discrete("Samples")+
    #geom_vline(xintercept=11, color = "black",linetype='dashed',size=0.5)
    
    p=grid.arrange(p1,nrow = 1)
    dev.off()
    
    kk = kk+1
    all_subj_placenta_probe$sig[[kk]] = placenta_probe$sig
    all_subj_placenta_probe$type[[kk]] = placenta_probe$type
    all_subj_placenta_probe$ts[[kk]] = placenta_probe$ts
    all_subj_placenta_probe$probe[[kk]] = placenta_probe$probe
    all_subj_placenta_probe$std_dev_igm[[kk]] = placenta_probe$sd_igm
    all_subj_placenta_probe$std_dev_igg[[kk]] = placenta_probe$sd_igg
    all_subj_placenta_probe$probe_name[[kk]] = placenta_probe$pn
    all_subj_placenta_probe$gene_name[[kk]] = placenta_probe$gn
    
    cat("subject ", subj, " is done!!\n")
  }
  return(all_subj_placenta_probe)
}

plot_top_varying_placenta_related_signal <- function(max_signals,all_subj_placenta_probe,
                                                  iso_ID_name_arr,subj_time_info,topK){
  
  
  top_varying_pn_gn = NULL
  i = 0
  for (subj in subj_time_info$all_subj){
    
    i = i+1
    
    id_sorted_igm = sort(all_subj_placenta_probe$std_dev_igm[[i]],
                         decreasing = T, index.return=T)
    
    id_sorted_igg = sort(all_subj_placenta_probe$std_dev_igg[[i]],
                         decreasing = T, index.return=T)
    
    
    id_tmp_igm = id_sorted_igm$ix[1:topK]
    id_tmp_igg = id_sorted_igg$ix[1:topK]
    
    #id_tmp_igm = which(all_subj_placenta_probe$std_dev_igm[[i]] > 4.0) ##
    #id_tmp_igg = which(all_subj_placenta_probe$std_dev_igg[[i]] > 4.0) ## 
    
    cat('number of highly varying sig: igm: ', length(id_tmp_igm),
        ' igg: ',length(id_tmp_igg),'\n')
    
    top_varying_pn_gn$pn_igm[[i]] = all_subj_placenta_probe$probe_name[[i]][id_tmp_igm]
    top_varying_pn_gn$pn_igg[[i]] = all_subj_placenta_probe$probe_name[[i]][id_tmp_igg]
    
    top_varying_pn_gn$gn_igm[[i]] = all_subj_placenta_probe$gene_name[[i]][id_tmp_igm]
    top_varying_pn_gn$gn_igg[[i]] = all_subj_placenta_probe$gene_name[[i]][id_tmp_igg]
    
  }
  
  kk = 0
  for (subj in subj_time_info$all_subj){
    
    ts_id = which(grepl(subj,max_signals$col_names) == TRUE)
    ts = max_signals$col_names[ts_id]
    
    IGM_sig = max_signals$mat_IGM_F_by_B[,c(ts_id)]
    IGG_sig = max_signals$mat_IGG_F_by_B[,c(ts_id)]
    
    kk = kk+1
    
    L2_igm = length(top_varying_pn_gn$pn_igm[[kk]])
    L2_igg = length(top_varying_pn_gn$pn_igg[[kk]])
    
    placenta_probe_igm = NULL
    placenta_probe_igg = NULL
    indx_set_igm = vector()
    indx_set_igg = vector()
    
    tmp_L = length(ts_id)
    k = 1
    
    for(j in 1:L2_igm){
      indx_set_igm[j] = which(max_signals$row_names ==  top_varying_pn_gn$pn_igm[[kk]][j])
      placenta_probe_igm$sig[k:(k+tmp_L-1)] = IGM_sig[indx_set_igm[j],]
      placenta_probe_igm$ts[k:(k+tmp_L-1)] = ts
      placenta_probe_igm$probe[k:(k+tmp_L-1)] = top_varying_pn_gn$pn_igm[[kk]][j]   
      k = k+tmp_L
    }
    
    k = 1
    for(j in 1:L2_igg){
      indx_set_igg[j] = which(max_signals$row_names ==  top_varying_pn_gn$pn_igg[[kk]][j])
      placenta_probe_igg$sig[k:(k+tmp_L-1)] = IGG_sig[indx_set_igg[j],]
      placenta_probe_igg$ts[k:(k+tmp_L-1)] = ts
      placenta_probe_igg$probe[k:(k+tmp_L-1)] = top_varying_pn_gn$pn_igg[[kk]][j]  
      k = k+tmp_L
    }
    
    placenta_probe_plot_igm = as.data.frame(placenta_probe_igm)
    placenta_probe_plot_igm$ts = factor(placenta_probe_plot_igm$ts,levels=c(ts))
    
    placenta_probe_plot_igg = as.data.frame(placenta_probe_igg)
    placenta_probe_plot_igg$ts = factor(placenta_probe_plot_igg$ts,levels=c(ts))
    
    fname12 = sprintf("./%s_IGM_IGG_placenta_top_varying_sig_topK_%d.png", subj,topK)
    
    title1 = sprintf("Top Varying Placenta probe signal (IGM) across all time point for subject %s",subj)
    title2 = sprintf("Top Varying Placenta probe signal (IGG) across all time point for subject %s",subj)
    
    png(file=fname12,width=3500,height=1500,res=270)
    
    p1 = ggplot(placenta_probe_plot_igm,aes(x=ts,y=sig,group=probe)) + 
      geom_line()+
      geom_point()+
      scale_x_discrete(name='Timepoints')+
      scale_y_continuous(name='Signal value (in log 10)',
                         trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))+
      #coord_trans(y="log10")+
      ggtitle(title1) + theme(plot.title = element_text(hjust = 0.5))
    
    p2 = ggplot(placenta_probe_plot_igg,aes(x=ts,y=sig,group=probe)) + 
      geom_line()+
      geom_point()+
      scale_x_discrete(name='Timepoints')+
      scale_y_continuous(name='Signal value (in log 10)',
                         trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))+
      #coord_trans(y="log10")+
      ggtitle(title2) + theme(plot.title = element_text(hjust = 0.5))
    
    p=grid.arrange(p1,p2,nrow = 1)
    dev.off()
    
  }
  
  return(top_varying_pn_gn)   
}

#############
### paper related gene 

plot_paper_gene_related_signal <-function(max_signals,
                                          paper_GN_list,
                                          iso_ID_name_arr,
                                          subj_time_info,top){
    

  all_subj_paper_probe = list()
  kk = 0
  for(subj in subj_time_info$all_subj){
    
    ts_id = which(grepl(subj,max_signals$col_names) == TRUE)
    ts = max_signals$col_names[ts_id]
    
    IGM_sig = max_signals$mat_IGM_F_by_B[,c(ts_id)]
    IGG_sig = max_signals$mat_IGG_F_by_B[,c(ts_id)]
    
    L1 = length(paper_GN_list)
    
    probe_names = NULL
    
    for(i in 1:L1){
      id_tmp = which(iso_ID_name_arr$gene_name == paper_GN_list[i])
      probe_names = c(probe_names,iso_ID_name_arr$iso_id[id_tmp])
    }
    
    L2 = length(probe_names)
    indx_set = vector()
    
    paper_probe = NULL
    indx_set = vector()
    
    tmp_L = length(ts_id)
    k = 1
    
    for(i in 1:L2){
      indx_set[i] = which(max_signals$row_names == probe_names[i])
      
      paper_probe$sig[k:(k+tmp_L-1)] = IGM_sig[indx_set[i],]
      paper_probe$type[k:(k+tmp_L-1)] = rep('IGM',tmp_L)
      paper_probe$ts[k:(k+tmp_L-1)] = ts
      paper_probe$probe[k:(k+tmp_L-1)] = probe_names[i]
      paper_probe$gene_name[k:(k+tmp_L-1)] = iso_ID_name_arr$gene_name[which(iso_ID_name_arr$iso_id == probe_names[i])]
      paper_probe$sd_igm[i] = sd(paper_probe$sig[k:(k+tmp_L-1)])
      
      k = k+tmp_L
      paper_probe$sig[k:(k+tmp_L-1)] = IGG_sig[indx_set[i],]
      paper_probe$type[k:(k+tmp_L-1)] = rep('IGG',tmp_L)
      paper_probe$ts[k:(k+tmp_L-1)] = ts
      paper_probe$probe[k:(k+tmp_L-1)] = probe_names[i]
      paper_probe$gene_name[k:(k+tmp_L-1)] = iso_ID_name_arr$gene_name[which(iso_ID_name_arr$iso_id == probe_names[i])]
      paper_probe$sd_igg[i] = sd(paper_probe$sig[k:(k+tmp_L-1)])
      
      paper_probe$pn[i] = probe_names[i]
      id_tmp = which(iso_ID_name_arr$iso_id == paper_probe$pn[i])
      if(length(id_tmp) > 1){
        id_tmp = id_tmp[1]
      }
      paper_probe$gn[i] = iso_ID_name_arr$gene_name[id_tmp]
      
      k = k+tmp_L
      
    }
    
    paper_probe_plot = as.data.frame(paper_probe)
    paper_probe_plot$ts = factor(paper_probe_plot$ts,levels=c(ts))
    paper_probe_plot$type = factor(paper_probe_plot$type,levels=c('IGM','IGG'))
    
    
    
    if(top == 1){
      title = sprintf("Paper (all) probe signal across all time point for subject %s",subj)
      fname1 = sprintf("./%s_IGM_IGG_paper_sig_color_by_gn.png", subj)
      png(file=fname1,width=3500,height=1500,res=270)
      p1 = ggplot(paper_probe_plot,aes(x=ts,y=sig,group=probe,color=gene_name)) 
    }else{
      title = sprintf("Paper (top) probe signal across all time point for subject %s",subj)
      fname1 = sprintf("./%s_IGM_IGG_paper_sig.png", subj)
      png(file=fname1,width=3500,height=1500,res=270)
      p1 = ggplot(paper_probe_plot,aes(x=ts,y=sig,group=probe)) 
    }
    p1 = p1+ geom_line()+
      geom_point()+
      facet_grid(.~type)+
      scale_x_discrete(name='Timepoints')+
      scale_y_continuous(name='Signal value (in log 10)',
                         trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(10^.x)))+
      #scale_y_continuous(name='Signal value')+
      #coord_trans(y="log10")+
      ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    #scale_color_discrete("Samples")+
    #geom_vline(xintercept=11, color = "black",linetype='dashed',size=0.5)
    
    p=grid.arrange(p1,nrow = 1)
    dev.off()
    
    kk = kk+1
    all_subj_paper_probe$sig[[kk]] = paper_probe$sig
    all_subj_paper_probe$type[[kk]] = paper_probe$type
    all_subj_paper_probe$ts[[kk]] = paper_probe$ts
    all_subj_paper_probe$probe[[kk]] = paper_probe$probe
    all_subj_paper_probe$std_dev_igm[[kk]] = paper_probe$sd_igm
    all_subj_paper_probe$std_dev_igg[[kk]] = paper_probe$sd_igg
    all_subj_paper_probe$probe_name[[kk]] = paper_probe$pn
    all_subj_paper_probe$gene_name[[kk]] = paper_probe$gn
    
    cat("subject ", subj, " is done!!\n")
  }
  return(all_subj_paper_probe)
}

##############

plot_correlation <- function(IGM_signals,IGG_signals,selected){
  
  pals=colorRampPalette(brewer.pal(6, "OrRd"))
  col = pals(20)
  #col =  colorRampPalette(c("blue", "white", "red"))(20)
  
  if(selected == 0){
    f_corr = sprintf("./Corr_IGM.png")
    f_heatmap = sprintf("./Corr_IGM_heatmap.png")
  }else{
    f_corr = sprintf("./Corr_IGM_sel.png")
    f_heatmap = sprintf("./Corr_IGM_heatmap_sel.png")
  }
  res_IGM = cor(IGM_signals)
  png(f_corr)
  corrplot(res_IGM, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
  dev.off()
  png(f_heatmap)
  heatmap(x = res_IGM, col = col, symm = TRUE)
  dev.off()
  
  if(selected == 0){
    f_corr = sprintf("./Corr_IGG.png")
    f_heatmap = sprintf("./Corr_IGG_heatmap.png")
  }else{
    f_corr = sprintf("./Corr_IGG_sel.png")
    f_heatmap = sprintf("./Corr_IGG_heatmap_sel.png")
  }
  res_IGG = cor(IGG_signals)
  png(f_corr)
  corrplot(res_IGG, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
  dev.off()
  png(f_heatmap)
  heatmap(x = res_IGG, col = col, symm = TRUE)
  dev.off()

}

plot_rank_correlation <- function(IGM_signals,IGG_signals,selected){
  
  pals=colorRampPalette(brewer.pal(6, "PuOr"))
  col = pals(20)
  #col =  colorRampPalette(c("blue", "white", "red"))(20)
  
  if(selected == 0){
    f_corr = sprintf("././rank_corr_IGM.png")
    f_heatmap = sprintf("./rank_corr_IGM_heatmap.png")
  }else{
    f_corr = sprintf("././rank_corr_IGM_sel.png")
    #f_heatmap = sprintf("./rank_corr_IGM_heatmap_sel.png")
    f_heatmap = sprintf("./rank_corr_IGM_heatmap_sel.pdf")
  }
  res_IGM = cor(IGM_signals,method = "spearman")
  png(f_corr,width=5000,height=5000,res=750)
  corrplot(res_IGM, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
  dev.off()
  #png(f_heatmap,width=5000,height=5000,res=750)
  pdf(f_heatmap,width=8,height=8)
  heatmap.2(x = res_IGM, col = col, symm = TRUE,density.info="none",trace="none")
  dev.off()
  
  
  if(selected == 0){
    f_corr = sprintf("././rank_corr_IGG.png")
    f_heatmap = sprintf("./rank_corr_IGG_heatmap.png")
  }else{
    f_corr = sprintf("././rank_corr_IGG_sel.png")
    #f_heatmap = sprintf("./rank_corr_IGG_heatmap_sel.png")
    f_heatmap = sprintf("./rank_corr_IGG_heatmap_sel.pdf")
  }
  res_IGG = cor(IGG_signals,method = "spearman")
  png(f_corr,width=5000,height=5000,res=750)
  corrplot(res_IGG, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
  dev.off()
  #png(f_heatmap,width=5000,height=5000,res=750)
  pdf(f_heatmap,width=8,height=8)
  heatmap.2(x = res_IGG, col = col, symm = TRUE,density.info="none",trace="none")
  dev.off()
  
}

ttest_paired <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y, paired = TRUE,alternative = "two.sided")
  results$p.value
}

wilcox_test_paired <- function(df, grp1, grp2, type) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  if(type == 0){
    results = wilcox.test(x, y, paired = TRUE,alternative = "two.sided")
  }
  if(type == 1){
    results = wilcox.test(x, y, paired = TRUE,alternative = "greater")
  }
  if(type == 2){
    results = wilcox.test(x, y, paired = TRUE,alternative = "less")
  }
  results$p.value
}

ttest_paired_greater <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y, paired = TRUE,alternative = "greater")
  results$p.value
}

ttest_paired_less <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y, paired = TRUE,alternative = "less")
  results$p.value
}

ttest_unpaired <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y, paired = FALSE, alternative = "two.sided")
  results$p.value
}


find_pre_post_case_ctrl_id <- function(pre_post_case_ctrl_samples_names,sample_names){

  k = unique(unlist(lapply(pre_post_case_ctrl_samples_names,function(x) strsplit(x,'_')[[1]][2])))
  
  subj_id = as.numeric(unlist(substr(lapply(pre_post_case_ctrl_samples_names,function(x) strsplit(x,'_')[[1]][1]),2,2)))
  
  case_sample_names = pre_post_case_ctrl_samples_names[subj_id <=3 ]
  ctrl_sample_names = pre_post_case_ctrl_samples_names[subj_id > 3 ]
  
  str1 = sprintf('_%s_',k[1])
  str2 = sprintf('_%s_',k[2])
  
  pre_sample_names = pre_post_case_ctrl_samples_names[which(grepl(str1,pre_post_case_ctrl_samples_names) == TRUE)]
  post_sample_names = pre_post_case_ctrl_samples_names[which(grepl(str2,pre_post_case_ctrl_samples_names) == TRUE)]
  
  case_sample_names = pre_post_case_ctrl_samples_names
  ctrl_sample_names = pre_post_case_ctrl_samples_names  
  
  grp_pre_id = unlist(lapply(pre_sample_names,function(x) which(grepl(x,sample_names) == TRUE)))
  grp_post_id = unlist(lapply(post_sample_names,function(x) which(grepl(x,sample_names) == TRUE)))
  
  grp_case_id = unlist(lapply(case_sample_names,function(x) which(grepl(x,sample_names) == TRUE)))
  grp_ctrl_id = unlist(lapply(ctrl_sample_names,function(x) which(grepl(x,sample_names) == TRUE)))
  
  
  id_set = NULL
  id_set$grp_pre_id = grp_pre_id
  id_set$grp_post_id = grp_post_id
  id_set$grp_case_id = grp_case_id
  id_set$grp_ctrl_id = grp_ctrl_id
  
  return(id_set)
  
}

find_case_ctrl_pre_post_id <- function(pre_post_samples_names,sample_names){
  
  pre_sample_names = pre_post_samples_names[which(grepl('_m_',pre_post_samples_names) == TRUE)]
  post_sample_names = pre_post_samples_names[which(grepl('_p_',pre_post_samples_names) == TRUE)]
  
  sample_id_fixed_pre = unlist(lapply(pre_sample_names,function(x) which(grepl(x,sample_names) == TRUE)))
  sample_id_fixed_post = unlist(lapply(post_sample_names,function(x) which(grepl(x,sample_names) == TRUE)))
  
  case_grp_pre_id = sample_id_fixed_pre[1:3] ### S1,S2,S3 are converted
  case_grp_post_id = sample_id_fixed_post[1:3] ### S4,S5,s6 are controls
  
  crtl_grp_pre_id = sample_id_fixed_pre[4:6]
  crtl_grp_post_id = sample_id_fixed_post[4:6]
  
  id_set = NULL
  id_set$case_grp_pre_id = case_grp_pre_id
  id_set$crtl_grp_pre_id = crtl_grp_pre_id
  
  id_set$case_grp_post_id = case_grp_post_id
  id_set$crtl_grp_post_id = crtl_grp_post_id
  
  return(id_set)
  
}


run_paired_t_test_case_ctrl <- function(probe_set,max_signals,iso_ID_name_arr,id_set,type,threshold){
  
  if(type == 'IGM'){
    F_by_B_max = max_signals$mat_IGM_F_by_B
    selected_probes = probe_set$IGM_probe_set_union
  }
  if(type == 'IGG'){
    F_by_B_max = max_signals$mat_IGG_F_by_B
    selected_probes = probe_set$IGG_probe_set_union
  }
  
  
  F_by_B_max_case_3_paired_sample = 
    F_by_B_max[selected_probes,c(id_set$case_grp_pre_id,id_set$case_grp_post_id)]
  
  F_by_B_max_crtl_3_paired_sample = 
    F_by_B_max[selected_probes,c(id_set$crtl_grp_pre_id,id_set$crtl_grp_post_id)]
  
  raw_p_val_3_sample_case_paired_t_test = apply(
    as.data.frame(F_by_B_max_case_3_paired_sample),1,ttest_paired, 
    grp1 = c(1:3), grp2 = c(4:6))
  
  raw_p_val_3_sample_crtl_paired_t_test = apply(
    as.data.frame(F_by_B_max_crtl_3_paired_sample),1,ttest_paired, 
    grp1 = c(1:3), grp2 = c(4:6))

  ### case  
  raw_p_val_3_sample_case_paired_t_test_sorted = sort(raw_p_val_3_sample_case_paired_t_test,index.return = TRUE)
  id_case_thresholded = which(raw_p_val_3_sample_case_paired_t_test_sorted$x < threshold)
  raw_p_val_3_sample_case_paired_t_test_sorted_thresholded = raw_p_val_3_sample_case_paired_t_test_sorted$x[id_case_thresholded]
  gene_name_case_thresholded = unlist(lapply(names(id_case_thresholded),
                                        function(x) iso_ID_name_arr$gene_name[which(iso_ID_name_arr$iso_id == x)]))
  gene_name_case_thresholded = unique(gene_name_case_thresholded)
  gene_name_case_thresholded = unlist(lapply(gene_name_case_thresholded,function(x){strsplit(x,split="[-|.|_]")[[1]][1]}))
  
  fName = sprintf("gene_name_%s_case_%f.txt",type,threshold)
  write.table(gene_name_case_thresholded, file = fName, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  ###Entrez Gene ID
  EID_case_genes = rep(-1,length(gene_name_case_thresholded))
  for (k in 1:length(gene_name_case_thresholded)){
    EID_case_genes[k] = unlist(mget(x=gene_name_case_thresholded[k],envir=org.Hs.egALIAS2EG,ifnotfound=NA))
  }
  EID_case_genes_rem_NA = EID_case_genes[!is.na(EID_case_genes)]
  fName = sprintf("EID_case_genes_%s_%f.txt",type,threshold)
  write.table(EID_case_genes_rem_NA, file = fName, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  
  ### control
  raw_p_val_3_sample_crtl_paired_t_test_sorted = sort(raw_p_val_3_sample_crtl_paired_t_test,index.return = TRUE)
  id_crtl_thresholded = which(raw_p_val_3_sample_crtl_paired_t_test_sorted$x < threshold)
  raw_p_val_3_sample_crtl_paired_t_test_sorted_thresholded = raw_p_val_3_sample_crtl_paired_t_test_sorted$x[id_crtl_thresholded]
  gene_name_crtl_thresholded = unlist(lapply(names(id_crtl_thresholded),
                                             function(x) iso_ID_name_arr$gene_name[which(iso_ID_name_arr$iso_id == x)]))
  gene_name_crtl_thresholded = unique(gene_name_crtl_thresholded)
  
  gene_name_crtl_thresholded = unlist(lapply(gene_name_crtl_thresholded,function(x){strsplit(x,split="[-|.|_]")[[1]][1]}))
  
  fName = sprintf("gene_name_%s_crtl_%f.txt",type,threshold)
  write.table(gene_name_crtl_thresholded, file = fName, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  ###Entrez Gene ID
  EID_crtl_genes = rep(-1,length(gene_name_crtl_thresholded))
  for (k in 1:length(gene_name_crtl_thresholded)){
    EID_crtl_genes[k] = unlist(mget(x=gene_name_crtl_thresholded[k],envir=org.Hs.egALIAS2EG,ifnotfound=NA))
  }
  EID_crtl_genes_rem_NA = EID_crtl_genes[!is.na(EID_crtl_genes)]
  fName = sprintf("EID_crtl_genes_%s_%f.txt",type,threshold)
  write.table(EID_crtl_genes_rem_NA, file = fName, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}

#### n-pre vs n-post
run_paired_t_test_pre_post <- function(analysis_name,probe_set,max_signals,
                                       iso_ID_name_arr,id_set,n,type,threshold,
                                       test_type,dir_name){
  
  n = as.integer(n)
  
  if(type == 'IGM'){
    F_by_B_max = max_signals$mat_IGM_F_by_B
    selected_probes = probe_set$IGM_probe_set_union
  }
  if(type == 'IGG'){
    F_by_B_max = max_signals$mat_IGG_F_by_B
    selected_probes = probe_set$IGG_probe_set_union
  }
  
  grp_id_1 = c(1:n)
  grp_id_2 = c((n+1):(2*n))
  
  F_by_B_max_pre_post_n_paired_sample = F_by_B_max[selected_probes,c(id_set$grp_pre_id,id_set$grp_post_id)]

  if(test_type == 0){
   raw_p_val_n_sample_pre_post_paired_t_test = apply(
    as.data.frame(F_by_B_max_pre_post_n_paired_sample),1,ttest_paired, 
      grp1 = grp_id_1, grp2 = grp_id_2)
   fName1 = sprintf("./%s/%s_gene_name_%s_pre_post_%f_two_sided_n_%d.txt",dir_name,analysis_name,type,threshold,n)
   fName2 = sprintf("./%s/%s_EID_pre_post_genes_%s_%f_two_sided_n_%d.txt",dir_name,analysis_name,type,threshold,n)
   
  }
  if(test_type == 1){
    raw_p_val_n_sample_pre_post_paired_t_test = apply(
     as.data.frame(F_by_B_max_pre_post_n_paired_sample),1,ttest_paired_greater, 
      grp1 = grp_id_1, grp2 = grp_id_2)
    fName1 = sprintf("./%s/%s_gene_name_%s_pre_post_%f_pre_greater_post_n_%d.txt",dir_name,analysis_name,type,threshold,n)
    fName2 = sprintf("./%s/%s_EID_pre_post_genes_%s_%f_pre_greater_post_n_%d.txt",dir_name,analysis_name,type,threshold,n)
    
  }
  if(test_type == 2){
    raw_p_val_n_sample_pre_post_paired_t_test = apply(
      as.data.frame(F_by_B_max_pre_post_n_paired_sample),1,ttest_paired_less, 
      grp1 = grp_id_1,grp2 = grp_id_2)
    fName1 = sprintf("./%s/%s_gene_name_%s_pre_post_%f_pre_less_post_n_%d.txt",dir_name,analysis_name,type,threshold,n)
    fName2 = sprintf("./%s/%s_EID_pre_post_genes_%s_%f_pre_less_post_n_%d.txt",dir_name,analysis_name,type,threshold,n)
  }
  
  #### gene name 
  raw_p_val_n_sample_pre_post_paired_t_test_sorted = sort(raw_p_val_n_sample_pre_post_paired_t_test,index.return = TRUE)
  id_pre_post_thresholded = which(raw_p_val_n_sample_pre_post_paired_t_test_sorted$x < threshold)
  raw_p_val_n_sample_pre_post_paired_t_test_sorted_thresholded = raw_p_val_n_sample_pre_post_paired_t_test_sorted$x[id_pre_post_thresholded]
  
  gene_name_pre_post_thresholded = unlist(lapply(names(id_pre_post_thresholded),
                                             function(x) iso_ID_name_arr$gene_name[which(iso_ID_name_arr$iso_id == x)]))
  gene_name_pre_post_thresholded = unique(gene_name_pre_post_thresholded)
  gene_name_pre_post_thresholded = unlist(lapply(gene_name_pre_post_thresholded,function(x){strsplit(x,split="[-|.|_]")[[1]][1]}))
  
  write.table(gene_name_pre_post_thresholded, file = fName1, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  ###Entrez Gene ID
  EID_pre_post_genes = rep(-1,length(gene_name_pre_post_thresholded))
  for (k in 1:length(gene_name_pre_post_thresholded)){
    EID_pre_post_genes[k] = unlist(mget(x=gene_name_pre_post_thresholded[k],envir=org.Hs.egALIAS2EG,ifnotfound=NA))
  }
  EID_pre_post_genes_rem_NA = EID_pre_post_genes[!is.na(EID_pre_post_genes)]
  write.table(EID_pre_post_genes_rem_NA, file = fName2, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}


find_case_crtl_sym_diff <-function(gene_name_case,gene_name_crtl,type,threshold){
  
  gene_name_case_crtl_sym_diff = sym_diff(gene_name_case,gene_name_crtl)
  fName = sprintf("gene_name_case_crtl_sym_diff_%s_%f.txt",type,threshold)
  write.table(gene_name_case_crtl_sym_diff, file = fName, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  EID_case_crtl_sym_diff = rep(-1,length(gene_name_case_crtl_sym_diff))
  for (k in 1:length(EID_case_crtl_sym_diff)){
    EID_case_crtl_sym_diff[k] = unlist(mget(x=gene_name_case_crtl_sym_diff[k],envir=org.Hs.egALIAS2EG,ifnotfound=NA))
  }
  
  EID_case_crtl_sym_diff_rem_NA = EID_case_crtl_sym_diff[!is.na(EID_case_crtl_sym_diff)]
  fName = sprintf("EID_case_crtl_sym_diff_%s_%f.txt",type,threshold)
  write.table(EID_case_crtl_sym_diff_rem_NA, file = fName, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
}

write_background_gene_EID <-function(iso_ID_name_arr){
  
  background_all_gene_name = unique(iso_ID_name_arr$gene_name)
  background_all_gene_name = unlist(lapply(background_all_gene_name,function(x){strsplit(x,split="[-|.|_]")[[1]][1]}))
  write.table(background_all_gene_name, file = "background_all_gene_name.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  EID_background_all_gene = rep(-1,length(background_all_gene_name))
  for (k in 1:length(EID_background_all_gene)){
    EID_background_all_gene[k] = unlist(mget(x=background_all_gene_name[k],envir=org.Hs.egALIAS2EG,ifnotfound=NA))
  }
  EID_background_all_gene_rem_NA = EID_background_all_gene[!is.na(EID_background_all_gene)]
  write.table(EID_background_all_gene_rem_NA, file = "EID_background_all_gene.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

}


#### S1 
#### TBD
run_unpaired_t_test_sample_based_pre_post <- function(){
  
}
