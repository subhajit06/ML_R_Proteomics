#remove(list = ls())

#plot_placenta <-function(sel_type){
  
sel_type = 'dup1'
setwd("~/Desktop/ToDo/Pablo/plasma_sample_data/pipeline")

source("~/Desktop/ToDo/Pablo/plasma_sample_data/code/GPR_data_input.R")

load_libs()

#IGM_IGG_med = read_GPR_files()
load('./IGM_IGG_med.Rdata')

#sample_F_by_B_signals = signal_extract_from_GPR_files(IGM_IGG_med)
load('./sample_F_by_B_signals.RData')

n_total_probes = length(sample_F_by_B_signals$row_names)
n_samples =  length(sample_F_by_B_signals$col_names)

iso_ID_name_arr = map_isoform_ID_gene_name(IGM_IGG_med) 
#write_background_gene_EID(iso_ID_name_arr)

#QN_signals = gen_QN_signals(sample_F_by_B_signals)
max_signals = gen_max_signals(sample_F_by_B_signals,sel_type)
probe_names = max_signals$row_names
sample_names = max_signals$col_names

all_subjects_sample_cnt = as.vector(table(unlist(lapply(sample_names,function(x) 
  sprintf('%s',strsplit(x,'_')[[1]][1])))))


subj_time_info = subject_and_time(sample_names)

sd_cutoff_factor = 10 ### I am not sure what is a good cut-off value
probe_set = get_threshold_probes(max_signals,sd_cutoff_factor)

IGM_signals = max_signals$mat_IGM_F_by_B
IGG_signals = max_signals$mat_IGG_F_by_B

IGM_signals_sel = max_signals$mat_IGM_F_by_B[probe_set$IGM_probe_set_union,]
IGG_signals_sel = max_signals$mat_IGG_F_by_B[probe_set$IGG_probe_set_union,]


overall_sig = overall_immune_signal(IGM_signals_sel,IGG_signals_sel,sample_names)
str1 = sprintf('all_%s',sel_type)
plot_overall_signal(overall_sig,str1)

##################################
######## Placenta related ########

EID_file_obj = read.table("~/Desktop/ToDo/Pablo/plasma_sample_data/pipeline/david_output/pre_post_IGM_0.05_two_sided_EID_list_placenta.txt",sep="\t")
L = dim(EID_file_obj)[1]

set1 = NULL
for(i in 2:L){
  set1 = union(set1,unlist(strsplit(
    as.vector(EID_file_obj$V6[[i]]),split=', ')))
}

placenta_GN_list = as.vector(unlist(mget(x=set1,envir=org.Hs.egSYMBOL,ifnotfound=NA)))

all_subj_placenta_probe = plot_placenta_related_signal(max_signals,placenta_GN_list,
                                                       iso_ID_name_arr,subj_time_info)

############################
### plot overall placenta signal 
cnt_subj = as.vector(table(overall_sig$subj))
n_placenta_probe = length(all_subj_placenta_probe$sig[[1]])/(cnt_subj[1]*2) ## 2 for IGM and IGG

placenta_overall = NULL
placenta_overall$subj = overall_sig$subj
placenta_overall$times = overall_sig$times

n_subj = length(unique(overall_sig$subj))
id_1 = c(1,cumsum(cnt_subj)+1)[1:n_subj]
id_2 = cumsum(cnt_subj)

placenta_all_probe_all_ts = NULL

for(i in 1:length(subj_time_info$all_subj)){
  m_tmp = NULL
  s_tmp = NULL
  beg_indx_IGM = seq(1,n_placenta_probe*cnt_subj[i]*2,by=cnt_subj[i]*2)
  end_indx_IGM = seq(cnt_subj[i],n_placenta_probe*cnt_subj[i]*2,by=cnt_subj[i]*2)
  
  for(j in 1:length(beg_indx_IGM)){
    s_tmp = c(s_tmp,all_subj_placenta_probe$sig[[i]][beg_indx_IGM[j]:end_indx_IGM[j]])
  }
  
  m_tmp = matrix(s_tmp,nrow = length(beg_indx_IGM))
  rownames(m_tmp) = all_subj_placenta_probe$probe_name[[i]]
  colnames(m_tmp) = all_subj_placenta_probe$ts[[i]][1:cnt_subj[i]]
  
  placenta_all_probe_all_ts$IGM_mat[[i]] = m_tmp
  v_mean_IGM = apply(m_tmp,2,mean)
  placenta_overall$IGM[id_1[i]:id_2[i]] = v_mean_IGM
  
  m_tmp = NULL
  s_tmp = NULL
  beg_indx_IGG = seq(cnt_subj[i]+1,n_placenta_probe*cnt_subj[i]*2,by=cnt_subj[i]*2)
  end_indx_IGG = seq(cnt_subj[i]*2,n_placenta_probe*cnt_subj[i]*2,by=cnt_subj[i]*2)
  
  for(j in 1:length(beg_indx_IGG)){
    s_tmp = c(s_tmp,all_subj_placenta_probe$sig[[i]][beg_indx_IGG[j]:end_indx_IGG[j]])
  }
  m_tmp = matrix(s_tmp,nrow = length(beg_indx_IGG))
  v_mean_IGG = apply(m_tmp,2,mean)
  placenta_overall$IGG[id_1[i]:id_2[i]] = v_mean_IGG
  
}

placenta_overall = as.data.frame(placenta_overall)
str1 = sprintf('placenta_%s',sel_type)
plot_overall_signal(placenta_overall,str1)
###############


#### 6 different plot
#for(i in 1:length(subj_time_info$all_subj)){
  df1 = NULL
  df1 = as.data.frame(placenta_all_probe_all_ts$IGM_mat[[1]])
  df1$probe_id =  all_subj_placenta_probe$probe_name[[1]]
  df_melted_1 = melt(df1,id.vars = 'probe_id')
  names(df_melted_1)[2] = 'timepoint'
  p1 = ggplot(data=df_melted_1, mapping=aes(x = timepoint, y = value, group = probe_id)) + geom_line(color = '#e41a1c')
  
  df2 = NULL
  df2 = as.data.frame(placenta_all_probe_all_ts$IGM_mat[[2]])
  df2$probe_id =  all_subj_placenta_probe$probe_name[[2]]
  df_melted_2 = melt(df2,id.vars = 'probe_id')
  names(df_melted_2)[2] = 'timepoint'
  p2 = ggplot(data=df_melted_2, mapping=aes(x = timepoint, y = value, group = probe_id)) + geom_line(color = '#377eb8')
  
  df3 = NULL
  df3 = as.data.frame(placenta_all_probe_all_ts$IGM_mat[[3]])
  df3$probe_id =  all_subj_placenta_probe$probe_name[[3]]
  df_melted_3 = melt(df3,id.vars = 'probe_id')
  names(df_melted_3)[2] = 'timepoint'
  p3 = ggplot(data=df_melted_3, mapping=aes(x = timepoint, y = value, group = probe_id)) + geom_line(color = '#4daf4a')
  
  df4 = NULL
  df4 = as.data.frame(placenta_all_probe_all_ts$IGM_mat[[4]])
  df4$probe_id =  all_subj_placenta_probe$probe_name[[4]]
  df_melted_4 = melt(df4,id.vars = 'probe_id')
  names(df_melted_4)[2] = 'timepoint'
  p4 = ggplot(data=df_melted_4, mapping=aes(x = timepoint, y = value, group = probe_id)) + geom_line(color = '#984ea3')
  
  df5 = NULL
  df5 = as.data.frame(placenta_all_probe_all_ts$IGM_mat[[5]])
  df5$probe_id =  all_subj_placenta_probe$probe_name[[5]]
  df_melted_5 = melt(df5,id.vars = 'probe_id')
  names(df_melted_5)[2] = 'timepoint'
  p5 = ggplot(data=df_melted_5, mapping=aes(x = timepoint, y = value, group = probe_id)) + geom_line(color = '#ff7f00')
  
  df6 = NULL
  df6 = as.data.frame(placenta_all_probe_all_ts$IGM_mat[[6]])
  df6$probe_id =  all_subj_placenta_probe$probe_name[[6]]
  df_melted_6 = melt(df6,id.vars = 'probe_id')
  names(df_melted_6)[2] = 'timepoint'
  p6 = ggplot(data=df_melted_6, mapping=aes(x = timepoint, y = value, group = probe_id)) + geom_line(color = '#a65628')
  
  #fName = sprintf("./S%d_IGM_placenta_probe_profile.png",i)
  fName = sprintf("./S1_S6_IGM_placenta_probe_profile_%s.pdf",sel_type)
  #png(file=fName,width=3500,height=3500,res=270)
  pdf(file=fName,width=8,height=8)
  p=grid.arrange(p1,p2,p3,p4,p5,p6,nrow = 3)
  dev.off()
  
  df_all = NULL
  df_all$df1 = df1
  df_all$df2 = df2
  df_all$df3 = df3
  df_all$df4 = df4
  df_all$df5 = df5
  df_all$df6 = df6
  
  fName = sprintf("./df_all_%s.RData",sel_type)
  save(df_all,file=fName)
  
#}
}


plot_placenta('dup1')
plot_placenta('dup2')


load("./df_all_dup1.RData")

l1 = df_all$df4[4:9,'probe_id']
l2 = df_all$df5[4:9,'probe_id']
l3 = df_all$df6[4:9,'probe_id']

d1 = intersect(intersect(l1,l2),l3)
save(d1,file="./d1.RData")

load("./df_all_dup2.RData")
l1 = df_all$df4[4:9,'probe_id']
l2 = df_all$df5[4:9,'probe_id']
l3 = df_all$df6[4:9,'probe_id']

d2 = intersect(intersect(l1,l2),l3)
save(d2,file="./d2.RData")

d12 = intersect(d1,d2)
id1 = unlist(lapply(d12,function(x) {which(iso_ID_name_arr$iso_id == x)}))
unique(iso_ID_name_arr$gene_name[id1])

#load("./df_all_dup1.RData")
load("./df_all_dup2.RData")
sel_type = 'dup2'

probe_info_1 = NULL
p_id = 4
probe_info_1$IGM = c(as.vector(unlist(df_all$df1[p_id,1:8])), 
                     as.vector(unlist(df_all$df2[p_id,1:4])),
                     as.vector(unlist(df_all$df3[p_id,1:4])),
                     as.vector(unlist(df_all$df4[p_id,1:3])),
                     as.vector(unlist(df_all$df5[p_id,1:3])),
                     as.vector(unlist(df_all$df6[p_id,1:3])))

id1 = which(iso_ID_name_arr$iso_id == df_all$df1[p_id,'probe_id'])
gn1 = iso_ID_name_arr$gene_name[id1]
probe_info_1$IGG = rep(0,n_samples)
probe_info_1$subj = overall_sig$subj
probe_info_1$times = overall_sig$times
probe_info_1 = as.data.frame(probe_info_1)
str1 = sprintf("%s_%s_%s",gn1,sel_type,df_all$df1[p_id,'probe_id'])
plot_one_gene_signal(probe_info_1,str1)

probe_info_2 = NULL
p_id = 5
probe_info_2$IGM = c(as.vector(unlist(df_all$df1[p_id,1:8])), 
                     as.vector(unlist(df_all$df2[p_id,1:4])),
                     as.vector(unlist(df_all$df3[p_id,1:4])),
                     as.vector(unlist(df_all$df4[p_id,1:3])),
                     as.vector(unlist(df_all$df5[p_id,1:3])),
                     as.vector(unlist(df_all$df6[p_id,1:3])))

id1 = which(iso_ID_name_arr$iso_id == df_all$df1[p_id,'probe_id'])
gn1 = iso_ID_name_arr$gene_name[id1]
probe_info_2$IGG = rep(0,n_samples)
probe_info_2$subj = overall_sig$subj
probe_info_2$times = overall_sig$times
probe_info_2 = as.data.frame(probe_info_2)
str1 = sprintf("%s_%s_%s",gn1,sel_type,df_all$df1[p_id,'probe_id'])
plot_one_gene_signal(probe_info_2,str1)

probe_info_3 = NULL
p_id = 6
probe_info_3$IGM = c(as.vector(unlist(df_all$df1[p_id,1:8])), 
                     as.vector(unlist(df_all$df2[p_id,1:4])),
                     as.vector(unlist(df_all$df3[p_id,1:4])),
                     as.vector(unlist(df_all$df4[p_id,1:3])),
                     as.vector(unlist(df_all$df5[p_id,1:3])),
                     as.vector(unlist(df_all$df6[p_id,1:3])))

id1 = which(iso_ID_name_arr$iso_id == df_all$df1[p_id,'probe_id'])
gn1 = iso_ID_name_arr$gene_name[id1]
probe_info_3$IGG = rep(0,n_samples)
probe_info_3$subj = overall_sig$subj
probe_info_3$times = overall_sig$times
probe_info_3 = as.data.frame(probe_info_3)
str1 = sprintf("%s_%s_%s",gn1,sel_type,df_all$df1[p_id,'probe_id'])
plot_one_gene_signal(probe_info_3,str1)

probe_info_4 = NULL
p_id = 7
probe_info_4$IGM = c(as.vector(unlist(df_all$df1[p_id,1:8])), 
                     as.vector(unlist(df_all$df2[p_id,1:4])),
                     as.vector(unlist(df_all$df3[p_id,1:4])),
                     as.vector(unlist(df_all$df4[p_id,1:3])),
                     as.vector(unlist(df_all$df5[p_id,1:3])),
                     as.vector(unlist(df_all$df6[p_id,1:3])))

id1 = which(iso_ID_name_arr$iso_id == df_all$df1[p_id,'probe_id'])
gn1 = iso_ID_name_arr$gene_name[id1]
probe_info_4$IGG = rep(0,n_samples)
probe_info_4$subj = overall_sig$subj
probe_info_4$times = overall_sig$times
probe_info_4 = as.data.frame(probe_info_4)
str1 = sprintf("%s_%s_%s",gn1,sel_type,df_all$df1[p_id,'probe_id'])
plot_one_gene_signal(probe_info_4,str1)

probe_info_5 = NULL
p_id = 8
probe_info_5$IGM = c(as.vector(unlist(df_all$df1[p_id,1:8])), 
                     as.vector(unlist(df_all$df2[p_id,1:4])),
                     as.vector(unlist(df_all$df3[p_id,1:4])),
                     as.vector(unlist(df_all$df4[p_id,1:3])),
                     as.vector(unlist(df_all$df5[p_id,1:3])),
                     as.vector(unlist(df_all$df6[p_id,1:3])))

id1 = which(iso_ID_name_arr$iso_id == df_all$df1[p_id,'probe_id'])
gn1 = iso_ID_name_arr$gene_name[id1]
probe_info_5$IGG = rep(0,n_samples)
probe_info_5$subj = overall_sig$subj
probe_info_5$times = overall_sig$times
probe_info_5 = as.data.frame(probe_info_5)
str1 = sprintf("%s_%s_%s",gn1,sel_type,df_all$df1[p_id,'probe_id'])
plot_one_gene_signal(probe_info_5,str1)

probe_info_6 = NULL
p_id = 9
probe_info_6$IGM = c(as.vector(unlist(df_all$df1[p_id,1:8])), 
                     as.vector(unlist(df_all$df2[p_id,1:4])),
                     as.vector(unlist(df_all$df3[p_id,1:4])),
                     as.vector(unlist(df_all$df4[p_id,1:3])),
                     as.vector(unlist(df_all$df5[p_id,1:3])),
                     as.vector(unlist(df_all$df6[p_id,1:3])))

id1 = which(iso_ID_name_arr$iso_id == df_all$df1[p_id,'probe_id'])
gn1 = iso_ID_name_arr$gene_name[id1]
probe_info_6$IGG = rep(0,n_samples)
probe_info_6$subj = overall_sig$subj
probe_info_6$times = overall_sig$times
probe_info_6 = as.data.frame(probe_info_6)
str1 = sprintf("%s_%s_%s",gn1,sel_type,df_all$df1[p_id,'probe_id'])
plot_one_gene_signal(probe_info_6,str1)

