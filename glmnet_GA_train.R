run_EN <- function(X,iso_ID_name_arr){
  
  sample_names = rownames(X)
  
  V1 = as.vector(unlist(lapply(sample_names, function(x) {ifelse(strsplit(x,"_")[[1]][2] == 'm',-1,1)}))) 
  V2 = as.vector(as.numeric(unlist(lapply(sample_names, function(x) {strsplit(x,"_")[[1]][3]}))))
  GA = V1*V2
  #X = cbind(X,GA)
  #Z = c(rep(1,16),rep(0,9))
  Y = GA#as.matrix(cbind(Z,GA))
  #dimnames(Y)[[1]] = c(1:25)
  cat(Y,"\n")
  #X = X[id_sel,]
  #Y = Y[id_sel]
  
  X1 = X[1:16,]
  X2 = X[17:25,]
  
  Y1 = Y[1:16]
  Y2 = Y[17:25]
  
  cat('dim X: ', dim(X)," length Y: ", length(Y),'\n')
  ### we have now data predictor X and response Y
  k = 0
  r = NULL
  e = NULL 
  list_genes = NULL
  
  for(a in seq(0,1.0,by=0.1)){
    
    a1 = a
    #fit1 = glmnet(X,Y,family='binomial',standardize = F,alpha=a1, nlambda = 200)
    #fit1 = cv.glmnet(X,Y,family='gaussian',standardize = F,alpha=a1, nlambda = 200)
    fit1 = cv.glmnet(X1, Y1, nfolds = nrow(X1),alpha = a1, family = "gaussian", grouped=FALSE)
    fit2 = cv.glmnet(X2, Y2, nfolds = nrow(X2),alpha = a1, family = "gaussian", grouped=FALSE)
    
    #fit1 = cv.glmnet(X,Y,nfolds = nrow(X), type.measure="class", 
    #                 alpha = a1, grouped = F,family = "binomial",standardize = F)
    #fit1 = cv.glmnet(X,Y,nfolds = nrow(X), alpha = a1, family = "mgaussian")
      
    coeff1 = coef(fit1, s="lambda.min")
    ###https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
    features1 = NULL
    coeffs1 = NULL
    features1 = coeff1@Dimnames[[1]][ which(coeff1 != 0 ) ]  #intercept included
    coeffs1   = coeff1              [ which(coeff1 != 0 ) ]  #intercept included
    
    res1 = NULL
    res1$features = features1[-1]
    res1$coefs = coeffs1[-1]
    res1$genes = unlist(lapply(res1$features,function(x){iso_ID_name_arr$gene_name[iso_ID_name_arr$iso_id == x]}))
    
    #########
    
    coeff2 = coef(fit2, s="lambda.min")
    ###https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
    features2 = NULL
    coeffs2 = NULL
    features2 = coeff2@Dimnames[[1]][ which(coeff2 != 0 ) ]  #intercept included
    coeffs2   = coeff2              [ which(coeff2 != 0 ) ]  #intercept included
    
    res2 = NULL
    res2$features = features2[-1]
    res2$coefs = coeffs2[-1]
    res2$genes = unlist(lapply(res2$features,function(x){iso_ID_name_arr$gene_name[iso_ID_name_arr$iso_id == x]}))
    
    #####
    
    k = k+1
    
    sym_diff_gene = setdiff( union(res1$genes, res2$genes), intersect(res1$genes, res2$genes))
    #cat(a,": length ", length(res1$features),"\n")
    #cat(res1$genes,"\n")
    #list_genes[[k]] = res1$genes
    
    cat(a,": length ", length(sym_diff_gene),"\n")
    #cat(sym_diff_gene,"\n")
    list_genes[[k]] = sym_diff_gene
    
  }
  
  return(list_genes)
}


run_EN_1 <- function(X,iso_ID_name_arr,type=1){
  
  sample_names = rownames(X)
  
  V1 = as.vector(unlist(lapply(sample_names, function(x) {ifelse(strsplit(x,"_")[[1]][2] == 'm',-1,1)}))) 
  V2 = as.vector(as.numeric(unlist(lapply(sample_names, function(x) {strsplit(x,"_")[[1]][3]}))))
  GA = V1*V2
  #X = cbind(X,GA)
  Z = c(rep(1,16),rep(0,9))
  Y = Z#as.matrix(cbind(Z,GA))
  #dimnames(Y)[[1]] = c(1:25)
  cat(Y,"\n")
  #X = X[id_sel,]
  #Y = Y[id_sel]
  if(type == 1){
    w_vec = rep(1,25)
  }
  if(type == 2){
    w_vec = c(rep(1/8,8),rep(1/4,8),rep(1/3,9))
  }
  cat('dim X: ', dim(X)," length Y: ", length(Y),'\n')
  ### we have now data predictor X and response Y
  k = 0
  r = NULL
  e = NULL 
  list_genes = NULL
  
  for(a in seq(0,1.0,by=0.1)){
    
    a1 = a
    #fit1 = glmnet(X,Y,family='binomial',standardize = F,alpha=a1, nlambda = 200)
    #fit1 = cv.glmnet(X,Y,family='gaussian',standardize = F,alpha=a1, nlambda = 200)
    #fit1 = cv.glmnet(X1, Y1, nfolds = nrow(X1),alpha = a1, family = "gaussian", grouped=FALSE)
    
    fit1 = cv.glmnet(X,Y,nfolds = nrow(X), type.measure="class",weights = w_vec,
                     alpha = a1, grouped = F,family = "binomial")
    
    #fit1 = cv.glmnet(X,Y,nfolds = nrow(X), type.measure="class", 
    #                 alpha = a1, grouped = F,family = "binomial",standardize = F)
    #fit1 = cv.glmnet(X,Y,nfolds = nrow(X), alpha = a1, family = "mgaussian")
    
    coeff1 = coef(fit1, s="lambda.min")
    ###https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
    features1 = NULL
    coeffs1 = NULL
    features1 = coeff1@Dimnames[[1]][ which(coeff1 != 0 ) ]  #intercept included
    coeffs1   = coeff1              [ which(coeff1 != 0 ) ]  #intercept included
    
    res1 = NULL
    res1$features = features1[-1]
    res1$coefs = coeffs1[-1]
    res1$genes = unlist(lapply(res1$features,function(x){iso_ID_name_arr$gene_name[iso_ID_name_arr$iso_id == x]}))
    
    #########
    k = k+1
    
    cat(a,": length ", length(res1$features),"\n")
    #cat(res1$genes,"\n")
    list_genes[[k]] = res1$genes

  }
  
  return(list_genes)
}
