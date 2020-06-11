run_EN <- function(X,iso_ID_name_arr){
    
  sample_names = rownames(X)
  
  V1 = as.vector(unlist(lapply(sample_names, function(x) {ifelse(strsplit(x,"_")[[1]][2] == 'm',-1,1)}))) 
  V2 = as.vector(as.numeric(unlist(lapply(sample_names, function(x) {strsplit(x,"_")[[1]][3]}))))
  GA = V1*V2
  #X = cbind(X,GA)
  Z = factor(c(rep(1,16),rep(0,9)))
  Y = Z
  
  ### we have now data predictor X and response Y
  k = 0
  r = NULL
  e = NULL 
  list_genes = NULL
  
  for(a in seq(0,1.0,by=0.1)){
  
    a1 = a
    #fit1 = glmnet(X,Y,family='binomial',standardize = F,alpha=a1, nlambda = 200)
    fit1 = cv.glmnet(X,Y,nfolds = nrow(X), type.measure="class",
                      alpha = a1, grouped = F,family = "binomial")#,standardize = T)
    
    coeff1 = coef(fit1, s="lambda.min")
    ###https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
    features = NULL
    coeffs = NULL
    features = coeff1@Dimnames[[1]][ which(coeff1 != 0 ) ]  #intercept included
    coeffs   = coeff1              [ which(coeff1 != 0 ) ]  #intercept included
  
    res1 = NULL
    res1$features = features[-1]
    res1$coefs = coeffs[-1]
    res1$genes = unlist(lapply(res1$features,function(x){iso_ID_name_arr$gene_name[iso_ID_name_arr$iso_id == x]}))
  
    #png("./fit1_g.png")
    #plot(fit1)
    #dev.off()
    # k=dim(fit1$beta)[2]
    # v1 = as.vector((which(fit1$beta[,k] !=0)))
    # 
    # gn1 = unique(iso_ID_name_arr$gene_name[v1])
    # id1 = unique(iso_ID_name_arr$iso_id[v1])
  
    # fit2 = glmnet(X[1:16,],Y[1:16],family='binomial',standardize = F,alpha=a1,nlambda = 200)
    # #png("./fit2_g.png")
    # #plot(fit2)
    # #dev.off()
    # k=dim(fit2$beta)[2]
    # v2 = as.vector((which(fit2$beta[,k] !=0)))
    # gn2 = unique(iso_ID_name_arr$gene_name[v2])
    # 
    # fit3 = glmnet(X[17:25,],Y[17:25],family='binomial',standardize = F,alpha=a1,nlambda = 200)
    # #png("./fit3_g.png")
    # #plot(fit3)
    # #dev.off()
    # k=dim(fit3$beta)[2]
    # v3 = as.vector((which(fit3$beta[,k] !=0)))
    # gn3 = unique(iso_ID_name_arr$gene_name[v3])
    # 
    k = k+1
    #r[[k]] = (setdiff(gn1,gn3)) ###  "MTMR1"  "CALHM1" "CTAGE5"
    #r[[k]] = gn1 ###  "MTMR1"  "CALHM1" "CTAGE5"
    #e[[k]] = id1
    #cat(a,": ",r[[k]],"\t","length: ",length(gn1),"\n")
    #cat(a,": ",e[[k]],"\t","length: ",length(id1),"\n")
    ## "CALHM1" is imp
    cat(a,": length ", length(res1$features),"\n")
    #cat(res1$genes,"\n")
    list_genes[[k]] = res1$genes
  }
  
  return(list_genes)
  #rf=randomForest(X, Y, proximity=TRUE, importance=FALSE, norm.votes=TRUE) 
  #plot(fit1)

}
