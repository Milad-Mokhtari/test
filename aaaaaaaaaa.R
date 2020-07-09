library(rlist)
library(Rmosek)
#library(profvis)
#library(bigmemory)
library(Matrix)
#profvis({
source("group_l1_train")

DATA_PATH <- "./data"
RESULT_PATH <- "./results_inf"
cohorts <- c("ACC","CHOL","DLBC","KICH","LAML","LGG")
args <- commandArgs(trailingOnly = TRUE)
cohorts <- cohorts[as.numeric(args[[1]]) %% length(cohorts) + 1]
##R <- as.numeric(args[[2]])

##load(sprintf("%s/TCGA-%s.RData", DATA_PATH, toupper(cohorts)))

#cohorts <- c("TCGA-ACC","TCGA-CHOL","TCGA-DLBC","TCGA-ESCA","TCGA-KICH",
#            "TCGA-LAML","TCGA-MESO","TCGA-PAAD","TCGA-PCPG","TCGA-READ",
#             "TCGA-SKCM","TCGA-TGCT","TCGA-THYM","TCGA-UCS","T CGA-UVM")
#DATA_PATH <- "/home/midaslab/tcga/primary_data"
#RESULT_PATH <- "/home/mmokhtaridoost/Milad/results"
#load(sprintf("%s/%s.RData",DATA_PATH,cohorts))
###load("TCGA-CHOL.RData")
###cohorts <- c("TCGA-CHOL")
##DATA_PATH <- "/home/users/gonen/lustre/mirna_mrna_regulatory_l1-minimization/data"
##RESULT_PATH <- "/home/users/gonen/data_transfer/Milad/results"
##load(sprintf("%s/%s.RData",DATA_PATH,cohorts[as.numeric(args[[1]])]))
if (file.exists(sprintf("%s/L1_minimization_%s_summary.RData", RESULT_PATH, cohorts[as.numeric(args[[1]])])) == FALSE) {
  load(sprintf("%s/TCGA-%s.RData", DATA_PATH, toupper(cohorts)))
  
  
  set.seed(6676)
  
  
  common_patident<-intersect(
    rownames(TCGA$mirna),rownames(TCGA$mrna))
  X<-TCGA$mirna[common_patident,]
  
  Y<-TCGA$mrna[common_patident,]
  
  X<-X[,which(colMeans(X>0)>=0.5)]
  X<-(scale(log2(X+1)))
  
  Y<-Y[,which(colMeans(Y>0)>=0.5)]
  Y<-(scale(log2(Y+1)))
#####################################

  maxiter <- 10000
  tol <- 1e-06
  
  
#####################################
  
  l_set<-c(0.5,0.25,0.1,0.05,0.01,0.005,0.001)
  fold_count <- 4
  train_ratio <- 0.8
  
######################################

  lambda_matrix <- matrix(NA, nrow = fold_count, ncol = length(l_set), dimnames = list(1:fold_count, sprintf("%g", l_set)))
  
  train_sample_indices <- sample(Y, ceiling(train_ratio * length(Y)))
  sample_allocation <- sample(rep(1:fold_count, ceiling(length(train_negative_indices) / fold_count)), length(train_negative_indices))
  ####################################
  for (fold in 1:fold_count) {
    train_indices <- train_sample_indices[which(sample_allocation != fold)]
    test_indices <- train_sample_indices[which(sample_allocation == fold)]
    
    X_train <- X[train_indices,]
    X_test <- X[test_indices,]
    X_train <- scale(X_train)
    X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
    
    #N_train <- nrow(X_train)
    #N_test <- nrow(X_test)

    Y_train <- Y[train_indices]
    Y_test <- Y[test_indices]
  ########################
    R<- 10
    n<-nrow(X_train) # Number of data instances (n)
    t<-ncol(Y_train) # Number of response variables
    d <- ncol(X_train) # Number of predictor variables  
    dimension <- list()
    dimension$n <- n
    dimension$t <- t
    dimension$d <- d  
  ##########################
    for (lam in l_set) {
      print(sprintf("running fold = %d, lambda = %g", fold, lam))
      parameters <- list()
      parameters$lam <- lam
      
      state <- group_l1_train(dimension, X_train, Y_train, lam)
      Y_hat <- X_test %*% state$W 
      
      lambda_matrix[fold, sprintf("%g", lam)] <- mean(sqrt(colMeans((Y_hat-Y_test)^2))/apply(Y_test,2,sd)) 
                                                
    }
  }
  ########################
  R<- 10
  n<-nrow(X) # Number of data instances (n)
  t<-ncol(Y) # Number of response variables
  d <- ncol(X) # Number of predictor variables  
  dimension <- list()
  dimension$n <- n
  dimension$t <- t
  dimension$d <- d  
  ##########################
  l_star_NRMSE <- l_set[max.col(t(colMeans(lambda_matrix)), ties.method = "last")]    
  lam <- l_star_NRMSE
  state <- group_l1_train(dimension, X, Y, lam)
  W <- state$W
}