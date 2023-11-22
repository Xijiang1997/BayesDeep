# ***README***
# The following script is used to run the BayesDeep method for spatial transcriptomics data
# proposed in the submitted manuscript titled "Reconstructing Spatial Transcriptomics at the Single-cell Resolution with BayesDeep".
# ***END***

# load packages
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(Matrix)

# load cpp functions
Rcpp::sourceCpp('R/bayesdeep_cpp_function_vs.cpp')

# Generate spot-cell spatial relationship matrix
# T_spot: a matrix with x and y coordinates of spots
# T_cell: a matrix with x and y coordinates of all identified cells
# spot_diameter: diameter of the barcoded area of spots for ST data
get.spatial.matrix <- function(T_spot, T_cell, spot_diameter){
  T_spot <- as.matrix(T_spot)
  T_cell <- as.matrix(T_cell)
  n_cell <- dim(T_cell)[1]
  n_spot <- dim(T_spot)[1]
  
  G <- Matrix(0, nrow = n_spot, ncol = n_cell, sparse = TRUE)
  
  cc <- 0
  # assign cells to spots
  for (i in 1:n_cell){
    dist_matrix <- vectorized_pdist(matrix(T_cell[i,], nrow = 1), as.matrix(T_spot))
    min_dist_cur <- min(dist_matrix)
    if (min_dist_cur <= spot_diameter/2){
      G[which(dist_matrix == min_dist_cur), i] <- 1
    }
    if (floor(i*100/n_cell) == cc)
    {
      print(paste0(cc, '% has been done'))
      cc <- cc + 10
    }
  }
  return(G)
}


# Generate dummy variables for cell type
# Z: a vector with categorical feature
# class_name: categories of vector Z; will be the column names of dummy covariate matrix
get.dummy <- function(Z, class_name = NULL){
  if (is.null(class_name)){
  class_name <- unique(Z)
  }
  Z_dummy <- matrix(0, ncol = length(class_name), nrow = length(Z))
  colnames(Z_dummy) <- class_name
  for (i in 1:length(class_name)){
    Z_dummy[Z == class_name[i], i] <- 1
  }
  return(Z_dummy)
}

# Get the estimated size factor
# count: the n*p matrix for gene expression counts, where n is the number of spots and p is the number of genes
# norm_method: the method to calculate size factor. Choices include 'tss', 'q75', and 'rle'. For details, please refer paper https://onlinelibrary.wiley.com/doi/10.1002/sim.9530
get.size.factor <- function(count, norm_method = 'tss'){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  
  count_rowsum <- rowSums(count)
  
  if(norm_method == "tss")
  {
    ## TSS(Total Sum Scaling)
    ### scale-factors
    raw_s_factors <- rowSums(count)/mean(rowSums(count))
  }
  else if(norm_method == "q75")
  {
    ## Q75(Upper Quantile normalization)
    ### scale factors
    count_q75 <- apply(count, 1,function(x){quantile(x[x>0],0.75)} )
    count_N <- rowSums(count)/nrow(count)
    raw_s_factors <- count_q75/count_N
  }
  else if(norm_method == "rle")
  {
    ## RLE(Relative Log Expression normalization)
    ### scale_factors
    ### function for calculating the geometric mean
    geo_mean <- function(x){
      exp(sum(log(x[x>0]))/length(x))
    }
    ### function for calculating non-zero median
    non_zero_median <- function(x){
      median(x[x>0])
    }
    ref_sample <- apply(count, 2, geo_mean)
    norm_rle_1 <- sweep(count, 2, ref_sample, FUN = "/")
    raw_s_factors <- apply(as.matrix(norm_rle_1), 1, non_zero_median)
  }
  else{
    stop("Please choose a valid normalization method")
  }
  return(raw_s_factors)
}

# run BayesDeep
# y: expression counts vector for one gene
# X: covariate matrix
# G: spot-cell spatial relationship matrix
# s: estimated library size factor
# normalize: If TRUE, return the normalized gene expression levels
run.BayesDeep <- function(y, X, G, s, normalize = TRUE){
  n_covariate <- ncol(X)
  X_all <- X
  G <- as.matrix(G)
  # select cells within spots
  y <- y[rowSums(G) != 0]
  s <- s[rowSums(G) != 0]
  X <- X[colSums(G) == 1, ]
  G <- G[rowSums(G) != 0, ]
  G <- G[, colSums(G) == 1]
  
  # Set parameters for BayesDeep model
  sigma_beta <- 1
  a_phi <- 0.1
  b_phi <- 0.1
    
  beta_initial <- rep(0, n_covariate)
  phi_initial <- 5
    
  # Fit the BayesDeep model
  gamma_initial <- sample(c(0, 1), n_covariate, replace = TRUE)
  
  res <- BayesDeep_beta_vs(y, X, s, G, rowSums(G), a_phi, b_phi, sigma_beta, beta_initial, phi_initial, gamma_initial)
  
  iter <- dim(res$beta)[1]
  # get beta
  beta <- apply(res$beta[(iter/2 + 1):iter, ], 2, mean)
  names(beta) <- colnames(X)
  
  # get gamma
  gamma <- apply(res$gamma[(iter/2 + 1):iter, ], 2, mean)
  names(gamma) <- colnames(X)
  
  # set beta = 0 if gamma = 0
  beta[gamma < 0.05] <- 0
  
  # get single-cell-resolution gene expression
  expr <-  exp(X_all %*% matrix(beta, ncol = 1))
  
  if (normalize == TRUE){
  expr <- expr/mean(expr)
  # for extreme values
  expr[expr > quantile(expr, 0.95)] <- quantile(expr, 0.95)
  }
  return(list(expr = expr, beta = beta, gamma = gamma))
}

# plot ST data
# expr: expression of gene for plotting
# loc: x and y coordinate of spots
plot.st <- function(expr, loc, point_size = 1, gene_name = NULL, main = NULL) {
  if (is.null(main)){main <- gene_name}
  data <- data.frame(expr = as.vector(expr), x = loc[, 1], y = loc[, 2]);
  ggplot(data) + geom_point(mapping = aes(x = x, y = y, color = expr), size = point_size) + 
    coord_fixed(ratio = 1) + scale_color_distiller(palette = "Spectral") + 
    theme_classic() + labs(color = "", title = main) + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(hjust = 0.5),
          legend.position="right")
}



  
  
  
  
  
  
  
#############################################

vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  m = nrow(A)
  n = nrow(B)
  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B) )
}


