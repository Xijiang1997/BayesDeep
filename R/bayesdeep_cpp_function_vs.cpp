#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]
using namespace Rcpp;

Rcpp::List BayesDeep_beta_vs(NumericVector y, arma::mat X, NumericVector s,  arma::mat F, NumericVector m,double a_phi, double b_phi, double sigma_beta, NumericVector beta_initial, double phi_initial, NumericVector gamma_initial);
double dnorm_functions_single(double r, double mu, double sigma);

// [[Rcpp::export]]
Rcpp::List BayesDeep_vs(NumericVector y, arma::mat X, arma::mat Z, NumericVector s,  arma::mat F, NumericVector m,double a_phi, double b_phi, double sigma_alpha, double sigma_beta, NumericVector alpha_initial, NumericVector beta_initial, double phi_initial, NumericVector gamma_initial){
  // Read data information
  int M = X.n_rows;
  int Q = X.n_cols;
  int N = F.n_rows;
  int K = Z.n_cols;
  
  double tau_alpha = 1; 
  double tau_beta = 1; 
  double tau_phi = 1;
  double hastings = 0;
  double a_gamma = 0.5;
  double b_gamma = 0.5; 
  
  // Set algorithm settings
  int iter = 2000;
  int burn = iter/2;
  
  int i, q, qq, qqq, qqqq, l, h, count;
  int count_2 = 10;
  
  // Set storage space
  NumericMatrix alpha_store(iter, K);
  NumericMatrix beta_store(iter, Q);
  NumericMatrix gamma_store(iter, Q);
  NumericMatrix phi_store(iter, 1);
  NumericVector accept_alpha(K);
  NumericVector accept_beta(Q);
  NumericVector accept_gamma(Q);
  int accept_phi = 0;
  
  // Set initial values
  arma::mat beta(Q, 1);
  arma::mat alpha(K, 1);
  NumericVector gamma(Q);
  double phi = phi_initial; 
  
  // Initialization
  for(q = 0; q < K; q++)
  {
    alpha(q, 0) = alpha_initial(q);
  }
  
  for(q = 0; q < Q; q++)
  {
    beta(q, 0) = beta_initial(q);
    gamma(q) = gamma_initial(q);
  }
  
  arma::mat log_mu(M, 1);
  arma::mat log_mu_temp(M, 1);
  arma::mat mu(M, 1);
  arma::mat mu_temp(M, 1);
  arma::mat lambda(N, 1);
  arma::mat lambda_temp(N, 1);
  arma::mat alpha_temp(K, 1);
  arma::mat beta_temp(Q, 1);
  int gamma_temp = 0;
  
 
  
  log_mu = Z * alpha + X * beta;
  
  for(i = 0; i < M; i++)
  {
    mu(i, 0) = exp(log_mu(i, 0));
  }
  
  lambda = F * mu;
  
  for(i = 0; i < N; i++)
  {
    lambda(i, 0) = lambda(i, 0)/m(i)*s(i);
  }
 
  
  // MCMC
  for (i = 0; i < iter; i++) {
    
    // Update phi
    double phi_temp = exp(rnorm(1, log(phi), tau_phi)(0));
    
    hastings = 0;
    for(q = 0; q < N; q++) {
      hastings = hastings + lgamma(y(q) + phi_temp) - lgamma(phi_temp) + phi_temp * log(phi_temp) - (phi_temp + y(q)) * log(lambda(q, 0) + phi_temp);
      hastings = hastings - (lgamma(y(q) + phi) - lgamma(phi) + phi * log(phi) - (phi + y(q)) * log(lambda(q, 0) + phi));
    }
    hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
    hastings = hastings - ((a_phi - 1)*log(phi) - b_phi*phi);
    
    if(hastings >= log(double(rand()%10001)/10000)) {
      phi = phi_temp;
      if (i >= burn) {
        accept_phi++;
      }
    }
    
    
    // Update alpha
    
    for (int k = 0; k < K; k++) {
      
      
      
      for(q = 0; q < K; q++)
      {
        alpha_temp(q, 0) = alpha(q, 0);
      }
      
      alpha_temp(k, 0) = rnorm(1, alpha(k, 0), tau_alpha)(0);
      log_mu_temp = log_mu + Z * (alpha_temp - alpha);
      
      for(q = 0; q < M; q++)
      {
        mu_temp(q, 0) = exp(log_mu_temp(q, 0));
      }
      
      lambda_temp = F * mu_temp;
      
      for(q = 0; q < N; q++)
      {
        lambda_temp(q, 0) = lambda_temp(q, 0)/m(q)*s(q);
      }
      
      hastings = 0;
      for(q = 0; q < N; q++) {
        hastings = hastings + (y(q)*log(lambda_temp(q)) - (phi + y(q))*log(lambda_temp(q) + phi)) - (y(q)*log(lambda(q)) - (phi + y(q))*log(lambda(q) + phi));
      }
      hastings = hastings + (-alpha_temp(k, 0)*alpha_temp(k, 0)/2/sigma_alpha/sigma_alpha);
      hastings = hastings - (-alpha(k, 0)*alpha(k, 0)/2/sigma_alpha/sigma_alpha);
      
      if(hastings >= log(double(rand()%10001)/10000)) {
          alpha(k, 0) = alpha_temp(k, 0);
        if (i >= burn){
          accept_alpha(k) = accept_alpha(k) + 1;
        }
          for(q = 0; q < M; q++)
          {
            mu(q, 0) = mu_temp(q, 0);
            log_mu(q, 0) = log_mu_temp(q, 0);
          }
          
          for(q = 0; q < N; q++)
          {
            lambda(q, 0) = lambda_temp(q, 0);
          }
        }
    }
    
   // Update Beta with variable selection
    for (q = 0; q < Q; q++) {
      
      for(qq = 0; qq < Q; qq++)
      {
        beta_temp(qq, 0) = beta(qq, 0);
      }
      
      if (gamma(q) == 0){
        gamma_temp = 1;
        beta_temp(q, 0) = rnorm(1, 0, tau_beta)(0);
      } else {
        gamma_temp = 0;
        beta_temp(q, 0) = 0.0;
      }
      
      log_mu_temp = log_mu + X * (beta_temp - beta);
      
      for(qq = 0; qq < M; qq++)
      {
        mu_temp(qq, 0) = exp(log_mu_temp(qq, 0));
      }
      
      lambda_temp = F * mu_temp;
      for(qq = 0; qq < N; qq++)
      {
        lambda_temp(qq, 0) = lambda_temp(qq, 0)/m(qq)*s(qq);
      }
      
      hastings = 0;
      for(qq = 0; qq < N; qq++) {
        hastings = hastings + (y(qq)*log(lambda_temp(qq, 0)) - (phi + y(qq))*log(lambda_temp(qq, 0) + phi)) - (y(qq)*log(lambda(qq, 0)) - (phi + y(qq))*log(lambda(qq, 0) + phi));
      }
      
      if (gamma(q) == 0){
        hastings = hastings + log(dnorm_functions_single(beta_temp(q, 0), 0.0, sigma_beta));
        
        hastings = hastings + (lgamma(a_gamma + gamma_temp) + lgamma(b_gamma + 2 - gamma_temp) - lgamma(1 + gamma_temp) - lgamma(2 - gamma_temp));
        hastings = hastings - (lgamma(a_gamma + gamma(q)) + lgamma(b_gamma + 2 - gamma(q)) - lgamma(1 + gamma(q)) - lgamma(2 - gamma(q)));
        hastings = hastings - log(dnorm_functions_single(beta_temp(q, 0), 0.0, tau_beta));
      } else{
          hastings = hastings - log(dnorm_functions_single(beta(q, 0), 0.0, sigma_beta));
          
          hastings = hastings + (lgamma(a_gamma + gamma_temp) + lgamma(b_gamma + 2 - gamma_temp) - lgamma(1 + gamma_temp) - lgamma(2 - gamma_temp));
          hastings = hastings - (lgamma(a_gamma + gamma(q)) + lgamma(b_gamma + 2 - gamma(q)) - lgamma(1 + gamma(q)) - lgamma(2 - gamma(q)));
          hastings = hastings + log(dnorm_functions_single(beta(q, 0), 0.0, tau_beta));
      }
        if(hastings > log(unif_rand()))
        {
          accept_gamma(q)++;
          gamma(q) = gamma_temp;
        }
      // further update beta
        if (gamma(q) == 0){
          beta(q, 0) = 0.0;
          // update mu and lambda
          log_mu = Z * alpha + X * beta;
          for(qq = 0; qq < M; qq++)
          {
            mu(qq, 0) = exp(log_mu(qq, 0));
          }
          lambda = F * mu;
          for(qq = 0; qq < N; qq++)
          {
            lambda(qq, 0) = lambda(qq, 0)/m(qq)*s(qq);
          }
        } else{
          beta_temp(q, 0) = rnorm(1,beta(q, 0), tau_beta/2)(0);
          
          log_mu_temp = log_mu + X * (beta_temp - beta);
          
          for(qq = 0; qq < M; qq++)
          {
            mu_temp(qq, 0) = exp(log_mu_temp(qq, 0));
          }
          
          lambda_temp = F * mu_temp;
          for(qq = 0; qq < N; qq++)
          {
            lambda_temp(qq, 0) = lambda_temp(qq, 0)/m(qq)*s(qq);
          }
          
          hastings = 0;
          for(qq = 0; qq < N; qq++) {
            hastings = hastings + (y(qq)*log(lambda_temp(qq, 0)) - (phi + y(qq))*log(lambda_temp(qq, 0) + phi)) - (y(qq)*log(lambda(qq, 0)) - (phi + y(qq))*log(lambda(qq, 0) + phi));
          }
          hastings = hastings + (-beta_temp(q, 0)*beta_temp(q, 0)/2/sigma_beta/sigma_beta);
          hastings = hastings - (-beta(q, 0)*beta(q, 0)/2/sigma_beta/sigma_beta);
          if(hastings > log(unif_rand()))
          {
            if (i >= burn){
            accept_beta(q)++;}
            beta(q, 0) = beta_temp(q, 0);
            for(qq = 0; qq < M; qq++)
            {
              mu(qq, 0) = mu_temp(qq, 0);
              log_mu(qq, 0) = log_mu_temp(qq, 0);
            }
            
            for(qq = 0; qq < N; qq++)
            {
              lambda(qq, 0) = lambda_temp(qq, 0);
            }
          }
          
        }
    }
    
    // store results
    for(q = 0; q < Q; q++)
    {
      beta_store(i, q) = beta(q, 0);
      gamma_store(i, q) = gamma(q);
    }
    for(q = 0; q < K; q++)
    {
      alpha_store(i, q) = alpha(q, 0);
    }
    phi_store(i, 0) = phi; 
    
    
    // Monitor the process
    // Monitor the process
    if (i*100/iter == count_2)
    {
      Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
}
  for(q = 0; q < Q; q++)
  {
    accept_beta(q) = accept_beta(q) / (iter - burn);
    accept_gamma(q) = accept_gamma(q) / (iter - burn);
  }
  for(q = 0; q < K; q++)
  {
    accept_alpha(q) = accept_alpha(q) / (iter - burn);
  }
  accept_phi = accept_phi / (iter - burn); 
  

  return Rcpp::List::create(Rcpp::Named("beta") = beta_store, Rcpp::Named("gamma") = gamma_store, Rcpp::Named("alpha") = alpha_store, Rcpp::Named("phi") = phi_store, Rcpp::Named("accept_phi") = accept_phi, Rcpp::Named("accept_beta") = accept_beta,Rcpp::Named("accept_gamma") = accept_gamma, Rcpp::Named("accept_alpha") = accept_alpha);
}



// [[Rcpp::export]]
Rcpp::List BayesDeep_beta_vs(NumericVector y, arma::mat X, NumericVector s,  arma::mat F, NumericVector m,double a_phi, double b_phi, double sigma_beta, NumericVector beta_initial, double phi_initial, NumericVector gamma_initial){
  // Read data information
  int M = X.n_rows;
  int Q = X.n_cols;
  int N = F.n_rows;
  
  double tau_beta = 1; 
  double tau_phi = 1;
  double hastings = 0;
  double a_gamma = 0.5;
  double b_gamma = 0.5; 
  
  // Set algorithm settings
  int iter = 2000;
  int burn = iter/2;
  
  int i, q, qq, qqq, qqqq, l, h, count;
  int count_2 = 10;
  
  // Set storage space
  NumericMatrix beta_store(iter, Q);
  NumericMatrix gamma_store(iter, Q);
  NumericMatrix phi_store(iter, 1);
  NumericVector accept_beta(Q);
  NumericVector accept_gamma(Q);
  int accept_phi = 0;
  
  // Set initial values
  arma::mat beta(Q, 1);
  NumericVector gamma(Q);
  double phi = phi_initial; 
  
  // Initialization
  for(q = 0; q < Q; q++)
  {
    beta(q, 0) = beta_initial(q);
    gamma(q) = gamma_initial(q);
  }
  
  arma::mat log_mu(M, 1);
  arma::mat log_mu_temp(M, 1);
  arma::mat mu(M, 1);
  arma::mat mu_temp(M, 1);
  arma::mat lambda(N, 1);
  arma::mat lambda_temp(N, 1);
  arma::mat beta_temp(Q, 1);
  int gamma_temp = 0;
  
  
  
  log_mu = X * beta;
  
  for(i = 0; i < M; i++)
  {
    mu(i, 0) = exp(log_mu(i, 0));
  }
  
  lambda = F * mu;
  
  for(i = 0; i < N; i++)
  {
    lambda(i, 0) = lambda(i, 0)/m(i)*s(i);
  }
  
  
  // MCMC
  for (i = 0; i < iter; i++) {
    
    // Update phi
    double phi_temp = exp(rnorm(1, log(phi), tau_phi)(0));
    
    hastings = 0;
    for(q = 0; q < N; q++) {
      hastings = hastings + lgamma(y(q) + phi_temp) - lgamma(phi_temp) + phi_temp * log(phi_temp) - (phi_temp + y(q)) * log(lambda(q, 0) + phi_temp);
      hastings = hastings - (lgamma(y(q) + phi) - lgamma(phi) + phi * log(phi) - (phi + y(q)) * log(lambda(q, 0) + phi));
    }
    hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
    hastings = hastings - ((a_phi - 1)*log(phi) - b_phi*phi);
    
    if(hastings >= log(double(rand()%10001)/10000)) {
      phi = phi_temp;
      if (i >= burn) {
        accept_phi++;
      }
    }
    
    
    
    // Update Beta with variable selection
    for (q = 0; q < Q; q++) {
      
      for(qq = 0; qq < Q; qq++)
      {
        beta_temp(qq, 0) = beta(qq, 0);
      }
      
      if (gamma(q) == 0){
        gamma_temp = 1;
        beta_temp(q, 0) = rnorm(1, 0, tau_beta)(0);
      } else {
        gamma_temp = 0;
        beta_temp(q, 0) = 0.0;
      }
      
      log_mu_temp = log_mu + X * (beta_temp - beta);
      
      for(qq = 0; qq < M; qq++)
      {
        mu_temp(qq, 0) = exp(log_mu_temp(qq, 0));
      }
      
      lambda_temp = F * mu_temp;
      for(qq = 0; qq < N; qq++)
      {
        lambda_temp(qq, 0) = lambda_temp(qq, 0)/m(qq)*s(qq);
      }
      
      hastings = 0;
      for(qq = 0; qq < N; qq++) {
        hastings = hastings + (y(qq)*log(lambda_temp(qq, 0)) - (phi + y(qq))*log(lambda_temp(qq, 0) + phi)) - (y(qq)*log(lambda(qq, 0)) - (phi + y(qq))*log(lambda(qq, 0) + phi));
      }
      
      if (gamma(q) == 0){
        hastings = hastings + log(dnorm_functions_single(beta_temp(q, 0), 0.0, sigma_beta));
        
        hastings = hastings + (lgamma(a_gamma + gamma_temp) + lgamma(b_gamma + 2 - gamma_temp) - lgamma(1 + gamma_temp) - lgamma(2 - gamma_temp));
        hastings = hastings - (lgamma(a_gamma + gamma(q)) + lgamma(b_gamma + 2 - gamma(q)) - lgamma(1 + gamma(q)) - lgamma(2 - gamma(q)));
        hastings = hastings - log(dnorm_functions_single(beta_temp(q, 0), 0.0, tau_beta));
      } else{
        hastings = hastings - log(dnorm_functions_single(beta(q, 0), 0.0, sigma_beta));
        
        hastings = hastings + (lgamma(a_gamma + gamma_temp) + lgamma(b_gamma + 2 - gamma_temp) - lgamma(1 + gamma_temp) - lgamma(2 - gamma_temp));
        hastings = hastings - (lgamma(a_gamma + gamma(q)) + lgamma(b_gamma + 2 - gamma(q)) - lgamma(1 + gamma(q)) - lgamma(2 - gamma(q)));
        hastings = hastings + log(dnorm_functions_single(beta(q, 0), 0.0, tau_beta));
      }
      if(hastings > log(unif_rand()))
      {
        accept_gamma(q)++;
        gamma(q) = gamma_temp;
      }
      // further update beta
      if (gamma(q) == 0){
        beta(q, 0) = 0.0;
        // update mu and lambda
        log_mu = X * beta;
        for(qq = 0; qq < M; qq++)
        {
          mu(qq, 0) = exp(log_mu(qq, 0));
        }
        lambda = F * mu;
        for(qq = 0; qq < N; qq++)
        {
          lambda(qq, 0) = lambda(qq, 0)/m(qq)*s(qq);
        }
      } else{
        beta_temp(q, 0) = rnorm(1,beta(q, 0), tau_beta/2)(0);
        
        log_mu_temp = log_mu + X * (beta_temp - beta);
        
        for(qq = 0; qq < M; qq++)
        {
          mu_temp(qq, 0) = exp(log_mu_temp(qq, 0));
        }
        
        lambda_temp = F * mu_temp;
        for(qq = 0; qq < N; qq++)
        {
          lambda_temp(qq, 0) = lambda_temp(qq, 0)/m(qq)*s(qq);
        }
        
        hastings = 0;
        for(qq = 0; qq < N; qq++) {
          hastings = hastings + (y(qq)*log(lambda_temp(qq, 0)) - (phi + y(qq))*log(lambda_temp(qq, 0) + phi)) - (y(qq)*log(lambda(qq, 0)) - (phi + y(qq))*log(lambda(qq, 0) + phi));
        }
        hastings = hastings + (-beta_temp(q, 0)*beta_temp(q, 0)/2/sigma_beta/sigma_beta);
        hastings = hastings - (-beta(q, 0)*beta(q, 0)/2/sigma_beta/sigma_beta);
        if(hastings > log(unif_rand()))
        {
          if (i >= burn){
            accept_beta(q)++;}
          beta(q, 0) = beta_temp(q, 0);
          for(qq = 0; qq < M; qq++)
          {
            mu(qq, 0) = mu_temp(qq, 0);
            log_mu(qq, 0) = log_mu_temp(qq, 0);
          }
          
          for(qq = 0; qq < N; qq++)
          {
            lambda(qq, 0) = lambda_temp(qq, 0);
          }
        }
        
      }
    }
    
    // store results
    for(q = 0; q < Q; q++)
    {
      beta_store(i, q) = beta(q, 0);
      gamma_store(i, q) = gamma(q);
    }
   
    phi_store(i, 0) = phi; 
    
    // Monitor the process
    if (i*100/iter == count_2)
    {
      Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
  }
  for(q = 0; q < Q; q++)
  {
    accept_beta(q) = accept_beta(q) / (iter - burn);
    accept_gamma(q) = accept_gamma(q) / (iter - burn);
  }
  
  accept_phi = accept_phi / (iter - burn); 
  
  
  return Rcpp::List::create(Rcpp::Named("beta") = beta_store, Rcpp::Named("gamma") = gamma_store,  Rcpp::Named("phi") = phi_store, Rcpp::Named("accept_phi") = accept_phi, Rcpp::Named("accept_beta") = accept_beta,Rcpp::Named("accept_gamma") = accept_gamma);
}



// [[Rcpp::export]]
double dnorm_functions_single(double r, double mu, double sigma){
  arma::mat rr(1,1);
  rr(0, 0) = r; 
  NumericVector mu1(1);
  arma::mat sigma1(1, 1);
  mu1(0) = mu;
  sigma1(0, 0) = sigma*sigma;
  arma::vec qq = dmvnorm(rr, mu1, sigma1);
    
return(qq(0, 0));
}

