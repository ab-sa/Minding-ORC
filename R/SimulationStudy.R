## Required packages
library(pROC)
library(tidyr)
library(dplyr)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(tidypaleo)

## Required functions
solve_or <- function(p, q, v)
{
  
  if (v > p * (1 - p)) return(NA)
  else {
    require(RConics)
    A <- p^3 - q*p^3
    B <- 3*q*p^3 - 2*p^3 - 3*q*p^2 + 2*p^2 - v
    C <- p^3 - 3*q*p^3 + 6*q*p^2 - 2*p^2 + p-3*q*p + v
    D <- q*p^3 - 3*q*p^2 + 3*q*p - q
    
    res <- cubic(c(A,B,C,D))
    res <- Re(res[which(Im(res)==0)]) #Remove non-reals
    res <- res[which(res>0)] #Remove negatives
    res <- res[which(sign(log(res)) == sign(log(q/p)))] #Removes ORs in the wrong direction
    return(res)
  }
}

EstOR_compare <- function(q_changeRate = 10, mu = 0.5, pi_val = NULL, mu_val = NULL,
                          n_sim = 10000, v_length = 50) {
  
  ## Var range
  v_vals <- seq(0, (mu) * (1 - mu), length.out = v_length)
  
  res <- data.frame(Var = v_vals,
                    q = NA,
                    AUC = NA,
                    Exact_OR = NA,
                    Exact_q = NA,
                    Naive_OR = NA,
                    Naive_q = NA,
                    TaylorApprox_OR = NA,
                    TaylorApprox_q = NA)
  
  
  for (i in 1 : length(v_vals)) {
    
    v <- v_vals[i]
    
    if (v == 0 | v >= (mu * (1 - mu) - 0.01)) res[i , c(3 : 9)] <- NA
    else {
      
      ## exact estimate (simulation from dist. of pi)
      pi_sim <- rbeta(n = n_sim,
                      shape1 = (mu * (1 - mu) / v - 1) * mu,
                      shape2 = (mu * (1 - mu) / v - 1) * (1 - mu))
      
      pi_sim[pi_sim == 1] <- 0.99
      pi_sim[pi_sim == 0] <- 0.01
      n_sim_adj <- length(pi_sim)
      
      l_OR_sim <- log(pi_sim / (1 - pi_sim))
      pi_mean <- mean(pi_sim)
      pi_var <- var(pi_sim)
      q <- mean(pi_sim) + mean(pi_sim) * q_changeRate / 100
      res$q[i] <- q
      
      y_temp <- c(rep(1, round(n_sim_adj * q)), rep(0, n_sim_adj - round(n_sim_adj * q)))
      glm_temp <- glm(y_temp ~ offset(l_OR_sim), family = binomial(link = "logit"))
      res$Exact_OR[i] <- exp(coef(glm_temp))
      res$Exact_q[i] <- mean(pi_sim * exp(coef(glm_temp)) / (1 + pi_sim * exp(coef(glm_temp)) - pi_sim))
      
      pi_sim_adj <- ((pi_sim / (1 - pi_sim)) * res$Exact_OR[i]) / (1 + ((pi_sim / (1 - pi_sim)) * res$Exact_OR[i]))
      y_sim <- rbinom(n = n_sim, size = 1, prob = pi_sim_adj)
      res$AUC[i] <- ifelse(mean(y_sim) > 0, try(auc(y_sim, pi_sim_adj), silent = T), NA)
      
      ## prevalance estimate (Janssen's suggestion)
      res$Naive_OR[i] <- (q / (1 - q)) / (pi_mean / (1 - pi_mean))
      res$Naive_q[i] <- mean(pi_sim * res$Naive_OR[i] / (1 + pi_sim * res$Naive_OR[i] - pi_sim))
      
      ## our estimate (Taylor approximation)
      res$TaylorApprox_OR[i] <- solve_or(p = pi_mean, q = q, v = pi_var)
      res$TaylorApprox_q[i] <- mean(pi_sim * res$TaylorApprox_OR[i] / (1 + pi_sim * res$TaylorApprox_OR[i] - pi_sim))
    }
    
  }
  
  return(res)
}

ggp_Creator <- function(data = simRes_long,
                        mu = 0.1, q = 10, AUC_max = 0.95,
                        y_range = NULL) {
  
  data_temp <- data[data$mu == mu & data$Delta == q & !is.na(data$AUC) &
                      data$AUC <= AUC_max, ]
  data_temp %<>%
    filter((Method == "Exact adjustment" & Var_sim > 0.0019 & mu == 0.1) |
             (Method == "Exact adjustment" & Var_sim > 0.0039 & mu == 0.25) |
             (Method == "Exact adjustment" & Var_sim > 0.0052 & mu == 0.50) |
             (Method != "Exact adjustment"))
  
  data_temp$AUC_Dlag <- data_temp$AUC - lag(data_temp$AUC)
  data_temp$AUC_Dlag[is.na(data_temp$AUC_Dlag)] <- 0
  data_temp <- data_temp[data_temp$AUC_Dlag >= 0 , ]
  
  model_temp <- age_depth_model(
    depth = data_temp$AUC,
    age = data_temp$Var_sim)
  
  ggp_temp <- ggplot(data = data_temp,
                     aes(x = Var_sim, y = q_error, group = Method, color = Method)) +
    geom_line() +
    labs(x = "Var", y = "Relative bias (%)") +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 15, face = "bold"),
          axis.text.x.top = element_text(size = 8, face = "bold", angle = 0),
          axis.text.x.bottom = element_text(size = 11, face = "bold")) +
    facet_wrap(~ p0_lab + p1_lab, labeller = label_wrap_gen(multi_line = FALSE)) +
    scale_x_age_depth(model_temp, depth_name = "AUC")
  if (! is.null(y_range)) ggp_temp <- ggp_temp + lims(y = y_range)
  
  return(ggp_temp)
}


#### Simulation study
mu_vals <- c(0.1, 0.25, 0.50)
q_vals <- c(-10, -25, -50, 10, 25, 50)

simRes <- data.frame(Run = rep(c(1 : (length(q_vals) * length(mu_vals))), each = 50),
                     mu = rep(mu_vals, each = 50 * length(q_vals)),
                     Delta = rep(rep(q_vals, each = 50), length(mu_vals)), AUC = NA,
                     Var_sim = NA, q_sim = NA,
                     Exact_OR = NA, Exact_q = NA,
                     Naive_OR = NA, Naive_q =  NA,
                     TaylorApprox_OR = NA, TaylorApprox_q = NA)

set.seed(1)

for (i in 1 : length(mu_vals)) {
  
  for (j in 1 : length(q_vals)) {
    
    simRes[((i - 1) * length(q_vals) * 50 + (j - 1) * 50 + 1) : ((i - 1) * length(q_vals) * 50 + (j) * 50) ,
           c("Var_sim", "q_sim", "AUC", "Exact_OR", "Exact_q",
             "Naive_OR", "Naive_q", "TaylorApprox_OR", "TaylorApprox_q")] <-
      EstOR_compare(q_changeRate = q_vals[j], mu = mu_vals[i])
    
  }
  
}

simRes <- simRes[! is.na(simRes$q_sim) , ]
simRes$Exact_q_error_rel <- (simRes$q_sim - simRes$Exact_q) / simRes$q_sim * 100
simRes$Naive_q_error_rel <- (simRes$q_sim - simRes$Naive_q) / simRes$q_sim * 100
simRes$TaylorApprox_q_error_rel <- (simRes$q_sim - simRes$TaylorApprox_q) / simRes$q_sim * 100

simRes %>%
  dplyr::select(Run, mu, Delta, AUC, Var_sim, q_sim,
                Exact_q_error_rel, Naive_q_error_rel, TaylorApprox_q_error_rel) %>%
  tidyr::gather(Estimator, q_error, -c(Run, mu, Delta, AUC, Var_sim, q_sim)) -> simRes_long

simRes_long$Method <- ifelse(simRes_long$Estimator == "Exact_q_error_rel", "Exact adjustment",
                             ifelse(simRes_long$Estimator == "Naive_q_error_rel", "Simple adjustment",
                                    "Taylor adjustment"))
simRes_long$Method <- factor(simRes_long$Method,
                             levels = c("Exact adjustment",
                                        "Taylor adjustment",
                                        "Simple adjustment"))
simRes_long$p0_lab <- paste0("p0=", simRes_long$mu)
simRes_long$p1_lab <- paste0("Delta=", simRes_long$Delta, "%")


### Visualizations
ggp_Creator(data = simRes_long, mu = mu_vals[1], q = q_vals[1])
ggp_Creator(data = simRes_long, mu = mu_vals[2], q = q_vals[4])

