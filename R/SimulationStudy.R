## Required packages
library(pROC)
library(tidyr)
library(dplyr)
library(ggplot2)
library(magrittr)
library(ggpubr)
library(tidypaleo)


## Solve cubic equation: excast
odds_adjust <- function(p0, p1, v)
{
  
  if (v > p0 * (1 - p0)) return(NA)
  else {
    require(RConics)
    A <- p0^3 - p1*p0^3
    B <- 3*p1*p0^3 - 2*p0^3 - 3*p1*p0^2 + 2*p0^2 - v
    C <- p0^3 - 3*p1*p0^3 + 6*p1*p0^2 - 2*p0^2 + p0-3*p1*p0 + v
    D <- p1*p0^3 - 3*p1*p0^2 + 3*p1*p0 - p1
    
    res <- cubic(c(A,B,C,D))
    res <- Re(res[which(Im(res)==0)]) #Remove non-reals
    res <- res[which(res>0)] #Remove negatives
    res <- res[which(sign(log(res)) == sign(log(p1/p0)))] #Removes ORs in the wrong direction
    return(res)
  }
}

## Solve cubic equation: numerical
odds_adjust_numerical <- function (p0, p1, v)
{
  f <- function(x,p0,p1,v)
  {
    x*p0/(1-(1-x)*p0)+(1-x)*x/(1-(1-x)*p0)^3*v-p1
  }
  
  if(p1<p0) inerval <- c(0,1) else interval  <- c(1,10)
  
  uniroot(f,interval,p0,p1,v)$root
}

## Simulation study runs function
run_sims <- function(p0 = 0.5,    # outcome risk in development population
                     delta = 10,  # difference between outcome risk in development and target populations
                     n = 10000)   # number of simulated outcome risk from a beta distribution
  {

  ## Required functions
  EstOR_compare <- function(p0 = p0, delta = delta,
                            n = n, v_length = 50) {
    
    ## Var range
    v_vals <- seq(0, (p0) * (1 - p0), length.out = v_length)
    
    res <- data.frame(Var = v_vals,
                      p1 = NA,
                      AUC = NA,
                      Exact_OR = NA,
                      Exact_p = NA,
                      Naive_OR = NA,
                      Naive_p = NA,
                      TaylorApprox_OR = NA,
                      TaylorApprox_p = NA)
    
    
    for (i in 1 : length(v_vals)) {
      
      v <- v_vals[i]
      
      if (v == 0 | v >= (p0 * (1 - p0) - 0.01)) res[i , c(3 : 9)] <- NA
      else {
        
        ## exact estimate (sip0lation from dist. of pi)
        pi_sim <- rbeta(n = n,
                        shape1 = (p0 * (1 - p0) / v - 1) * p0,
                        shape2 = (p0 * (1 - p0) / v - 1) * (1 - p0))
        
        pi_sim[pi_sim == 1] <- 0.99
        pi_sim[pi_sim == 0] <- 0.01
        n_sim_adj <- length(pi_sim)
        
        l_OR_sim <- log(pi_sim / (1 - pi_sim))
        pi_mean <- mean(pi_sim)
        pi_var <- var(pi_sim)
        p1 <- mean(pi_sim) + mean(pi_sim) * delta / 100
        res$p1[i] <- p1
        
        y_temp <- c(rep(1, round(n_sim_adj * p1)), rep(0, n_sim_adj - round(n_sim_adj * p1)))
        glm_temp <- glm(y_temp ~ offset(l_OR_sim), family = binomial(link = "logit"))
        res$Exact_OR[i] <- exp(coef(glm_temp))
        res$Exact_p[i] <- mean(pi_sim * exp(coef(glm_temp)) / (1 + pi_sim * exp(coef(glm_temp)) - pi_sim))
        
        pi_sim_adj <- ((pi_sim / (1 - pi_sim)) * res$Exact_OR[i]) / (1 + ((pi_sim / (1 - pi_sim)) * res$Exact_OR[i]))
        y_sim <- rbinom(n = n, size = 1, prob = pi_sim_adj)
        res$AUC[i] <- ifelse(mean(y_sim) > 0, try(auc(y_sim, pi_sim_adj), silent = T), NA)
        
        ## prevalance estimate (Janssen's suggestion)
        res$Naive_OR[i] <- (p1 / (1 - p1)) / (pi_mean / (1 - pi_mean))
        res$Naive_p[i] <- mean(pi_sim * res$Naive_OR[i] / (1 + pi_sim * res$Naive_OR[i] - pi_sim))
        
        ## our estimate (Taylor approximation)
        res$TaylorApprox_OR[i] <- odds_adjust(p0 = pi_mean, p1 = p1, v = pi_var)
        res$TaylorApprox_p[i] <- mean(pi_sim * res$TaylorApprox_OR[i] / (1 + pi_sim * res$TaylorApprox_OR[i] - pi_sim))
      }
      
    }
    
    return(res)
  }
  
  ggp_Creator <- function(data, p0, delta, AUC_max = 0.95, y_lim = NULL) {
    
    data_temp <- data[data$p0 == p0 & data$Delta == delta & !is.na(data$AUC) &
                        data$AUC <= AUC_max, ]
    data_temp %<>%
      filter((Method == "Exact adjustment" & Var_sim > 0.0019 & p0 == 0.1) |
               (Method == "Exact adjustment" & Var_sim > 0.0039 & p0 == 0.25) |
               (Method == "Exact adjustment" & Var_sim > 0.0052 & p0 == 0.50) |
               (Method != "Exact adjustment"))
    
    data_temp$AUC_Dlag <- data_temp$AUC - lag(data_temp$AUC)
    data_temp$AUC_Dlag[is.na(data_temp$AUC_Dlag)] <- 0
    data_temp <- data_temp[data_temp$AUC_Dlag >= 0 , ]
    
    model_temp <- age_depth_model(
      depth = data_temp$AUC,
      age = data_temp$Var_sim)
    
    ggp_temp <- ggplot(data = data_temp,
                       aes(x = Var_sim, y = p_error, group = Method, color = Method)) +
      geom_line() +
      labs(x = "Var", y = "Relative bias (%)") +
      theme_bw() +
      theme(legend.position = "bottom", text = element_text(size = 15, face = "bold"),
            axis.text.x.top = element_text(size = 8, face = "bold", angle = 0),
            axis.text.x.bottom = element_text(size = 11, face = "bold")) +
      facet_wrap(~ p0_lab + p1_lab, labeller = label_wrap_gen(multi_line = FALSE)) +
      scale_x_age_depth(model_temp, depth_name = "AUC")
    if (! is.null(y_lim)) ggp_temp <- ggp_temp + lims(y = y_lim)
    
    return(ggp_temp)
  }
  

  
  #### Simulation study
  simRes <- data.frame(p0 = rep(p0, 50),
                       Delta = rep(delta, 50),
                       AUC = NA, Var_sim = NA, p_sim = NA,
                       Exact_OR = NA, Exact_p = NA,
                       Naive_OR = NA, Naive_p =  NA,
                       TaylorApprox_OR = NA, TaylorApprox_p = NA)
  
  simRes[ , c("Var_sim", "p_sim", "AUC", "Exact_OR", "Exact_p",
              "Naive_OR", "Naive_p", "TaylorApprox_OR", "TaylorApprox_p")] <-
    EstOR_compare(p0 = p0, delta = delta, n = n)

  simRes <- simRes[! is.na(simRes$p_sim) , ]
  simRes$Exact_p_error_rel <- (simRes$p_sim - simRes$Exact_p) / simRes$p_sim * 100
  simRes$Naive_p_error_rel <- (simRes$p_sim - simRes$Naive_p) / simRes$p_sim * 100
  simRes$TaylorApprox_p_error_rel <- (simRes$p_sim - simRes$TaylorApprox_p) / simRes$p_sim * 100
  
  simRes %>%
    dplyr::select(p0, Delta, AUC, Var_sim, p_sim,
                  Exact_p_error_rel, Naive_p_error_rel, TaylorApprox_p_error_rel) %>%
    tidyr::gather(Estimator, p_error, -c(p0, Delta, AUC, Var_sim, p_sim)) -> simRes_long
  
  simRes_long$Method <- ifelse(simRes_long$Estimator == "Exact_p_error_rel", "Exact adjustment",
                               ifelse(simRes_long$Estimator == "Naive_p_error_rel", "Simple adjustment",
                                      "Taylor adjustment"))
  simRes_long$Method <- factor(simRes_long$Method,
                               levels = c("Exact adjustment",
                                          "Taylor adjustment",
                                          "Simple adjustment"))
  simRes_long$p0_lab <- paste0("p0=", simRes_long$p0)
  simRes_long$p1_lab <- paste0("Delta=", simRes_long$Delta, "%")
  
  
  return(ggp_Creator(data = simRes_long, p0 = p0, delta = delta))
    
}


##### Examples
## simulation runs
run_sims(p0 = 0.1, delta = -10)  ## top-left panel in the manuscript
run_sims(p0 = 0.5, delta = 50)   ## bottom-right panel in the manuscript


## Adjust odds ratio
# exact
odds_adjust(p0 = 0.1, p1 = 0.12, v = 0.01)
# numerical
odds_adjust_numerical(p0 = 0.1, p1 = 0.12, v = 0.01)

