fit_models <- function(model_func, outcome_vars, rhs_formula, data, ...) {
  models <- lapply(outcome_vars, function(outcome_var) {
    formula <- reformulate(rhs_formula, response = outcome_var)
    model <- do.call(model_func, list(formula = formula, data = data, ...))
    return(model)
  })
  names(models) <- outcome_vars
  return(models)
}

outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}


p_sobel <- function(a, s_a, b, s_b){
  z = (a*b)/(sqrt(b^2 * s_a^2 + a^2 * s_b^2))
  return (1- pnorm(z))*2
}

sobel_test <- function(Hypothesis, model_med, vcov_med, model_out, vcov_out, predictor, mediator, Predictor, Mediator, stars=FALSE, method, decimal_places = 3, model_orig, vcov_orig){
  #replace vcov by cluster robust version
  model_med@vcov_beta <- as.matrix(vcov_med) 
  model_out@vcov_beta <- as.matrix(vcov_out)
  model_orig@vcov_beta <- as.matrix(vcov_orig)
  
  parameters_med <- data.frame(summary(model_med)$coefficients, ddf = method) %>% rownames_to_column()
  a <- parameters_med %>% filter(rowname == predictor) %>% .[[1,"Estimate"]]
  s_a <- parameters_med %>% filter(rowname == predictor) %>% .[[1,"Std..Error"]]
  
  parameters_out <- data.frame(summary(model_out)$coefficients, ddf = method) %>% rownames_to_column()
  b <- parameters_out %>% filter(rowname == mediator) %>% .[[1,"Estimate"]]
  s_b <- parameters_out %>% filter(rowname == mediator) %>% .[[1,"Std..Error"]]
  c_dash <- parameters_out %>% filter(rowname == predictor) %>% .[[1,"Estimate"]]
  
  parameters_orig <- data.frame(summary(model_orig)$coefficients, ddf = method) %>% rownames_to_column()
  c <- parameters_orig %>% filter(rowname == predictor) %>% .[[1,"Estimate"]]
  
  if(stars){
    p_a <- parameters_med %>% filter(rowname == predictor) %>% .[[1,"Pr...t.."]]
    p_b <- parameters_out %>% filter(rowname == mediator) %>% .[[1,"Pr...t.."]]
    p_c_dash <- parameters_out %>% filter(rowname == predictor) %>% .[[1,"Pr...t.."]]
    p_c <- parameters_orig %>% filter(rowname == predictor) %>% .[[1,"Pr...t.."]]
    a_star <- paste0(round(a,decimal_places), ifelse(p_a<0.1,
                              ifelse(p_a<0.05,
                                     ifelse(p_a<0.01,
                                            ifelse(p_a<0.001,
                                                       "***",
                                                       "**"
                                                ),
                                            "*"),
                                     "+"),
                              ""))
    b_star <- paste0(round(b,decimal_places), ifelse(p_b<0.1,
                               ifelse(p_b<0.05,
                                      ifelse(p_b<0.01,
                                             ifelse(p_b<0.001,
                                                    "***",
                                                    "**"
                                             ),
                                             "*"),
                                      "+"),
                               ""))
    c_dash_star <- paste0(round(c_dash,decimal_places), ifelse(p_c_dash<0.1,
                               ifelse(p_c_dash<0.05,
                                      ifelse(p_c_dash<0.01,
                                             ifelse(p_c_dash<0.001,
                                                    "***",
                                                    "**"
                                             ),
                                             "*"),
                                      "+"),
                               ""))
    c_star <- paste0(round(c,decimal_places), ifelse(p_c<0.1,
                                                     ifelse(p_c<0.05,
                                                            ifelse(p_c<0.01,
                                                                   ifelse(p_c<0.001,
                                                                          "***",
                                                                          "**"
                                                                   ),
                                                                   "*"),
                                                            "+"),
                                                     ""))
    return(data.frame("Hypothesis" = Hypothesis, "Predictor" = Predictor, "Mediator" = Mediator, "c" = c_star, "a" = a_star, "b" = b_star, "c'" = c_dash_star, "p_Sobel" = p_sobel(a, s_a, b, s_b)))
  }
  
  return(data.frame("Hypothesis" = Hypothesis, "Predictor" = Predictor, "Mediator" = Mediator, "c" = c, "a" = a, "b" = b, "c'" = c_dash, "p_Sobel" = p_sobel(a, s_a, b, s_b)))
}


mediation_table <- function(list_mediation_models, stars = FALSE, method, decimal_places = 3, model_orig, vcov_orig){
  
  
  sobel_results_list <- lapply(list_mediation_models, 
                               FUN = function(x){ 
                                 sobel_test(x$Hypothesis, x$model_med, x$vcov_med, x$model_out, x$vcov_out, x$predictor, x$mediator, x$Predictor, x$Mediator, stars, method, decimal_places, model_orig, vcov_orig)
                                 })
  
  
  
  return(bind_rows(sobel_results_list))
}


# 
significant_ranges <- function (obj, varnames, vcov = NULL, 
                                other_interactions = NULL,
                                alpha = 0.05,
                                find_insignificant = FALSE){
  
  
  if(is.null(vcov)) vcov = vcov(obj)
  
  
  if(class(obj) == "lm"){
    
    if(!is.null(other_interactions)) print("Warning: so far only one interaction term per effect in lm")
    
    if (!("model" %in% names(obj))) {
      obj <- update(obj, model = T)
    }
    
    v1 <- varnames[1]
    v2 <- varnames[2]

    b1.pos <- grep(pattern = v1, x = names(coef(obj)), fixed = TRUE)[1]
    b3.pos <- c(grep(pattern = paste(v1,":",v2, sep = ""), x = names(coef(obj)), fixed = TRUE), 
                grep(pattern = paste(v2,":",v1, sep = ""), x = names(coef(obj)), fixed = TRUE))[1]
    b1 <- as.numeric(obj$coef[b1.pos])
    b3 <- as.numeric(obj$coef[b3.pos])
    var1 <- vcov[b1.pos, b1.pos]
    var3 <- vcov[b3.pos, b3.pos]
    cov13 <- vcov[b1.pos, b3.pos]
    

    
    low <- function(x){
      (b1 + b3 * x)  - qt(1-alpha/2, df = obj$df.residual) * sqrt(var1 + x^2 * var3 + 2 * x * cov13)
    }
    
    up <- function(x){
      (b1 + b3 * x)  + qt(1-alpha/2, df = obj$df.residual) * sqrt(var1 + x^2 * var3 + 2 * x * cov13)
    }
    
    #slack <- 0.1*(max(model.matrix(obj)[, v2]) - min(model.matrix(obj)[, v2]))
    x_range <- c(min(model.matrix(obj)[, v2]), max(model.matrix(obj)[, v2]))
    
    if(find_insignificant){
      return(find_zero_intervals(low,up,x_range))
    } else{
    return(find_non_zero_intervals(low,up,x_range))}
    
  }
  
  if(class(obj) == "lmerMod"){
    
    print("Warning: so far only implemented for lm.")
    
    

    
  }

}

find_non_zero_intervals <- function(lower_bound_func, upper_bound_func, x_range, tol = .Machine$double.eps^0.25) {
  library(purrr)
  
  # Define a function to find roots within a given interval
  find_roots <- function(func, x_range) {
    roots <- c()
    for (i in seq_along(x_range)[-length(x_range)]) {
      x1 <- x_range[i]
      x2 <- x_range[i + 1]
      if (func(x1) * func(x2) < 0) {
        root <- uniroot(func, c(x1, x2), tol = tol)$root
        roots <- c(roots, root)
      }
    }
    return(roots)
  }
  
  # Find zero points for lower and upper bound functions
  x_values <- seq(min(x_range), max(x_range), length.out = 1000)
  lower_zeros <- find_roots(lower_bound_func, x_values)
  upper_zeros <- find_roots(upper_bound_func, x_values)
  
  # Combine and sort all zero points
  all_zeros <- sort(unique(c(lower_zeros, upper_zeros)))
  

  # Initialize variables to store the results
  non_zero_intervals <- list()
  
  # If no zero points are found, check the entire range
  if (length(all_zeros) == 0) {
    mid_x <- mean(x_range)
    if (lower_bound_func(mid_x) > 0 || upper_bound_func(mid_x) < 0) {
      non_zero_intervals <- list(x_range)
    }
  } else {
    # Evaluate intervals between zero points
    interval_borders <- c(x_range[1], all_zeros, x_range[2])
    for (i in seq_along(interval_borders)[-length(interval_borders)]) {
      x1 <- interval_borders[i]
      x2 <- interval_borders[i + 1]
      mid_x <- (x1 + x2) / 2
      if (lower_bound_func(mid_x) > 0 || upper_bound_func(mid_x) < 0) {
        non_zero_intervals <- append(non_zero_intervals, list(sort(c(x1, x2))))
      }
    }
  }
  
  return(non_zero_intervals)
}


find_zero_intervals <- function(lower_bound_func, upper_bound_func, x_range, tol = .Machine$double.eps^0.25) {
  library(purrr)
  
  # Define a function to find roots within a given interval
  find_roots <- function(func, x_range) {
    roots <- c()
    for (i in seq_along(x_range)[-length(x_range)]) {
      x1 <- x_range[i]
      x2 <- x_range[i + 1]
      if (func(x1) * func(x2) < 0) {
        root <- uniroot(func, c(x1, x2), tol = tol)$root
        roots <- c(roots, root)
      }
    }
    return(roots)
  }
  
  # Find zero points for lower and upper bound functions
  x_values <- seq(min(x_range), max(x_range), length.out = 1000)
  lower_zeros <- find_roots(lower_bound_func, x_values)
  upper_zeros <- find_roots(upper_bound_func, x_values)
  
  # Combine and sort all zero points
  all_zeros <- sort(unique(c(lower_zeros, upper_zeros)))
  
  
  # Initialize variables to store the results
  zero_intervals <- list()
  
  # If no zero points are found, check the entire range
  if (length(all_zeros) == 0) {
    mid_x <- mean(x_range)
    if (lower_bound_func(mid_x) < 0 & upper_bound_func(mid_x) > 0) {
      zero_intervals <- list(x_range)
    }
  } else {
    # Evaluate intervals between zero points
    interval_borders <- c(x_range[1], all_zeros, x_range[2])
    for (i in seq_along(interval_borders)[-length(interval_borders)]) {
      x1 <- interval_borders[i]
      x2 <- interval_borders[i + 1]
      mid_x <- (x1 + x2) / 2
      if (lower_bound_func(mid_x) < 0 & upper_bound_func(mid_x) > 0) {
        zero_intervals <- append(zero_intervals, list(sort(c(x1, x2))))
      }
    }
  }
  
  #print(zero_intervals)
  
  return(zero_intervals)
}



# R function to plot conditional effects for linear (and multilevel linear) regression
# Johannes Karreth, customized by Lea Tamberg

# This function was inspired by David A. Armstrong's DAintfun2
#  function (in the DAmisc package).
# It uses ggplot2, adds custom variable names for axis labels (varlabs)
#  and puts the conditional effect of X1/X2 in the subplot title to save space.
# The function also automatically adjusts the plot type for binary vs. continuous
#  moderators.
# obj: lm or lmer object
# varnames: character vector of the constitutive terms, e.g. c("x1", "x2")
# varlabs: character vector of length 2 with the desired variable names for
#  x1 and x2 to show up as labels in the plot.



# theme_jk <- function(base_size = 12, base_family = "sans") #"Latin Modern Roman 10"
# {
#   theme_bw(base_size = base_size, base_family = base_family) %+replace%
#     theme(
#       strip.background = element_rect(fill = NA, color = "black", size = 0.75),
#       strip.text.x = element_text(size = rel(1), margin = margin(t = base_size/2 * 0.6, b = base_size/2 * 0.6), face = "bold"),
#       panel.border = element_rect(fill = NA, color = "black", size = 0.75),
#       plot.title = element_text(size = rel(1.2), margin = margin(b = base_size/2 * 1.2), face = "bold"),
#       panel.grid.major = element_line(colour = "grey85", size = 0.4, linetype = 3), 
#       panel.grid.minor = element_line(colour = "grey90", size = 0.2, linetype = 3),
#       axis.text = element_text(color = "black")
#     )
# }


ggintplot <- function (obj, varnames, vcov = NULL, varlabs, ref = NULL, 
                       alt_ref = NULL, ref_lab = NULL,
          other_interactions = NULL,
          data = NULL,
          category_name = NULL,
          alpha = 0.05,
          title = FALSE,
          subtitle = NULL,
          rug = FALSE, 
          twoways = FALSE, 
          rugwidth = 0.1, 
          jitter_factor = 0,
          max_x = NULL,
          unlog_x = FALSE,
          exp_y = FALSE,
          log_scale_x_axis = FALSE,
          font_family = "sans",
          font_size = 12,
          color_rug = FALSE,
          color_lines = color_rug,
          color_band = FALSE,
          color_rug_v1 = "#004488",
          color_rug_ref = "#DDAA33",
          legend = TRUE,
          rug_alpha = 1){
  
  require(ggplot2); require(gridExtra)
  
  if(is.null(vcov)) vcov = vcov(obj)
  
  if(!is.null(other_interactions) & twoways) print("Warning: other_interactions will only be considered in the first marginal effects plot.")
  if(!is.null(alt_ref) & twoways) print("Warning: alt_ref will only be considered in the first marginal effects plot.")
  
  
  v1 <- varnames[1]
  v2 <- varnames[2]
  vlab1 <- varlabs[1]
  vlab2 <- varlabs[2]
  if(is.null(max_x)) max_x <- max(model.matrix(obj)[, v2])
  
  if(length(unique(model.matrix(obj)[, paste(v1)])) == 2){
    s1 <- c(min(model.matrix(obj)[, paste(v1)]), max(model.matrix(obj)[, paste(v1)]))
  }
  
  if(length(unique(model.matrix(obj)[, paste(v1)])) > 2){
    s1 <- seq(from = min(model.matrix(obj)[, paste(v1)]), to = max(model.matrix(obj)[, paste(v1)]), length.out = 25)
  }
  
  if(length(unique(model.matrix(obj)[, paste(v2)])) == 2){
    s2 <- c(min(model.matrix(obj)[, paste(v2)]), max(model.matrix(obj)[, paste(v2)]))
  }
  
  if(length(unique(model.matrix(obj)[, paste(v2)])) > 2){
    s2 <- seq(from = min(model.matrix(obj)[, v2]), to = max_x, length.out = 25)
  }
  
  
  
  if(class(obj) == "lm"){
    
    if(!is.null(other_interactions)) print("Warning: so far only one interaction term per effect in lm")
    
    if (!("model" %in% names(obj))) {
      obj <- update(obj, model = T)
    }
    
    
    
    # relevant parameter values and standard deviations for FIRST marginal effects plot
    b1.pos <- grep(pattern = v1, x = names(coef(obj)), fixed = TRUE)[1]
    b3.pos <- c(grep(pattern = paste(v1,":",v2, sep = ""), x = names(coef(obj)), fixed = TRUE), 
                grep(pattern = paste(v2,":",v1, sep = ""), x = names(coef(obj)), fixed = TRUE))[1]
    b1 <- as.numeric(obj$coef[b1.pos])
    b3 <- as.numeric(obj$coef[b3.pos])
    var1 <- vcov[b1.pos, b1.pos]
    var3 <- vcov[b3.pos, b3.pos]
    cov13 <- vcov[b1.pos, b3.pos]
    
    # calculate marginal effects and confidence band for FIRST marginal effects plot
    if(is.null(alt_ref)){ #standard case
      eff1 <- b1 + b3 * s2
      var.eff1 <- var1 + s2^2 * var3 + 2 * s2 * cov13
      se.eff1 <- sqrt(var.eff1)
      low1 <- eff1 - qt(1-alpha/2, df = obj$df.residual) * se.eff1
      up1 <- eff1 + qt(1-alpha/2, df = obj$df.residual) * se.eff1
    } else { # case with alternative reference category (works only for categorical moderator)
      # more parameters and standard deviations are needed
      b1_ref.pos <- grep(pattern = alt_ref, x = names(coef(obj)), fixed = TRUE)[1]
      b3_ref.pos <- c(grep(pattern = paste(alt_ref,":",v2, sep = ""), x = names(coef(obj)), fixed = TRUE), 
                  grep(pattern = paste(v2,":", alt_ref, sep = ""), x = names(coef(obj)), fixed = TRUE))[1]
      
      b1_ref <- as.numeric(obj$coef[b1_ref.pos])
      b3_ref <- as.numeric(obj$coef[b3_ref.pos])

      var1_ref <- vcov[b1_ref.pos, b1_ref.pos]
      var3_ref <- vcov[b3_ref.pos, b3_ref.pos]
      cov1_ref_3_ref <- vcov[b1_ref.pos, b3_ref.pos]
      cov13_ref <- vcov[b1.pos, b3_ref.pos]
      cov1_ref_3 <- vcov[b1_ref.pos, b3.pos]
      cov1_1_ref <- vcov[b1.pos, b1_ref.pos]
      cov3_3_ref <- vcov[b3.pos, b3_ref.pos]
      
      
      eff1 <- b1 - b1_ref + (b3 - b3_ref)* s2
      var.eff1 <- var1 + var1_ref  +  s2^2 * var3 + s2^2 * var3_ref - 2 * cov1_1_ref + 2 * s2 * cov13 - 2 * s2 * cov13_ref - 2 * s2 * cov1_ref_3 + 2 * s2 * cov1_ref_3_ref - 2 * s2^2 * cov3_3_ref 
                   
                  
    }
    
    se.eff1 <- sqrt(var.eff1)
    low1 <- eff1 - qt(1-alpha/2, df = obj$df.residual) * se.eff1
    up1 <- eff1 + qt(1-alpha/2, df = obj$df.residual) * se.eff1
    
    
    # relevant parameter values and standard deviations for SECOND marginal effects plot
    b2.pos <- grep(pattern = v2, x = names(coef(obj)), fixed = TRUE)[1]
    b2 <- as.numeric(obj$coef[b2.pos])
    var2 <- vcov[b2.pos, b2.pos]
    cov23 <- vcov[b2.pos, b3.pos]
    
    # calculate marginal effects and confidence band for SECOND marginal effects plot
    eff2 <- b2 + b3 * s1
    var.eff2 <- var2 + s1^2 * var3 + 2 * s1 * cov23
    se.eff2 <- sqrt(var.eff2)
    low2 <- eff2 - qt(1-alpha/2, df = obj$df.residual) * se.eff2
    up2 <- eff2 + qt(1-alpha/2, df = obj$df.residual) * se.eff2
    
    
    }
  
  if(class(obj) == "lmerMod" | class(obj) == "lmerModLmerTest"){
    
    if(!is.null(alt_ref)) print("Warning: So far, alternative reference category only implemented in lm!")
    
    b1.pos <- grep(v1, names(fixef(obj)))[1]
    b3.pos <- c(grep(pattern = paste(v1,":",v2, sep = ""), x = names(fixef(obj)), fixed = TRUE), 
                grep(pattern = paste(v2,":",v1, sep = ""), x = names(fixef(obj)), fixed = TRUE))[1]
    b1 <- as.numeric(fixef(obj)[b1.pos])
    b3 <- as.numeric(fixef(obj)[b3.pos])
    var1 <- vcov[b1.pos, b1.pos]
    var3 <- vcov[b3.pos, b3.pos]
    cov13 <- vcov[b1.pos, b3.pos]
    

    
    eff1 <- b1 + b3 * s2
    var.eff1 <- var1 + s2^2 * var3 + 2 * s2 * cov13
    
    if(!is.null(other_interactions)){
      positions_values <- lapply(other_interactions, 
                    FUN = function(v){
                      pos = c(grep(pattern = paste(v1,":",v, sep = ""), x = names(fixef(obj)), fixed = TRUE), 
                        grep(pattern = paste(v,":",v1, sep = ""), x = names(fixef(obj)), fixed = TRUE))[1]
                      if(length(unique(model.matrix(obj)[, paste(v)])) == 2){
                        value = 1
                      } else{
                        value = mean(model.matrix(obj)[, v])
                      }
                      return(c(pos,value))
                    }
                    )
      
      for(pos_value in positions_values){
        pos = pos_value[1]
        value = pos_value[2]
        eff1 <- eff1 + value * as.numeric(fixef(obj)[pos])
        var.eff1 <- var.eff1 + value^2 * vcov[pos,pos] + 2 * value * vcov[b1.pos, pos]  + 2 * s2 * value * vcov[b3.pos, pos]
        for(pos_value2 in positions_values){
          pos2 = pos_value2[1]
          value2 = pos_value2[2]
          if(pos < pos2){
            var.eff1 <- var.eff1 + 2 * value * value2 * vcov[pos,pos2]
          }
        }
      }
      
    }
    
    se.eff1 <- sqrt(var.eff1)
    low1 <- eff1 - qnorm(1-alpha/2) * se.eff1
    up1 <- eff1 + qnorm(1-alpha/2) * se.eff1
    
    b2.pos <- grep(v2, names(fixef(obj)))[1]
    b2 <- as.numeric(fixef(obj)[b2.pos])
    var2 <- vcov[b2.pos, b2.pos]

    cov23 <- vcov[b2.pos, b3.pos]
    
    eff2 <- b2 + b3 * s1
    var.eff2 <- var2 + s1^2 * var3 + 2 * s1 * cov23
    se.eff2 <- sqrt(var.eff2)
    low2 <- eff2 - qnorm(1-alpha/2) * se.eff2
    up2 <- eff2 + qnorm(1-alpha/2) * se.eff2
    
    
  }
  
  if(rug){
    # data points for rug at bottom of plots
    rug2.dat <- data.frame(x = model.matrix(obj)[, paste(v2)]) %>% 
      filter(x < max_x)
    
    if(color_rug){
      if(is.null(ref) & is.null(alt_ref) ){
        print("Warning: Either ref or alt_ref needs to be provided to color the rug. Proceeding without color information")
      } else {
        if(!is.null(alt_ref)){ # case: alternative reference category
          rug2.dat <- data.frame(x = model.matrix(obj)[, paste(v2)], 
                                 model.matrix(obj)[, paste(v1)],
                                 model.matrix(obj)[, paste(alt_ref)]
          ) %>% 
            mutate(category = factor(2 *model.matrix.obj....paste.v1.. + model.matrix.obj....paste.alt_ref..)) %>% 
            select(x, category)
        } else { # case: standard reference category
          if(is.null(data) | is.null(category_name)){
            print("Warning: dataset and category_name have to be provided to color the standard reference category. Proceeding without color information.")
          } else {
            rug2.dat <- data %>% 
              mutate(x = model.matrix(obj)[, paste(v2)], category = factor(case_match(!!sym(category_name), str_remove(v1, category_name) ~ 2, ref ~ 1, .default = 0))) %>% 
              select(x, category)
          }
        }
        rug2.dat <- rug2.dat %>% filter(x < max_x) %>% arrange(category)
      }
    }
  
    if(unlog_x){
      rug2.dat <- rug2.dat %>% mutate(x = exp(x))
    }
  
    rug1.dat <- data.frame(x = model.matrix(obj)[, paste(v1)])
  }
  
  
  if(unlog_x){
    s2 <-  exp(s2)
  }
  
  if(exp_y){
    eff1 <-  exp(eff1)
    low1 <- exp(low1)
    up1 <- exp(up1)
  }
  
  plot1.dat <- data.frame(s2, eff1, low1, up1)
  
  plot2.dat <- data.frame(s1, eff2, low2, up2)  
  
  
  # Generate the plot(s)

  if(length(unique(model.matrix(obj)[, paste(v2)])) == 2){
    
    p1 <- ggplot(data = plot1.dat, aes(x = factor(s2), y = eff1), environment = environment()) +     # env. important for rug
      geom_hline(yintercept = 0, color ="black", linetype = "dashed", size = 0.5) + 
      geom_segment(aes(x = factor(s2), xend = factor(s2), y = low1, yend = up1), color = "black") +
      geom_point() + 
      xlab(vlab2) + ylab(paste("Effect of ", vlab1, sep = ""))
  }
  
  if(length(unique(model.matrix(obj)[, paste(v2)])) > 2){
    
    p1 <- ggplot(data = plot1.dat, aes(x = s2, y = eff1), environment = environment()) +     # env. important for rug
      geom_segment(aes(y = ifelse(exp_y,1,0), yend = ifelse(exp_y,1,0), x = min(s2), xend = max(s2)), color = ifelse(color_lines, color_rug_ref, "black"), linetype = "dashed", size = 0.5) + 
      geom_ribbon(aes(x = s2, ymin = low1, ymax = up1), alpha = 0.25, fill = ifelse(color_band, hex(mixcolor(0.5, hex2RGB(color_rug_v1), hex2RGB(color_rug_ref))), "grey20")) +
      geom_line(color = ifelse(color_lines, color_rug_v1, "black")) + 
      xlab(vlab2) + ylab(paste("Effect of ", vlab1, sep = "")) + scale_x_continuous(label = comma)
  }
  
  if(title == TRUE){
    p1 <- p1 + labs(title = paste("Conditional effect of \n", vlab1, sep = ""),
                    subtitle = subtitle)
  }  
  
  if(rug == TRUE & length(unique(model.matrix(obj)[, paste(v2)])) > 2){
    if(length(unique(model.matrix(obj)[, paste(v1)])) == 2 & color_rug){
      p1 <- p1 + 
          geom_rug(data = rug2.dat, 
                   aes(x = jitter(x, factor = jitter_factor), y = 0, color = category), sides = "b", linewidth = rugwidth, alpha = rug_alpha, show.legend = legend) +
          scale_color_manual(breaks = c("2", "1"), values = c("2" = color_rug_v1, "1" = color_rug_ref, "0" = "grey"), labels = c(vlab1, paste0(ref_lab, " (reference)"))) #, "other"))
    } else{
      p1 <- p1 + 
        geom_rug(data = rug2.dat, 
                 aes(x = jitter(x, factor = jitter_factor), y = 0), sides = "b", linewidth = rugwidth, alpha = rug_alpha, show.legend = FALSE)
    }
  }     
  
  p1 <- p1 + theme_light() +
    theme(text = element_text(size=font_size, family=font_family), legend.position="top", legend.title = element_blank())
  
  if(log_scale_x_axis){
    p1 <- p1 + scale_x_continuous(trans='log10', label = comma)
  }
    
    
  # Make second plot if needed
  if (twoways == TRUE) {
    
    print("starting to make second plot")
    
    if(length(unique(model.matrix(obj)[, paste(v1)])) == 2){
      p2 <- ggplot(data = plot2.dat, aes(x = factor(s1), y = eff2), environment = environment()) +     # env. important for rug
        geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) + 
        geom_segment(aes(x = factor(s1), xend = factor(s1), y = low2, yend = up2), color = "black") +
        geom_point() + 
        xlab(paste(vlab1)) + ylab(paste("Effect of ", vlab2, sep = ""))
    }
    
    if(length(unique(model.matrix(obj)[, paste(v1)])) > 2){
      p2 <- ggplot(data = plot2.dat, aes(x = s1, y = eff2), environment = environment()) +     # env. important for rug
        geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) + 
        geom_ribbon(aes(x = s1, ymin = low2, ymax = up2), alpha = 0.25, color = NA) +
        geom_line() + 
        xlab(paste(vlab1)) + ylab(paste("Effect of ", vlab2, sep = "")) +
        scale_x_continuous(label = comma)
      
    }  
    
    p2 <- p2 + labs(title = paste("Conditional effect of \n", vlab2, sep = ""),
                    subtitle = subtitle) + theme_light() +
      theme(text = element_text(size=font_size, family=font_family))
    
    if(rug == TRUE & length(unique(model.matrix(obj)[, paste(v1)])) > 2){
      if(length(unique(model.matrix(obj)[, paste(v2)])) == 2 & color_rug){
        p2 <- p2 + 
          geom_rug(data = rug1.dat, 
                   aes(x = jitter(x, factor = jitter_factor), y = 0, color = category), sides = "b", linewidth = rugwidth, alpha = rug_alpha, show.legend = FALSE) +
          scale_color_manual(values = c("black", color_rug))
      }else{
        p2 <- p2 + 
          geom_rug(data = rug1.dat, 
                   aes(x = jitter(x, factor = jitter_factor), y = 0), sides = "b", linewidth = rugwidth, alpha = rug_alpha, show.legend = FALSE)
      }
    }     
    
    
    if(log_scale_x_axis){
      p1 <- p1 + scale_x_continuous(trans='log', label = comma)
    }
    
    grid.arrange(p1, p2, ncol = 2)
  } else {
    return(p1)
  }
}





ggintplot_categorical_moderator <- function (obj, pred, categorical_moderator, levels, vcov = NULL, varlabs, other_interactions = NULL,
                       alpha = 0.05,
                       title = FALSE,
                       subtitle = NULL, 
                       font_family = "sans",
                       font_size = 7, verbose = FALSE){
  
  require(ggplot2); require(gridExtra)
  
  if(is.null(vcov)) vcov = vcov(obj)
  
  if(class(obj) == "lmerMod"){
    
    if(!is.null(other_interactions)) print("Warning: so far only one interaction term per effect")
    
    

    # build the names of the dummy variables (first one needs none because it is the reference level)
    modnames <- lapply(levels[-1], function(level) paste(categorical_moderator, level, sep=""))
    
    
    
    # the effect of the reference level is simply the one reported for the predictor 
    eff_reference <- as.numeric(fixef(obj)[pred])
    var_reference <- vcov[pred,pred]
    se_reference <- sqrt(var_reference)
    
    # for each other level, we have to add the interaction coefficient
    effects <- 
      lapply(modnames,
             function(mod){
               pos = c(grep(pattern = paste(mod,":",pred, sep = ""), x = names(fixef(obj)), fixed = TRUE), 
                       grep(pattern = paste(pred,":",mod, sep = ""), x = names(fixef(obj)), fixed = TRUE))[1]
               eff_reference + as.numeric(fixef(obj)[pos])
             })
    
    # same for standard errors
    standard_errors <- 
      lapply(modnames,
             function(mod){
               pos = c(grep(pattern = paste(mod,":",pred, sep = ""), x = names(fixef(obj)), fixed = TRUE), 
                       grep(pattern = paste(pred,":",mod, sep = ""), x = names(fixef(obj)), fixed = TRUE))[1]
               var <- var_reference + vcov[pos, pos] + 2 * vcov[pred, pos]
               if(var < 0) {
                 print("warning: Negative variance for group")
                 print(mod)
               }
               if(verbose){
                 print(mod)
                 print(pos)
                 print(var)
               }
               return(sqrt(var))
             })
    
    
    
    effects <- append(effects, eff_reference, after = 0)
    standard_errors <- append(standard_errors, se_reference, after = 0)
    
    plot.dat <- 
      tibble(moderator = unlist(levels), effect = unlist(effects), standard_error = unlist(standard_errors)) %>% 
      mutate(low = effect - qnorm(1-alpha/2) * standard_error,
             up = effect + qnorm(1-alpha/2) * standard_error
      )
    
    
    
    
  }
  
  
  if(class(obj) == "lm"){
    
    if(!is.null(other_interactions)) print("Warning: so far only one interaction term per effect")
    
    if (!("model" %in% names(obj))) {
      obj <- update(obj, model = T)
    }
    

    # build the names of the dummy variables (first one needs none because it is the reference level)
    modnames <- lapply(levels[-1], function(level) paste(categorical_moderator, level, sep=""))
    
    # the effect of the reference level is simply the one reported for the predictor 
    eff_reference <- as.numeric(obj$coef[pred])
    var_reference <- vcov[pred,pred]
    se_reference <- sqrt(var_reference)
    
    # for each other level, we have to add the interaction coefficient
    effects <- 
      lapply(modnames,
             function(mod){
               pos = c(grep(pattern = paste(mod,":",pred, sep = ""), x = names(coef(obj)), fixed = TRUE), 
                       grep(pattern = paste(pred,":",mod, sep = ""), x = names(coef(obj)), fixed = TRUE))[1]
               eff_reference + as.numeric(obj$coef[pos])
    })
    
    # same for standard errors
    standard_errors <- 
      lapply(modnames,
             function(mod, verbose){
               pos = c(grep(pattern = paste(mod,":",pred, sep = ""), x = names(coef(obj)), fixed = TRUE), 
                       grep(pattern = paste(pred,":",mod, sep = ""), x = names(coef(obj)), fixed = TRUE))[1]
               var <- var_reference + vcov[pos, pos] + 2 * vcov[pred, pos]
               if(var < 0) {
                 print("warning: Negative variance for group")
                 print(mod)
               }
               if(verbose){
                 print(mod)
                 print(pos)
                 print(var)
                 }
               return(sqrt(var))
             },
             verbose)
    
    
    
    effects <- append(effects, eff_reference, after = 0)
    standard_errors <- append(standard_errors, se_reference, after = 0)
    
    plot.dat <- 
      tibble(moderator = unlist(levels), effect = unlist(effects), standard_error = unlist(standard_errors)) %>% 
      mutate(low = effect - qt(1-alpha/2, df = obj$df.residual) * standard_error,
             up = effect + qt(1-alpha/2, df = obj$df.residual) * standard_error
      )
      
    if(verbose) print(plot.dat)
    
    
  }
  
    p1 <- ggplot(data = plot.dat, aes(x = forcats::fct_inorder(moderator), y = effect)) +   
      geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.5) + 
      geom_segment(aes(x = forcats::fct_inorder(moderator), xend = moderator, y = low, yend = up), color = "black") +
      geom_point() +
      xlab(sub("_", " ", categorical_moderator)) + 
      ylab(paste("Effect of ", sub("_", " ", pred), sep = "")) +
      scale_x_discrete(labels = sub("_", " ", levels))
 
    

    
    if(title == TRUE){
      p1 <- p1 + labs(title = paste("Conditional effect of \n", sub("_", " ", pred), sep = ""),
                      subtitle = subtitle)
    }  
    
    p1 <- p1 + theme_light() +
      theme(text = element_text(size=font_size, family=font_family))
    

    
    return(p1)
  
  
  
  
  
  
}




get_insignificant_groups <- function (obj, pred, categorical_moderator, levels, vcov = NULL, other_interactions = NULL,
                                      alpha = 0.05, negative_variance_insignificant = TRUE){
  
  
  if(is.null(vcov)) vcov = vcov(obj)
  
  if(class(obj) == "lmerMod"){
    
    if(!is.null(other_interactions)) print("Warning: so far only one interaction term per effect")
    
    # find name of predictor including transformations
    pred <- str_subset(names(obj$coef), pattern = paste("^(?!.*:).*",pred,".*", sep = ""))[[1]]
    
    # build the names of the dummy variables (first one needs none because it is the reference level)
    modnames <- lapply(levels[-1], function(level) paste(categorical_moderator, level, sep=""))
    
    
    
    # the effect of the reference level is simply the one reported for the predictor 
    eff_reference <- as.numeric(fixef(obj)[pred])
    var_reference <- vcov[pred,pred]
    se_reference <- sqrt(var_reference)
    
    # for each other level, we have to add the interaction coefficient
    effects <- 
      lapply(modnames,
             function(mod){
               pos = c(grep(pattern = paste(mod,":",pred, sep = ""), x = names(fixef(obj)), fixed = TRUE), 
                       grep(pattern = paste(pred,":",mod, sep = ""), x = names(fixef(obj)), fixed = TRUE))[1]
               eff_reference + as.numeric(fixef(obj)[pos])
             })
    
    # same for standard errors
    standard_errors <- 
      lapply(modnames,
             function(mod){
               pos = c(grep(pattern = paste(mod,":",pred, sep = ""), x = names(fixef(obj)), fixed = TRUE), 
                       grep(pattern = paste(pred,":",mod, sep = ""), x = names(fixef(obj)), fixed = TRUE))[1]
               var <- var_reference + vcov[pos, pos] + 2 * vcov[pred, pos]
               if(var < 0){
                 print("warning: Negative variance for group")
                 print(mod)
               }
               return(sqrt(var))
             })
    
    
    
    effects <- append(effects, eff_reference, after = 0)
    standard_errors <- append(standard_errors, se_reference, after = 0)
    
    p_values <- 
      tibble(level = unlist(levels), effect = unlist(effects), standard_error = unlist(standard_errors)) %>% 
      mutate(p = 2*(1-pnorm(abs(effect/ standard_error))))
      

    
    
    
  }
  
  
  if(class(obj) == "lm"){
    
    if(!is.null(other_interactions)) print("Warning: so far only one interaction term per effect")
    
    if (!("model" %in% names(obj))) {
      obj <- update(obj, model = T)
    }
    
    
    # find name of predictor including transformations
    pred <- str_subset(names(obj$coef), pattern = paste("^(?!.*:).*",pred,".*", sep = ""))[[1]]
    
    
    
    # build the names of the dummy variables (first one needs none because it is the reference level)
    modnames <- lapply(levels[-1], function(level) paste(categorical_moderator, level, sep=""))
    
    # the effect of the reference level is simply the one reported for the predictor 
    eff_reference <- as.numeric(obj$coef[pred])
    var_reference <- vcov[pred,pred]
    se_reference <- sqrt(var_reference)
    
    # for each other level, we have to add the interaction coefficient
    effects <- 
      lapply(modnames,
             function(mod){
               pos = c(grep(pattern = paste(mod,":",pred, sep = ""), x = names(coef(obj)), fixed = TRUE), 
                       grep(pattern = paste(pred,":",mod, sep = ""), x = names(coef(obj)), fixed = TRUE))[1]
               eff_reference + as.numeric(obj$coef[pos])
             })
    
    # same for standard errors
    standard_errors <- 
      lapply(modnames,
             function(mod){
               pos = c(grep(pattern = paste(mod,":",pred, sep = ""), x = names(coef(obj)), fixed = TRUE), 
                       grep(pattern = paste(pred,":",mod, sep = ""), x = names(coef(obj)), fixed = TRUE))[1]
               var <- var_reference + vcov[pos, pos] + 2 * vcov[pred, pos]
               if(var < 0){
                 print("warning: Negative variance for group")
                 print(mod)
               }
               return(sqrt(var))             
               })
    
    
    
    effects <- append(effects, eff_reference, after = 0)
    standard_errors <- append(standard_errors, se_reference, after = 0)
    
    p_values <- 
      tibble(level = unlist(levels), effect = unlist(effects), standard_error = unlist(standard_errors)) %>% 
      mutate(p = 2*(1-pt(abs(effect/ standard_error), df = obj$df.residual)))
    

    
    
    
    
  }
  
  insignificant_groups <- p_values %>% filter(p>alpha | (negative_variance_insignificant & is.na(p))) %>% pull(level)
  
 return(insignificant_groups)
  
  
}


insertLayer <- function(P, after=0, ...) {
  #  P     : Plot object
  # after  : Position where to insert new layers, relative to existing layers
  #  ...   : additional layers, separated by commas (,) instead of plus sign (+)
  
  if (after < 0)
    after <- after + length(P$layers)
  
  if (!length(P$layers))
    P$layers <- list(...)
  else 
    P$layers <- append(P$layers, list(...), after)
  
  return(P)
}



#' Plot interaction effects in regression models
#'
#' \code{interact_plot} plots regression lines at user-specified levels of a
#'  moderator variable to explore interactions. The plotting is done with
#'  \code{ggplot2} rather than base graphics, which some similar functions use.
#'
#' @param model A regression model. The function is tested with \code{lm},
#'   \code{glm}, \code{\link[survey]{svyglm}}, \code{\link[lme4]{merMod}},
#'   \code{\link[quantreg]{rq}}, \code{\link[brms]{brmsfit}},
#'   \code{stanreg} models.
#'   Models from other classes may work as well but are not officially
#'   supported. The model should include the interaction of interest.
#'
#' @param pred The name of the predictor variable involved
#'  in the interaction. This can be a bare name or string. Note that it
#'  is evaluated using `rlang`, so programmers can use the `!!` syntax
#'  to pass variables instead of the verbatim names.
#'
#' @param modx The name of the moderator variable involved
#'  in the interaction. This can be a bare name or string. The same
#'  `rlang` proviso applies as with `pred`.
#'
#' @param mod2 Optional. The name of the second moderator
#'  variable involved in the interaction. This can be a bare name or string.
#'  The same `rlang` proviso applies as with `pred`.
#'
#' @param modx.values For which values of the moderator should lines be
#'   plotted? There are two basic options:
#'
#'   * A vector of values (e.g., `c(1, 2, 3)`)
#'   * A single argument asking to calculate a set of values. See details
#'   below.
#'
#'   Default is \code{NULL}. If \code{NULL} (or `mean-plus-minus`),
#'   then the customary +/- 1 standard
#'   deviation from the mean as well as the mean itself are used for continuous
#'   moderators. If \code{"plus-minus"}, plots lines when the moderator is at
#'   +/- 1 standard deviation without the mean. You may also choose `"terciles"`
#'   to split the data into equally-sized groups and choose the point at the
#'   mean of each of those groups.
#'
#'   If the moderator is a factor variable and \code{modx.values} is
#'   \code{NULL}, each level of the factor is included. You may specify
#'   any subset of the factor levels (e.g., `c("Level 1", "Level 3")`) as long
#'   as there is more than 1. The levels will be plotted in the order you
#'   provide them, so this can be used to reorder levels as well.
#'
#' @param mod2.values For which values of the second moderator should the plot
#'   be
#'   facetted by? That is, there will be a separate plot for each level of this
#'   moderator. Defaults are the same as \code{modx.values}.
#'
#' @param centered A vector of quoted variable names that are to be
#'   mean-centered. If `"all"`, all non-focal predictors are centered. You
#'   may instead pass a character vector of variables to center. User can
#'   also use "none" to base all predictions on variables set at 0.
#'   The response variable, `pred`, `modx`, and `mod2` variables are never
#'   centered.
#'
#' @param data Optional, default is NULL. You may provide the data used to
#'   fit the model. This can be a better way to get mean values for centering
#'   and can be crucial for models with variable transformations in the formula
#'   (e.g., `log(x)`) or polynomial terms (e.g., `poly(x, 2)`). You will
#'   see a warning if the function detects problems that would likely be
#'   solved by providing the data with this argument and the function will
#'   attempt to retrieve the original data from the global environment.
#'
#' @param plot.points Logical. If \code{TRUE}, plots the actual data points as
#'   a scatterplot on top of the interaction lines. The color of the dots will
#'   be based on their moderator value.
#'
#' @param interval Logical. If \code{TRUE}, plots confidence/prediction
#'   intervals around the line using \code{\link[ggplot2]{geom_ribbon}}.
#'
#' @param int.type Type of interval to plot. Options are "confidence" or
#'  "prediction". Default is confidence interval.
#'
#' @param int.width How large should the interval be, relative to the standard
#'   error? The default, .95, corresponds to roughly 1.96 standard errors and
#'   a .05 alpha level for values outside the range. In other words, for a
#'   confidence interval, .95 is analogous to a 95% confidence interval.
#'
#' @param outcome.scale For nonlinear models (i.e., GLMs), should the outcome
#'   variable be plotted on the link scale (e.g., log odds for logit models) or
#'   the original scale (e.g., predicted probabilities for logit models)? The
#'   default is \code{"response"}, which is the original scale. For the link
#'   scale, which will show straight lines rather than curves, use
#'   \code{"link"}.
#'
#' @param linearity.check For two-way interactions only. If `TRUE`, plots a
#'   pane for each level of the moderator and superimposes a loess smoothed
#'   line (in gray) over the plot. This enables you to see if the effect is
#'   linear through the span of the moderator. See Hainmueller et al. (2016) in
#'   the references for more details on the intuition behind this. It is
#'   recommended that you also set `plot.points = TRUE` and use
#'   `modx.values = "terciles"` with this option.
#'
#' @inheritParams jtools::summ.lm
#'
#' @param vcov Optional. You may supply the variance-covariance matrix of the
#'  coefficients yourself. This is useful if you are using some method for
#'  robust standard error calculation not supported by the \pkg{sandwich}
#'  package.
#'
#' @param set.offset For models with an offset (e.g., Poisson models), sets an
#'   offset for the predicted values. All predicted values will have the same
#'   offset. By default, this is set to 1, which makes the predicted values a
#'   proportion. See details for more about offset support.
#'
#' @param x.label A character object specifying the desired x-axis label. If
#'   \code{NULL}, the variable name is used.
#'
#' @param y.label A character object specifying the desired x-axis label. If
#'   \code{NULL}, the variable name is used.
#'
#' @param pred.labels A character vector of 2 labels for the predictor if it is
#'   a 2-level factor or a continuous variable with only 2 values. If
#'   \code{NULL}, the default, the factor labels are used.
#'
#' @param modx.labels A character vector of labels for each level of the
#'   moderator values, provided in the same order as the \code{modx.values}
#'   argument. If \code{NULL}, the values themselves are used as labels unless
#'   \code{modx,values} is also \code{NULL}. In that case, "+1 SD" and "-1 SD"
#'   are used.
#'
#' @param mod2.labels A character vector of labels for each level of the 2nd
#'   moderator values, provided in the same order as the \code{mod2.values}
#'   argument. If \code{NULL}, the values themselves are used as labels unless
#'   \code{mod2.values} is also \code{NULL}. In that case, "+1 SD" and "-1 SD"
#'   are used.
#'
#' @param main.title A character object that will be used as an overall title
#'   for the plot. If \code{NULL}, no main title is used.
#'
#' @param legend.main A character object that will be used as the title that
#'   appears above the legend. If \code{NULL}, the name of the moderating
#'   variable is used.
#'
#' @param colors See [jtools_colors] for details on the types of arguments
#'    accepted. Default is "CUD Bright" for factor
#'    moderators, "Blues" for +/- SD and user-specified \code{modx.values}
#'    values.
#'
#' @param line.thickness How thick should the plotted lines be? Default is 1.
#'
#' @param vary.lty Should the resulting plot have different shapes for each
#'   line in addition to colors? Defaults to \code{TRUE}.
#'
#' @param jitter How much should `plot.points` observed values be "jittered"
#'    via [ggplot2::position_jitter()]? When there are many points near each
#'    other, jittering moves them a small amount to keep them from
#'    totally overlapping. In some cases, though, it can add confusion since
#'    it may make points appear to be outside the boundaries of observed
#'    values or cause other visual issues. Default is 0, but try various
#'    small values (e.g., 0.1) and increase as needed if your points are
#'    overlapping too much. If the argument is a vector with two values,
#'    then the first is assumed to be the jitter for width and the second
#'    for the height.
#'
#' @param rug Show a rug plot in the margins? This uses [ggplot2::geom_rug()]
#'    to show the distribution of the predictor (top/bottom) and/or
#'    response variable (left/right) in the original data. Default is
#'    FALSE.
#'
#' @param rug.sides On which sides should rug plots appear? Default is "b",
#'    meaning bottom. "t" and/or "b" show the distribution of the predictor
#'    while "l" and/or "r" show the distribution of the response. "bl" is
#'    a good option to show both the predictor and response.
#'
#' @param point.size What size should be used for observed data when
#'   `plot.points` is TRUE? Default is 1.5.
#'
#' @param facet.modx Create separate panels for each level of the moderator?
#'   Default is FALSE, except when `linearity.check` is TRUE.
#'
#' @param robust Should robust standard errors be used to find confidence
#'   intervals for supported models? Default is FALSE, but you should specify
#'   the type of sandwich standard errors if you'd like to use them (i.e.,
#'   `"HC0"`, `"HC1"`, and so on). If `TRUE`, defaults to `"HC3"` standard
#'   errors.
#'
#' @param cluster For clustered standard errors, provide the column name of
#'   the cluster variable in the input data frame (as a string). Alternately,
#'   provide a vector of clusters.
#'
#' @param ... extra arguments passed to `make_predictions`
#'
#' @inheritParams cat_plot
#' @inheritParams jtools::effect_plot
#'
#' @details This function provides a means for plotting conditional effects
#'   for the purpose of exploring interactions in regression models.
#'
#'   The function is designed for two and three-way interactions. For
#'   additional terms, the \pkg{effects} package may be better suited to the
#'   task.
#'
#'   This function supports nonlinear and generalized linear models and by
#'   default will plot them on their original scale
#'   (`outcome.scale = "response"`). To plot them on the linear scale,
#'   use "link" for `outcome.scale`.
#'
#'   While mixed effects models from \code{lme4} are supported, only the fixed
#'   effects are plotted. \code{lme4} does not provide confidence intervals,
#'   so they are not supported with this function either.
#'
#'   Note: to use transformed predictors, e.g., \code{log(variable)},
#'   put its name in quotes or backticks in the argument.
#'
#'   \emph{Details on how observed data are split in multi-pane plots}:
#'
#'   If you set `plot.points = TRUE` and request a multi-pane (facetted) plot
#'   either with a second moderator, `linearity.check = TRUE`, or
#'   `facet.modx = TRUE`, the observed
#'   data are split into as many groups as there  are panes and plotted
#'   separately. If the moderator is a factor, then the way this happens will
#'   be very intuitive since it's obvious which values go in which pane. The
#'   rest of this section will address the case of continuous moderators.
#'
#'   My recommendation is that you use `modx.values = "terciles"` or
#'   `mod2.values = "terciles"` when you want to plot observed data on
#'   multi-pane
#'   plots. When you do, the data are split into three approximately
#'   equal-sized groups with the lowest third, middle third, and highest third
#'   of the data split accordingly. You can replicate this procedure using
#'   [Hmisc::cut2()] with `g = 3` from the `Hmisc` package. Sometimes, the
#'   groups will not be equal in size because the number of observations is
#'   not divisible by 3 and/or there are multiple observations with the same
#'   value at one of the cut points.
#'
#'   Otherwise, a more ad hoc procedure is used to split the data. Quantiles
#'   are found for each `mod2.values` or `modx.values` value. These are not the
#'   quantiles used to split the data, however, since we want the plotted lines
#'   to represent the slope at a typical value in the group. The next step,
#'   then, is to take the mean of each pair of neighboring quantiles and use
#'   these as the cut points.
#'
#'   For example, if the `mod2.values` are at the 25th, 50th, and 75th
#'   percentiles
#'   of the distribution of the moderator, the data will be split at the
#'   37.5th and and 62.5th percentiles. When the variable is
#'   normally distributed, this will correspond fairly closely to using
#'   terciles.
#'
#'   \emph{Info about offsets:}
#'
#'   Offsets are partially supported by this function with important
#'   limitations. First of all, only a single offset per model is supported.
#'   Second, it is best in general to specify offsets with the offset argument
#'   of the model fitting function rather than in the formula. You are much
#'   more likely to have success if you provide the data used to fit the model
#'   with the `data` argument.
#'
#'
#' @return The functions returns a \code{ggplot} object, which can be treated
#'   like a user-created plot and expanded upon as such.
#'
#' @author Jacob Long \email{jacob.long@@sc.edu}
#'
#' @seealso \code{\link[rockchalk]{plotSlopes}} from \pkg{rockchalk} performs a
#'   similar function, but
#'   with R's base graphics---this function is meant, in part, to emulate
#'   its features.
#'
#'   \code{\link{sim_slopes}} performs a simple slopes analysis with a similar
#'   argument syntax to this function.
#'
#' @references
#'
#' Bauer, D. J., & Curran, P. J. (2005). Probing interactions in fixed and
#'  multilevel regression: Inferential and graphical techniques.
#'  \emph{Multivariate Behavioral
#'  Research}, \emph{40}(3), 373-400.
#'  \doi{10.1207/s15327906mbr4003_5}
#'
#' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied
#' multiple
#' regression/correlation analyses for the behavioral sciences} (3rd ed.).
#' Mahwah, NJ: Lawrence Erlbaum Associates, Inc.
#'
#' Hainmueller, J., Mummolo, J., & Xu, Y. (2016). How much should we trust
#'   estimates from multiplicative interaction models? Simple tools to improve
#'   empirical practice. SSRN Electronic Journal.
#'   \doi{10.2139/ssrn.2739221}
#'
#' @examples
#' # Using a fitted lm model
#' states <- as.data.frame(state.x77)
#' states$HSGrad <- states$`HS Grad`
#' fit <- lm(Income ~ HSGrad + Murder * Illiteracy, data = states)
#' interact_plot(model = fit, pred = Murder, modx = Illiteracy)
#'
#' # Using interval feature
#' fit <- lm(accel ~ mag * dist, data = attenu)
#' interact_plot(fit, pred = mag, modx = dist, interval = TRUE,
#'   int.type = "confidence", int.width = .8)
#'
#' # Using second moderator
#' fit <- lm(Income ~ HSGrad * Murder * Illiteracy, data = states)
#' interact_plot(model = fit, pred = Murder, modx = Illiteracy, mod2 = HSGrad)
#'
#' # With svyglm
#' if (requireNamespace("survey")) {
#' library(survey)
#' data(api)
#' dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw,
#'                     data = apistrat, fpc = ~fpc)
#' regmodel <- svyglm(api00 ~ ell * meals, design = dstrat)
#' interact_plot(regmodel, pred = ell, modx = meals)
#' }
#'
#' # With lme4
#' \dontrun{
#' library(lme4)
#' data(VerbAgg)
#' mv <- glmer(r2 ~ Anger * mode + (1 | item), data = VerbAgg,
#'             family = binomial,
#'             control = glmerControl("bobyqa"))
#' interact_plot(mv, pred = Anger, modx = mode)
#' }
#'
#' @importFrom stats coef coefficients lm predict sd qnorm getCall model.offset
#' @importFrom stats median ecdf quantile
#' @import ggplot2
#' @import rlang
#' @export interact_plot

interact_plot <- function(model, pred, modx, resp = NULL, saturation_level = NULL, reverse = FALSE, modx.values = NULL, mod2 = NULL,
                          mod2.values = NULL, centered = "all", data = NULL,
                          at = NULL,
                          plot.points = FALSE, interval = FALSE,
                          int.type = c("confidence", "prediction"),
                          int.width = .95, outcome.scale = "response",
                          linearity.check = FALSE, facet.modx = FALSE,
                          robust = FALSE, cluster = NULL, vcov = NULL,
                          set.offset = 1,
                          x.label = NULL, y.label = NULL,
                          pred.labels = NULL, modx.labels = NULL,
                          mod2.labels = NULL, main.title = NULL,
                          legend.main = NULL, colors = NULL,
                          line.thickness = 1, vary.lty = TRUE,
                          point.size = 1.5, point.shape = FALSE,
                          jitter = 0, rug = FALSE, rug.sides = "b",
                          partial.residuals = FALSE, point.alpha = 0.6,
                          color.class = NULL,
                          restrict_lines_data = FALSE,
                          insignificant_lines = "normal",
                          show_guide_line = TRUE,
                          trans_x = NA,
                          x_max = NA,
                          slack_factor = 0.005,
                          verbose = FALSE,  ...) {
  
  
  # Capture extra arguments
  dots <- list(...)
  if (length(dots) > 0) {
    if ("modxvals" %in% names(dots)) {
      modx.values <- dots$modxvals
    }
    if ("mod2vals" %in% names(dots)) {
      mod2.values <- dots$mod2vals
    }
    # If it's a categorical predictor, I want a different default for the
    # geom argument than cat_plot() uses so it looks more like what you'd 
    # expect from this function
    if ("geom" %nin% names(dots)) {
      geom <- "line"
    } else {
      geom <- dots$geom
    }
  } else {
    geom <- "line"
  }
  
  if (!is.null(color.class)) {
    colors <- color.class
    msg_wrap("The color.class argument is deprecated. Please use 'colors'
             instead.")
  }
  
  # Evaluate the modx, mod2, pred args
  pred <- as_name(enquo(pred))
  modx <- enquo(modx)
  modx <- if (quo_is_null(modx)) {NULL} else {as_name(modx)}
  mod2 <- enquo(mod2)
  mod2 <- if (quo_is_null(mod2)) {NULL} else {as_name(mod2)}
  
  if (any(c(pred, modx, mod2) %in% centered)) {
    warn_wrap("You cannot mean-center the focal predictor or moderators with
              this function.")
    centered <- centered %not% c(pred, modx, mod2)
    if (length(centered) == 0) {centered <- "none"}
  }
  
  # Defining "global variables" for CRAN
  modxvals2 <- mod2vals2 <- NULL
  
  # Pulling the name of the response variable 
  if(is.null(resp)) resp <- jtools::get_response_name(model, ...)
  
  if(verbose) print(resp)
  
  # Change facet.modx to TRUE if linearity.check is TRUE
  if (linearity.check == TRUE) {facet.modx <- TRUE}
  
  if (is.null(data)) {
    d <- get_data(model, warn = TRUE, ...)
  } else {
    d <- data
  }
  weights <- get_weights(model, d)$weights_name
  
  # Check for variables in the data
  if (any(c(pred, modx, mod2) %nin% names(d))) {
    missed_vars <- c(pred, modx, mod2) %not% names(d)
    stop_wrap(paste(missed_vars, collapse = " and "),
              ifelse(length(missed_vars) > 1, yes = " were ", no = " was "),
              "not found in the data. If you are using a transformed variable,
              like 'log(x)', use the non-transformed name ('x') as the input to
              this function.")
  }
  
  # If modx.values is named, use the names as labels
  if (is.null(modx.labels) & !is.null(names(modx.values))) {
    modx.labels <- names(modx.values)
  }
  # If mod2.values is named, use the names as labels
  if (is.null(mod2.labels) & !is.null(names(mod2.values))) {
    mod2.labels <- names(mod2.values)
  }
  
  pred_out <- prep_data(model = model, pred = pred, resp = resp, saturation_level = saturation_level, reverse = reverse, modx = modx,
                        modx.values = modx.values, mod2 = mod2,
                        mod2.values = mod2.values, centered = centered,
                        interval = interval, int.type = int.type,
                        int.width = int.width, outcome.scale = outcome.scale,
                        linearity.check = linearity.check, robust = robust,
                        cluster = cluster, vcov = vcov, set.offset = set.offset,
                        modx.labels = modx.labels, mod2.labels = mod2.labels,
                        facet.modx = facet.modx, d = d,
                        survey = "svyglm" %in% class(model), weights = weights,
                        preds.per.level = 500,
                        trans_x = trans_x,
                        x_max = x_max,
                        slack_factor = slack_factor,
                        partial.residuals = partial.residuals, at = at, restrict_lines_data = restrict_lines_data, insignificant_lines = insignificant_lines, verbose = verbose, ...)
  
  # These are the variables created in the helper functions
  meta <- attributes(pred_out)
  # This function attaches all those variables to this environment
  lapply(names(meta), function(x, env) {env[[x]] <- meta[[x]]},
         env = environment())
  
  # Putting these outputs into separate objects
  pm <- pred_out$predicted
  d <- pred_out$original
  
  # Check for factor predictor and send to plot_cat() if so
  if (!is.numeric(d[[pred]])) {
    # Warn users that this is kinda janky
    msg_wrap("Detected factor predictor. Plotting with cat_plot() instead.
              Consult cat_plot() documentation (?cat_plot) for full details
              on how to specify models with categorical predictors. If you 
              experience errors or unexpected results, try using cat_plot() 
              directly.")
    # Gather arguments for plot_cat()
    args <- list(predictions = pm, pred = pred, modx = modx, mod2 = mod2,
                 data = d, modx.values = modxvals2, mod2.values = mod2vals2,
                 interval = interval,
                 plot.points = plot.points | partial.residuals,
                 point.shape = point.shape, vary.lty = vary.lty,
                 pred.labels = pred.labels, modx.labels = modx.labels,
                 mod2.labels = mod2.labels, x.label = x.label, y.label = y.label,
                 main.title = main.title, legend.main = legend.main,
                 colors = colors, weights = weights, resp = resp,
                 point.size = point.size, line.thickness = line.thickness,
                 jitter = jitter, point.alpha = point.alpha, geom = geom)
    # Deal with cat_plot() arguments provided via ...
    if (length(dots) > 0) {
      # Make sure it's not geom which I've already handled
      if (length(dots %not% "geom") > 0) {
        # Append to this list
        args <- c(args, dots %not% "geom")
      }
    }
    # Call internal plotting function
    return(do.call("plot_cat", args))
    # Using plot_cat turns out to be more robust than cat_plot. I'm doing 
    # something dumb and/or badly with environments that makes calling 
    # plot_cat() within a function fail inside of testthat (and inside of
    # any arbitrary function that calls interact_plot() with a categorical
    # predictor even in the "normal" environment).
  } else {
    # Send to internal plotting function
    plot_mod_continuous(model = model, predictions = pm, pred = pred, modx = modx, modx.values = modx.values, resp = resp,
                        mod2 = mod2, data = d, vcov = vcov,
                        plot.points = plot.points | partial.residuals,
                        interval = interval, linearity.check = linearity.check,
                        x.label = x.label, y.label = y.label,
                        pred.labels = pred.labels, modx.labels = modx.labels,
                        mod2.labels = mod2.labels, main.title = main.title,
                        legend.main = legend.main, colors = colors,
                        line.thickness = line.thickness,
                        vary.lty = vary.lty, jitter = jitter,
                        modxvals2 = modxvals2, mod2vals2 = mod2vals2,
                        weights = weights, rug = rug, rug.sides = rug.sides,
                        point.size = point.size, point.shape = point.shape,
                        facet.modx = facet.modx, point.alpha = point.alpha,
                        insignificant_lines = insignificant_lines, show_guide_line = show_guide_line,
                        verbose = verbose)
  }
  
}

# Workhorse plotting function
plot_mod_continuous <- function(model, predictions, pred, modx, modx.values, resp, mod2 = NULL,
                                data = NULL, vcov = NULL,
                                plot.points = FALSE,
                                interval = FALSE, linearity.check = FALSE,
                                x.label = NULL, y.label = NULL,
                                pred.labels = NULL, modx.labels = NULL,
                                mod2.labels = NULL, main.title = NULL,
                                legend.main = NULL, colors = NULL,
                                line.thickness = 1.1, vary.lty = TRUE,
                                jitter = 0, modxvals2 = NULL,
                                mod2vals2 = NULL, weights = NULL, rug = FALSE,
                                rug.sides = "b",
                                point.shape = FALSE, point.size = 2,
                                facet.modx = FALSE, point.alpha = 0.6,
                                insignificant_lines="normal", 
                                show_guide_line = TRUE,
                                verbose = FALSE) {
  
  if(verbose) print(resp)
  
  d <- data
  pm <- predictions
  
  # Setting default for colors
  if (is.null(colors) && (facet.modx == TRUE | linearity.check == TRUE)) {
    colors <- rep("black", times = length(modxvals2))
    vary.lty <- FALSE
    point.shape <- FALSE
  }
  if (is.factor(d[[modx]])) {
    facmod <- TRUE
    if (is.null(colors)) {
      colors <- "CUD Bright"
    }
  } else {
    facmod <- FALSE
    if (is.null(colors)) {
      colors <- "blue"
    }
  }
  
  # If only 1 jitter arg, just duplicate it
  if (length(jitter) == 1) {jitter <- rep(jitter, 2)}
  
  # If no user-supplied legend title, set it to name of moderator
  if (is.null(legend.main)) {
    legend.main <- modx
  }
  
  if (is.null(x.label)) {
    x.label <- pred
  }
  
  if (is.null(y.label)) {
    y.label <- resp
  }
  
  if (is.null(modxvals2)) {
    modxvals2 <- unique(pm[[modx]])
  }
  
  if (!is.null(mod2) && is.null(mod2vals2)) {
    mod2vals2 <- unique(pm[[mod2]])
  }
  
  gradient <- is.numeric(d[[modx]])
  
  # Checking if user provided the colors his/herself
  colors <- suppressWarnings(get_colors(colors, length(modx.labels),
                                        gradient = gradient))
  
  # Manually set linetypes
  types <- c("solid", "4242", "2222", "dotdash", "dotted", "twodash",
             "12223242", "F282", "F4448444", "224282F2", "F1")
  ltypes <- types[seq_along(modxvals2)]
  
  # Reverse the order of the linetypes to make thick line go to biggest value
  if (is.numeric(modxvals2) & all(sort(modxvals2) == modxvals2)) {
    ltypes <- rev(ltypes)
  } else if (!is.null(mod2) & !(is.numeric(modxvals2) & !all(sort(modxvals2) == modxvals2))) { # also flip for factor second moderators
    ltypes <- rev(ltypes)
  }
  
  if (gradient == FALSE) {
    names(colors) <- modx.labels
  }
  names(ltypes) <- modx.labels
  
  if(verbose) print(resp)
  
  # Prepare names for tidy evaluation
  pred <- sym(pred)
  resp <- sym(resp)
  if (!is.null(modx)) {modx <- sym(modx)}
  if (!is.null(mod2)) {mod2 <- sym(mod2)}
  if (!is.null(weights)) {weights <- sym(weights)}
  
  if(verbose) print(pm)
  
  lty <- if (vary.lty) sym("modx_group") else{ if (insignificant_lines == "dashed") sym("significance_line") else NULL}
  if(verbose) print(lty)
  # Don't use 'modx_group' if I don't have to since it makes it harder for
  # users to make changes after the fact
  grp <- if (vary.lty | facet.modx) sym("modx_group") else modx
  
  p <- ggplot(pm, aes(x = !! pred, y = !! resp, colour = !! modx,
                      group = !! grp, linetype = !! lty))
  
  
  
  if(insignificant_lines == "dashed"){
    p <- p + geom_path(data = pm, mapping = aes(linewidth = !!lty), show.legend = !facet.modx, show_guide = show_guide_line)
  } else p <- p + geom_path(data = pm, linewidth = line.thickness, show.legend = !facet.modx, show_guide = show_guide_line)
  
  # Plot intervals if requested
  if (interval == TRUE) {
    p <- p + geom_ribbon(data = pm,
                         aes(x = !! pred, ymin = !! sym("ymin"),
                             ymax = !! sym("ymax"), fill = !! modx,
                             group = !! grp, colour = !! modx, linetype = NA),
                         alpha = 1/5, show.legend = FALSE,
                         inherit.aes = TRUE)
  }
  
  # If third mod, facet by third mod
  facet_form <- "~"
  modgroup <- NULL
  # First, decide whether we're faceting at all
  if (!is.null(mod2) || facet.modx == TRUE) {
    do_facets <- TRUE
  } else {do_facets <- FALSE}
  # If faceting by modx, add that to formula
  if (linearity.check == TRUE | facet.modx == TRUE) {
    facet_form <- paste(facet_form, "modx_group")
    modgroup <- "modx_group"
  }
  # If faceting by mod2, add that to formula
  if (!is.null(mod2)) {
    facet_form <- paste(facet_form,
                        ifelse(facet_form == "~", yes = "mod2_group",
                               no = "+ mod2_group"))
    if (!is.null(modgroup)) {
      modgroup <- "modgroup"
    } else {
      modgroup <- "mod2group"
    }
  }
  
  if (do_facets == TRUE) {
    if (!is.null(mod2) & (linearity.check == TRUE | facet.modx == TRUE)) {
      num_unique <- nrow(unique(pm[c("modx_group", "mod2_group")]))
      if (num_unique %in% c(3, 6, 9)) {
        # 1 x 3, 2 x 3, or (most commonly) 3 x 3
        num_cols <- 3
      } else if (num_unique %in% c(4)) {
        # 2 x 2
        num_cols <- 2
      } else { # let ggplot2 decide
        num_cols <- NULL
      }
    } else {num_cols <- NULL}
    p <- p + facet_wrap(as.formula(facet_form), ncol = num_cols)
  }
  
  if (linearity.check == TRUE) {
    p <- p + stat_smooth(data = d,
                         aes(x = !! pred, y = !! resp, group = !! grp),
                         method = "loess", linewidth = 1,
                         show.legend = FALSE, inherit.aes = FALSE,
                         se = FALSE, span = 2, geom = "line",
                         alpha = 0.6, color = "red")
  }
  
  # For factor vars, plotting the observed points
  # and coloring them by factor looks great
  if (plot.points == TRUE) {
    
    if (!is.numeric(d[[as_string(modx)]]) & point.shape) {
      shape_arg <- modx
      # Only show legend if there's a shape aesthetic
      show_legend <- TRUE
    } else {
      shape_arg <- NULL
      show_legend <- FALSE
    }
    constants <- list(alpha = point.alpha)
    if (is.null(weights)) {
      # Only use constant size if weights are not used
      constants$size <- point.size
    }
    # Need to use layer function to programmatically define constant aesthetics
    p <- p + layer(geom = "point", data = d, stat = "identity",
                   inherit.aes = FALSE, show.legend = show_legend,
                   mapping = aes(x = !! pred, y = !! resp, size = !! weights,
                                 group = !! grp, colour = !! modx,
                                 shape = !! shape_arg),
                   position = position_jitter(width = jitter[1],
                                              height = jitter[2]),
                   params = constants) +
      scale_size(range = c(1 * point.size, 5 * point.size), guide = "none")
    
  }
  
  # Rug plot for marginal distributions
  if (rug == TRUE) {
    p <- p + geom_rug(data = d, aes(linetype = NULL),
                      alpha = 0.6,
                      position = position_jitter(width = jitter[1],
                                                 height = jitter[2]),
                      sides = rug.sides, inherit.aes = TRUE)
  }
  
  # # Using theme_apa for theming...but using legend title and side positioning
  # if (is.null(mod2)) {
  #   p <- p + theme_nice(legend.pos = "right")
  # } else {
  #   # make better use of space by putting legend on bottom for facet plots
  #   p <- p + theme_nice(legend.pos = "bottom")
  # }
  
  # Using theme_apa for theming...but using legend title and side positioning
  if (is.null(mod2)) {
    p <- p + theme_light() + theme(legend.pos = "right")
  } else {
    # make better use of space by putting legend on bottom for facet plots
    p <- p + theme_light() + theme(legend.pos = "bottom")
  }
  p <- p + labs(x = x.label, y = y.label) # better labels for axes
  
  # Getting rid of tick marks for factor predictor
  if (length(unique(d[[pred]])) == 2) {
    # Predictor has only two unique values
    # Make sure those values are in increasing order
    brks <- sort(unique(d[[pred]]), decreasing = F)
    if (is.null(pred.labels)) {
      p <- p + scale_x_continuous(breaks = brks)
    } else {
      if (length(pred.labels) == 2) { # Make sure pred.labels has right length
        p <- p + scale_x_continuous(breaks = brks, labels = pred.labels)
      } else {
        warning("pred.labels argument has the wrong length. It won't be used")
        p <- p + scale_x_continuous(breaks = brks)
      }
    }
  }
  
  # Get scale colors, provide better legend title
  if (!is.numeric(d[[as_string(modx)]])) {
    p <- p + scale_colour_manual(name = legend.main, values = colors,
                                 breaks = names(colors),
                                 aesthetics = c("colour", "fill"))
  } else {
    limits <- quant(d[[modx]], probs = c(.1, .9))
    if (min2(modxvals2) < limits[1]) {limits[1] <- min2(modxvals2)}
    if (max2(modxvals2) > limits[2]) {limits[2] <- max2(modxvals2)}
    p <- p + scale_colour_gradientn(name = legend.main,
                                    breaks = modxvals2,
                                    labels = modx.labels,
                                    colors = colors,
                                    limits = limits,
                                    oob = squish,
                                    aesthetics = c("colour", "fill"),
                                    guide = "legend")
  }
  
  if (vary.lty == TRUE) { # Add line-specific changes
    p <- p + scale_linetype_manual(name = legend.main, values = ltypes,
                                   breaks = names(ltypes),
                                   na.value = "blank")
    # Need some extra width to show the linetype pattern fully
    p <- p + theme(legend.key.width = grid::unit(3, "lines"))
  }
  
  if(insignificant_lines == "dashed"){
    p <- p + scale_linetype_manual(name = "Significance of effect", values = c("solid","21"),
                                   breaks = c("significant", "insignificant"),
                                   labels = c("Significant", "Insignificant"),
                                   na.value = "blank") +
      scale_linewidth_manual(name = "Significance of effect", values = c(line.thickness,0.75*line.thickness),
                             breaks = c("significant", "insignificant"),
                             labels = c("Significant", "Insignificant"),
                             na.value = "blank")
  }
  
  # Give the plot the user-specified title if there is one
  if (!is.null(main.title)) {
    p <- p + ggtitle(main.title)
  }
  
  # make sure the legend lines have correct size
  p <- p + guides(color = guide_legend(override.aes = list(linewidth = line.thickness, size = point.size*1.2)),
                  #fill = guide_legend(override.aes = list(linewidth = line.thickness)),
                  linewidth = FALSE) 
  
  return(p)
  
}





#' @importFrom stats as.formula complete.cases df.residual model.frame pt
#' @importFrom stats residuals terms weighted.mean
prep_data <- function(model, d, pred, modx, mod2, pred.values = NULL, saturation_level = NULL, reverse = FALSE, resp = NULL,
                      modx.values, mod2.values, survey, pred.labels = NULL,
                      modx.labels, mod2.labels, wname, weights,
                      linearity.check, robust = FALSE,
                      cluster = NULL,
                      vcov = NULL,
                      interval, set.offset, facvars, centered,
                      preds.per.level, trans_x = NA, 
                      slack_factor = 0.005,
                      x_max = NA,
                      force.cat = FALSE, facet.modx = FALSE,
                      partial.residuals = FALSE, outcome.scale, at, restrict_lines_data = FALSE, insignificant_lines = "normal", verbose = FALSE,  ...) {
  # offset?
  offname <- jtools::get_offset_name(model)
  off <- !is.null(offname)
  
  if (!is.numeric(d[[pred]])) {
    facpred <- TRUE
    if (is.character(d[[pred]])) {d[[pred]] <- factor(d[[pred]])}
  } else if (force.cat == FALSE) {
    facpred <- FALSE
  } else {
    facpred <- TRUE
  }
  
  # Setting default for colors
  if (!is.null(modx) && !is.numeric(d[[modx]])) {
    facmod <- TRUE
    if (is.character(d[[modx]])) {d[[modx]] <- factor(d[[modx]])}
  } else if (force.cat == FALSE | is.null(modx)) {
    facmod <- FALSE
  } else if (!is.null(modx)) {
    facmod <- TRUE
  }
  
  # Treat binary numeric moderators as quasi-categorical
  if (!is.null(modx) && length(unique(d[[modx]])) == 2) {
    if (is.null(modx.values)) {modx.values <- sort(unique(d[[modx]]))}
  }
  
  # Fix character mod2 as well
  if (!is.null(mod2) && !is.numeric(d[[mod2]])) {
    facmod2 <- TRUE
    if (is.character(d[[mod2]])) {d[[mod2]] <- factor(d[[mod2]])}
  } else if (force.cat == FALSE | is.null(mod2)) {
    facmod2 <- FALSE
  } else if (!is.null(mod2)) {
    facmod2 <- TRUE
  }
  
  # Treat binary numeric moderators as quasi-categorical
  if (!is.null(mod2) && length(unique(d[[mod2]])) == 2) {
    if (is.null(mod2.values)) {mod2.values <- sort(unique(d[[mod2]]))}
  }
  
  
  # Get the formula from lm object if given
  formula <- get_formula(model, ...)
  
  # Pulling the name of the response variable for labeling
  if(is.null(resp)) resp <- jtools::get_response_name(model, ...)
  
  # define function to transform y back to original domain (relevant for saturation models)
  if(!is.null(saturation_level)){
    resp_trans <- as.character(deparse(formula[[2]]))
    if(reverse){ 
      inverse_trans <- function(x) exp(x) + saturation_level
    }else{
      inverse_trans <- function(x) saturation_level - exp(x)
    }
  }
  
  # Create a design object
  design <- if (inherits(model, "svyglm")) {
    model$survey.design
  } else {
    NULL
  }
  
  # Warn user if interaction term is absent
  if (!check_interactions(formula, c(pred, modx, mod2))) {
    warn_wrap(paste(c(pred, modx, mod2), collapse = " and "),
              " are not included in an interaction with one another in the
              model.")
  }
  
  ####Getting moderator values ##################################################
  
  if (facpred == TRUE) {
    
    pred.values <- mod_vals(d = d, modx = pred, modx.values = pred.values,
                            survey = survey, weights = weights,
                            design = design,
                            modx.labels = pred.labels, is.mod2 = TRUE,
                            facet.modx = facet.modx, add.varname = FALSE)
    pred.labels <- names(pred.values)
    
  }
  
  if (!is.null(modx)) {
    
    modxvals2 <- mod_vals(d = d, modx = modx, modx.values = modx.values,
                          survey = survey, weights = weights,
                          design = design,
                          modx.labels = modx.labels, any.mod2 = !is.null(mod2),
                          facet.modx = facet.modx, force.cat = force.cat)
    modx.labels <- names(modxvals2)
    
  } else {
    
    modxvals2 <- NULL
    
  }
  
  if (!is.null(mod2)) {
    
    mod2vals2 <- mod_vals(d = d, modx = mod2, modx.values = mod2.values,
                          survey = survey, weights = weights,
                          design = design,
                          modx.labels = mod2.labels, any.mod2 = !is.null(mod2),
                          is.mod2 = TRUE, force.cat = force.cat)
    mod2.labels <- names(mod2vals2)
    
  } else {
    
    mod2vals2 <- NULL
    
  }
  
  ### Drop unwanted factor levels ###############################################
  
  
  if (facpred == TRUE && !is.numeric(d[[pred]])) {
    
    d <- drop_factor_levels(d = d, var = pred, values = pred.values,
                            labels = pred.labels)
    
  }
  
  if (facmod == TRUE && !is.numeric(d[[modx]])) {
    
    d <- drop_factor_levels(d = d, var = modx, values = modxvals2,
                            labels = modx.labels)
    
  }
  
  if (facmod2 == TRUE && !is.numeric(d[[mod2]])) {
    
    d <- drop_factor_levels(d = d, var = mod2, values = mod2vals2,
                            labels = mod2.labels)
    
  }
  
  #### Creating predicted frame #################################################
  

  seq_here <- function(from, to, length.out) {
    if( !is.na(trans_x) & trans_x == "log10") {
      # logarithmic spaced sequence
      return(exp(seq(log(from), log(to), length.out = length.out)))
      } else{
    return(seq(from, to, length.out = length.out))
        }
  }
  

  
  
  if (facpred == TRUE) {
    pred.predicted <- levels(factor(d[[pred]]))
  } else {
    pred.predicted <- seq_here(from = min(d[[pred]], na.rm = TRUE),
                          to = max(d[[pred]], na.rm = TRUE),
                          length.out = preds.per.level)
  }
  
  if (!is.null(modx)) {
    num_combos <- length(modxvals2)
    combos <- expand.grid(modxvals2)
    names(combos) <- modx
  } else {
    num_combos <- 1
  }
  if (!is.null(mod2)) {
    num_combos <- nrow(expand.grid(modxvals2, mod2vals2))
    combos <- expand.grid(modxvals2, mod2vals2)
    names(combos) <- c(modx, mod2)
  }
  
  pms <- list()
  
  for (i in seq_len(num_combos)) {
    
    at_list <- list()
    if (!is.null(modx)) {
      at_list[[modx]] <- combos[i, modx]
      if(restrict_lines_data) {
        slack_diff <-  slack_factor*(min(c(max(d[[pred]], na.rm = TRUE), x_max), na.rm = TRUE) - min(d[[pred]], na.rm = TRUE))
        slack_ratio <- (min(c(max(d[[pred]], na.rm = TRUE), x_max), na.rm = TRUE)/ min(d[[pred]], na.rm = TRUE))^slack_factor
        # print("slack_ratio lines")
        # print(slack_ratio)
        x_min_category <- min(d[d[[modx]] == combos[i, modx],][[pred]], na.rm = TRUE)
        x_max_category <- max(d[d[[modx]] == combos[i, modx],][[pred]], na.rm = TRUE)
        # print("min max of category")
        # print(c(x_min_category, x_max_category))
        
        x_lower <- ifelse(!is.na(trans_x) & trans_x == "log10",
                        #max(min(d[[pred]], na.rm = TRUE), 
                        x_min_category / slack_ratio,
                        #max(min(d[[pred]], na.rm = TRUE), 
                        x_min_category - slack_diff)
        x_upper <- ifelse(!is.na(trans_x) & trans_x == "log10",
                        # min(max(d[[pred]], na.rm = TRUE),  
                            x_max_category * slack_ratio,
                        # min(max(d[[pred]], na.rm = TRUE), 
                        x_max_category + slack_diff)
        pred.predicted <- seq_here(from = x_lower,
                              to = x_upper,
                              length.out = preds.per.level)
        # print(c(x_lower, x_upper))
      }
    }
    if (!is.null(mod2)) {
      at_list[[mod2]] <- combos[i, mod2]
    }
    
    if (!is.null(at)) {
      at_list[names(at)] <- at
    }
    
    suppressMessages({pms[[i]] <- jtools::make_predictions(
      model = model, data = d, pred = pred, pred.values = pred.predicted,
      at = at_list, set.offset = set.offset, center = centered,
      interval = interval, outcome.scale = outcome.scale, robust = robust,
      cluster = cluster,
      vcov = vcov, ...
    ) 
    if(!is.null(saturation_level)){
      pms[[i]] <-  pms[[i]] %>% mutate(.,"{resp}" := inverse_trans(!!sym(resp_trans)),
                                       ymin = inverse_trans(ymin),
                                       ymax = inverse_trans(ymax))}})
    
    # only looking for completeness in these variables
    check_vars <- all.vars(get_formula(model, ...)) %just% names(pms[[i]])
    pms[[i]] <-
      pms[[i]][complete.cases(pms[[i]][check_vars]), ]
  }
  
  if (off == TRUE) {
    msg_wrap("Outcome is based on a total of ", set.offset, " exposures.")
  }
  
  pm <- do.call("rbind", pms)
  
  
  if(verbose){
    print(head(pm))
  } 
  
  # kick out predictions for those groups where the effect of the predictor is insignificant
  if(insignificant_lines != "normal"){
    insignificant_groups <- get_insignificant_groups(obj = model, pred = as_name(pred), categorical_moderator = modx, levels = modx.values, vcov = vcov, alpha = 0.05)
    if(insignificant_lines == "remove") pm <- pm %>% filter(!(!!sym(modx) %in% insignificant_groups))
    if(insignificant_lines == "dashed") pm <- pm %>% mutate(significance_line = as.factor(if_else(!!sym(modx) %in% insignificant_groups, "insignificant", "significant")))
  }
  
  if(verbose){
    print(pm)
  } 
  
  # Do partial residuals if requested
  if (partial.residuals == TRUE) {
    suppressMessages({
      d <- partialize(model, vars = c(pred, modx, mod2), center = centered,
                      data = d, scale = outcome.scale, set.offset = set.offset,
                      ...)
    })
  }
  
  ## Prep original data for splitting into groups ##
  if (!is.null(modx)) {
    d <- split_int_data(d = d, modx = modx, mod2 = mod2,
                        linearity.check = linearity.check | facet.modx,
                        modx.values = modx.values,
                        modxvals2 = modxvals2, mod2.values = mod2.values,
                        mod2vals2 = mod2vals2, facmod = facmod,
                        facmod2 = facmod2)
  }
  
  # Labels for values of moderator
  if (!is.null(modx) && !is.numeric(d[[modx]])) {
    pm[[modx]] <- factor(pm[[modx]], levels = modxvals2, labels = modx.labels)
  }
  if (facmod == TRUE) {
    d[[modx]] <- factor(d[[modx]], levels = modxvals2, labels = modx.labels)
  }
  if (!is.null(modx)) {
    if (is.numeric(d[[modx]])) {
      pm$modx_group <- factor(pm[[modx]], levels = modxvals2,
                              labels = modx.labels)
    } else {
      pm$modx_group <- factor(pm[[modx]], levels = modx.labels)
    }
  }
  
  # Setting labels for second moderator
  if (!is.null(mod2)) {
    
    # Convert character moderators to factor
    if (!is.numeric(d[[mod2]])) {
      d[[mod2]] <- factor(d[[mod2]], levels = mod2vals2, labels = mod2.labels)
      pm[[mod2]] <- factor(pm[[mod2]], levels = mod2vals2, labels = mod2.labels)
      pm$mod2_group <- pm[[mod2]]
    } else {
      pm$mod2_group <- factor(pm[[mod2]], levels = mod2vals2,
                              labels = mod2.labels)
    }
    
  }
  
  
  # Dealing with transformations of the dependent variable
  # Have to make sure not to confuse situations with brmsfit objects and
  # distributional DVs
  if (resp %nin% names(d) & "dpar" %nin% names(list(...))) {
    trans_name <- as.character(deparse(formula[[2]]))
    d[[trans_name]] <- eval(formula[[2]], d)
  }
  
  
  
  out <- list(predicted = pm, original = d)
  out <- structure(out, resp = resp, facmod = facmod,
                   pred.values = pred.values, pred.labels = pred.labels,
                   modxvals2 = modxvals2, modx.labels = modx.labels,
                   mod2vals2 = mod2vals2, mod2.labels = mod2.labels,
                   facet.modx = facet.modx)
  return(out)
  
}

split_int_data <- function(d, modx, mod2, linearity.check, modx.values,
                           modxvals2, mod2.values, mod2vals2, facmod, facmod2) {
  # For numeric, non-binary moderators...
  if (facmod == FALSE &
      !(length(unique(d[[modx]])) == 2 & length(modxvals2) == 2)) {
    
    # Use ecdf function to get quantile of the modxvals
    mod_val_qs <- ecdf(d[[modx]])(sort(modxvals2))
    
    # Now I am going to split the data in a way that roughly puts each modxval
    # in the middle of each group. mod_val_qs is a vector of quantiles for each
    # modxval, so I will now build a vector of the midpoint between each
    # neighboring pair of quantiles — they will become the cutpoints for
    # splitting the data into groups that roughly correspond to the modxvals
    cut_points <- c() # empty vector
    # Iterate to allow this to work regardless of number of modxvals
    for (i in 1:(length(modxvals2) - 1)) {
      
      cut_points <- c(cut_points, mean(mod_val_qs[i:(i + 1)]))
      
    }
    
    # Add Inf to both ends to encompass all values outside the cut points
    cut_points <- c(-Inf, quantile(d[[modx]], cut_points, na.rm = TRUE), Inf)
    
    # Create variable storing this info as a factor
    d["modx_group"] <- cut(d[[modx]], cut_points,
                           labels = names(sort(modxvals2)))
    
    if (!is.null(modx.values) && modx.values[1] == "terciles") {
      d$modx_group <- factor(cut2(d[[modx]], g = 3, levels.mean = TRUE),
                             labels = c(paste("Lower tercile of", modx),
                                        paste("Middle tercile of", modx),
                                        paste("Upper tercile of", modx)))
    }
    
  } else {
    
    d["modx_group"] <- factor(d[[modx]], levels = modxvals2,
                              labels = names(modxvals2))
    
  }
  
  if (!is.null(mod2)) {
    if (facmod2 == FALSE &
        !(length(unique(d[[mod2]])) == 2 & length(mod2vals2) == 2)) {
      
      mod_val_qs <- ecdf(d[[mod2]])(sort(mod2vals2))
      
      
      cut_points2 <- c()
      for (i in 1:(length(mod2vals2) - 1)) {
        
        cut_points2 <- c(cut_points2, mean(mod_val_qs[i:(i + 1)]))
        
      }
      
      cut_points2 <- c(-Inf, quantile(d[[mod2]], cut_points2, na.rm = TRUE),
                       Inf)
      
      d["mod2_group"] <- cut(d[[mod2]], cut_points2,
                             labels = names(sort(mod2vals2)))
      
      if (!is.null(mod2.values) && mod2.values[1] == "terciles") {
        d$mod2_group <- factor(cut2(d[[mod2]], g = 3, levels.mean = TRUE),
                               labels = c(paste("Lower tercile of", mod2),
                                          paste("Middle tercile of", mod2),
                                          paste("Upper tercile of", mod2)))
      }
      
    } else if (facmod2 == TRUE) {
      
      d["mod2_group"] <- factor(d[[mod2]], levels = mod2vals2,
                                labels = names(mod2vals2))
      
    }
  }
  
  return(d)
  
}

drop_factor_levels <- function(d, var, values, labels) {
  # Need to save the rownames because of tibble's stupidity
  the_row_names <- rownames(d)
  the_row_names <- the_row_names[d[[var]] %in% values]
  d <- d[d[[var]] %in% values,]
  d[[var]] <- factor(d[[var]], levels = values)
  # Can't use rowname assignment method because of tibble's stupidity
  attr(d, "row.names") <- the_row_names
  return(d)
  
}


any_interaction <- function(formula) {
  any(attr(terms(formula), "order") > 1)
}

get_interactions <- function(formula) {
  if (any_interaction(formula)) {
    ts <- terms(formula)
    labs <- paste("~", attr(ts, "term.labels"))
    forms <- lapply(labs, as.formula)
    forms <- forms[which(attr(ts, "order") > 1)]
    ints <- lapply(forms, all.vars)
    names(ints) <- attr(ts, "term.labels")[which(attr(ts, "order") > 1)]
    return(ints)
  } else {
    NULL
  }
}

check_interactions <- function(formula, vars) {
  vars <- vars[!is.null(vars)]
  if (any_interaction(formula)) {
    checks <- sapply(get_interactions(formula), function(x, vars) {
      if (all(vars %in% x)) TRUE else FALSE
    }, vars = vars)
    any(checks)
  } else {
    FALSE
  }
}


#### Non-exported functions ###################################################

## Utility function for getting values of moderator values for interaction
## functions

mod_vals <- function(d, modx, modx.values, survey, weights,
                     design = design, modx.labels = NULL,
                     any.mod2 = FALSE, is.mod2 = FALSE,
                     sims = FALSE, facet.modx = FALSE, force.cat = FALSE,
                     add.varname = TRUE) {
  
  # Get moderator mean
  if (survey == FALSE & is.numeric(d[[modx]])) {
    weights <- if (is.null(weights)) {
      rep(1, nrow(d))
    } else if (is.character(weights)) {
      d[[weights]]
    } else weights
    modmean <- weighted.mean(d[[modx]], weights, na.rm = TRUE)
    modsd <- wtd.sd(d[[modx]], weights)
    
  } else if (survey == TRUE & is.numeric(d[[modx]])) {
    
    modsd <- svysd(as.formula(paste("~", modx, sep = "")), design = design)
    # Have to construct the formula this way since the syntax for svymean
    # differs from mean
    modmean <- survey::svymean(as.formula(paste("~", modx, sep = "")),
                               design = design)
    
  }
  
  is_fac <- if (!is.numeric(d[[modx]]) | force.cat == TRUE) TRUE else FALSE
  
  # Testing whether modx.values refers to pre-defined arg or list of factor levels
  predefined_args <- c("mean-plus-minus", "plus-minus", "terciles")
  if (is.character(modx.values) & length(modx.values) == 1) {
    char1 <- if (modx.values %in% predefined_args) TRUE else FALSE
    if (is_fac == TRUE & char1 == TRUE) {
      stop_wrap(modx.values, " is not supported for a non-numeric moderator.")
    } else if (is_fac == FALSE & char1 == FALSE) {
      stop_wrap(modx.values, " is not a valid ",
                ifelse(is.mod2, yes = "mod2.values", no = "modx.values"),
                " argument for a numeric moderator.")
    }
  } else {char1 <- FALSE}
  
  user_specified <- length(modx.values) > 1
  
  # If using a preset, send to auto_mod_vals function
  if (is_fac == FALSE && (is.null(modx.values) | is.character(modx.values))) {
    
    modxvals2 <- auto_mod_vals(d, modx.values = modx.values, modx = modx,
                               modmean = modmean, modsd = modsd,
                               modx.labels = modx.labels,
                               mod2 = (is.mod2 | facet.modx),
                               sims = sims, add.varname = add.varname)
    
  }
  
  # For user-specified numbers or factors, go here
  if (is.null(modx.values) & is_fac == TRUE) {
    
    modxvals2 <- ulevels(d[[modx]])
    if (is.null(modx.labels)) {
      
      if ((is.mod2 | facet.modx) & add.varname == TRUE) {
        modx.labels <- paste(modx, "=", ulevels(d[[modx]]))
      } else {
        modx.labels <- ulevels(d[[modx]])
      }
      
    }
    names(modxvals2) <- modx.labels
    
  } else if (!is.null(modx.values) & char1 == FALSE) {
    # Use user-supplied values otherwise
    
    if (!is.null(modx.labels)) {
      # What I'm doing here is preserving the label order
      names(modx.values) <- modx.labels
      if (!is.mod2 & !is.factor(d[[modx]])) {
        modxvals2 <- rev(modx.values)
      } else {
        modxvals2 <- modx.values
      }
      modx.labels <- names(modxvals2)
      
    } else {
      
      names(modx.values) <- if ((is.mod2 | facet.modx) & add.varname == TRUE) {
        paste(modx, "=", modx.values)
      } else {
        modx.values
      }
      if (!is.mod2 & !is.factor(d[[modx]])) {
        modxvals2 <- rev(modx.values)
      } else {
        modxvals2 <- modx.values
      }
      modx.labels <- names(modxvals2)
      
    }
    
  }
  
  if (is.null(modx.labels)) {
    # Name the modx.labels object with modxvals2 names
    
    modx.labels <- if ((is.mod2 | facet.modx) & add.varname == TRUE) {
      paste(modx, "=", modxvals2)
    } else {
      names(modxvals2)
    }
    
  }
  
  # Hacky way to have a shorthand to drop NA
  range2 <- function(...) {
    range(..., na.rm = TRUE)
  }
  if (is_fac == FALSE & user_specified == FALSE) {
    # The proper order for interact_plot depends on presence of second moderator
    modxvals2 <- sort(modxvals2, decreasing = (!any.mod2 & !facet.modx))
    if (any(modxvals2 > range2(d[[modx]])[2])) {
      warn_wrap(paste(modxvals2[which(modxvals2 > range2(d[[modx]])[2])],
                      collapse = " and "), " is outside the observed range of ",
                modx)
    }
    if (any(modxvals2 < range2(d[[modx]])[1])) {
      warn_wrap(paste(modxvals2[which(modxvals2 < range2(d[[modx]])[1])],
                      collapse = " and "), " is outside the observed range of ",
                modx)
    }
  }
  
  return(modxvals2)
  
}

## Gets the preset values, e.g., mean plus/minus 1 SD

auto_mod_vals <-
  function(d, modx, modx.values, modmean, modsd, modx.labels = NULL,
           mod2 = FALSE, sims = FALSE, add.varname = TRUE) {
    
    # Default to +/- 1 SD unless modx is factor
    if ((is.null(modx.values) || modx.values == "mean-plus-minus") &
        length(unique(d[[modx]])) > 2) {
      
      modxvals2 <- c(modmean - modsd,
                     modmean,
                     modmean + modsd)
      if (mod2 == FALSE) {
        names(modxvals2) <- c("- 1 SD", "Mean", "+ 1 SD")
      } else {
        names(modxvals2) <- c(paste("Mean of", modx, "- 1 SD"),
                              paste("Mean of", modx),
                              paste("Mean of", modx, "+ 1 SD"))
      }
      
    } else if (!is.null(modx.values) && modx.values[1] == "plus-minus") { # No mean
      
      modxvals2 <- c(modmean - modsd, modmean + modsd)
      if (mod2 == FALSE) {
        names(modxvals2) <- c("- 1 SD", "+ 1 SD")
      } else {
        names(modxvals2) <- c(paste("Mean of", modx, "- 1 SD"),
                              paste("Mean of", modx, "+ 1 SD"))
      }
      
    } else if (!is.null(modx.values) && modx.values[1] == "terciles") {
      
      x_or_2 <- switch(as.character(mod2),
                       "TRUE" = "2",
                       "FALSE" = "x")
      group_name <- paste0("mod", x_or_2)
      d[[group_name]] <- cut2(d[[modx]], g = 3, levels.median = TRUE)
      modxvals2 <- as.numeric(levels(d[[group_name]]))
      msg_wrap("Medians of each tercile of ", modx, " are ",
               paste(modxvals2, collapse = ", "))
      
      if (mod2 == FALSE) {
        names(modxvals2) <- c("Lower tercile median", "Middle tercile median",
                              "Upper tercile median")
      } else {
        names(modxvals2) <- c(paste("Lower tercile of", modx),
                              paste("Middle tercile of", modx),
                              paste("Upper tercile of", modx))
      }
      
    } else if (is.null(modx.values) & length(unique(d[[modx]])) == 2) {
      
      modxvals2 <- as.numeric(levels(factor(d[[modx]])))
      if (!is.null(modx.labels)) {
        
        names(modxvals2) <- modx.labels
        
      } else {
        
        if (mod2 == TRUE & sims == FALSE & add.varname == TRUE) {
          names(modxvals2) <-
            sapply(modxvals2, FUN = function(x) {paste(modx, "=", round(x,3))})
        } else {
          names(modxvals2) <- modxvals2
        }
        
      }
      
    }
    
    if (!is.null(modx.labels)) {
      if (length(modx.labels) == length(modxvals2)) {
        names(modxvals2) <- modx.labels
      } else {
        warn_wrap("modx.labels or mod2.labels argument is not the same length
                  as the number of moderator values used. It will be ignored.")
      }
    }
    
    return(modxvals2)
    
  }


## Centering

center_ss <- function(d, weights, facvars = NULL, fvars, pred, resp, modx,
                      survey, design = NULL, mod2, wname, offname, centered) {
  
  # Just need to pick a helper function based on survey vs no survey
  if (survey == TRUE) {
    
    out <- center_ss_survey(d, weights, facvars, fvars, pred, resp, modx,
                            survey, design, mod2, wname, offname, centered)
    
  } else {
    
    out <- center_ss_non_survey(d, weights, facvars, fvars, pred, resp, modx,
                                mod2, wname, offname, centered)
    
  }
  
  return(out)
  
  
}

## If not svydesign, centering is fairly straightforward

center_ss_non_survey <- function(d, weights, facvars = NULL, fvars, pred,
                                 resp, modx, mod2, wname, offname, centered) {
  
  omitvars <- c(pred, resp, modx, mod2, wname, offname)
  
  # Dealing with two-level factors that aren't part of an interaction
  # /focal pred
  fv2 <- fvars[fvars %nin% omitvars]
  
  # Handling user-requested centered vars
  if (centered[1] != "all" && centered[1] != "none") {
    
    if (any(omitvars %in% centered)) {
      warning("Moderators, outcome variables, and weights/offsets",
              " cannot be centered.")
      centered <- centered[centered %nin% omitvars]
    }
    if (length(centered) > 0) {
      d <- gscale(vars = centered, data = d, center.only = TRUE,
                  weights = weights)
    }
    
    for (v in fv2) {
      
      if (is.factor(d[[v]]) &&
          length(unique(d[[v]])) == 2 && v %nin% centered) {
        
        facvars <- c(facvars, v)
        
      }
      
    }
    
  } else if (centered[1] == "all") {
    
    # saving all vars except response
    vars <- names(d)[names(d) %nin% omitvars]
    
    d <- gscale(vars = vars, data = d, center.only = TRUE,
                weights = weights)
    
  } else if (centered == "none") {
    
    # Dealing with two-level factors that aren't part
    # of an interaction/focal pred
    for (v in fv2) {
      if (is.factor(d[[v]]) & length(unique(d[[v]])) == 2) {
        
        facvars <- c(facvars, v)
        
      }
    }
    
  }
  
  # Fixes a data type error with predict() later
  d <- as.data.frame(d)
  
  out <- list(d = d, facvars = facvars, fvars = fvars, design = NULL)
  
  return(out)
  
}

## Svydesigns get their own function to make control flow easier to follow

center_ss_survey <- function(d, weights, facvars = NULL, fvars, pred, resp,
                             modx, survey, design, mod2, wname, offname,
                             centered) {
  
  omitvars <- c(pred, resp, modx, mod2, wname, offname)
  
  # Dealing with two-level factors that aren't part of an interaction
  # /focal pred
  fv2 <- fvars[fvars %nin% omitvars]
  
  # Handling user-requested centered vars
  if (centered[1] != "all" && centered[1] != "none") {
    
    if (any(omitvars %in% centered)) {
      warning("Moderators, outcome variables, and weights/offsets",
              " cannot be centered.")
      centered <- centered[centered %nin% omitvars]
    }
    design <- gscale(vars = centered, data = design, center.only = TRUE)
    d <- design$variables
    
    # Dealing with two-level factors that aren't part
    # of an interaction/focal pred
    for (v in fv2) {
      if (is.factor(d[[v]]) &&
          length(unique(d[[v]])) == 2 && v %nin% centered) {
        
        facvars <- c(facvars, v)
        
      }
    }
    
  } else if (centered == "none") {
    
    # Dealing with two-level factors that aren't part
    # of an interaction/focal pred
    for (v in fv2) {
      if (is.factor(d[[v]]) && length(unique(d[[v]])) == 2) {
        
        facvars <- c(facvars, v)
        
      }
    }
    
  } else if (centered == "all") {
    
    # Center all non-focal
    ndfvars <- fvars[fvars %nin% omitvars]
    
    if (length(ndfvars) > 0) {
      design <- gscale(vars = ndfvars, data = design, center.only = TRUE)
      d <- design$variables
    }
    
  }
  
  out <- list(d = d, design = design, facvars = facvars, fvars = fvars)
  
  return(out)
  
}
