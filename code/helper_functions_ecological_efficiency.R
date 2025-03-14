
# get average values for predictor variables
cross_section <- function(df, start_average, end_average, year_outcome, input_variable = "GDP_PPP", outcome_variable, aggregated_control = NULL){
  df %>% 
    dplyr::select(-!!sym(outcome_variable)) %>% 
    filter(year>=start_average & year <= end_average) %>% 
    group_by(country_code, country_name, regime_type, former_state) %>% 
    { if(!is.null(aggregated_control)) dplyr::summarise(., across(c(!!sym(input_variable), !!sym(aggregated_control)), ~ mean(.,na.rm = TRUE))) else dplyr::summarise(.,across(c(!!sym(input_variable)), ~ mean(.,na.rm = TRUE)))} %>% 
    full_join(df %>%
                filter(year == year_outcome) %>% 
                select(country_code, outcome_variable)) %>% 
    ungroup()
}


# reduce data set to former states
former_states <- function(df, aggregated_control= NULL){
  df %>% 
    { if(!is.null(aggregated_control)) mutate(., !!aggregated_control := case_match(
    country_name, 
    "USSR" ~ df %>% filter(former_state == "USSR") %>% pull(!!sym(aggregated_control)) %>% magrittr::extract(1),
    "Czechoslovakia" ~ df %>% filter(former_state == "CSK") %>% pull(!!sym(aggregated_control)) %>% magrittr::extract(1),
    "Yugoslavia" ~ df %>% filter(former_state == "Yugo") %>% pull(!!sym(aggregated_control)) %>% magrittr::extract(1),
    .default = !!sym(aggregated_control))) else .} %>% 
    filter(is.na(former_state)) %>% 
    select(-former_state)
}

# aggregate different regime types into categories and choose regime type pf Vietnam (socialist or revolutionary)
update_regime_types <- function(df, granularity, reference_regime, regime_vietnam){
  if(granularity == 4){
    df %>% mutate(regime_type = case_match(
      regime_type, c("Corporatist core", "Liberal core", "Social democratic core") ~ "capitalist_core", 
      "Capitalist periphery" ~ "capitalist_periphery", 
      "Revolutionary states" ~ "revolutionary", 
      "Socialist states" ~ "socialist")) %>%
      mutate(regime_type = case_match(country_name, "Vietnam" ~ regime_vietnam, .default = regime_type)) %>% 
      mutate(regime_type = relevel(factor(regime_type), ref = reference_regime))
  }
}


# customised border handling functions for ggplot
oob_squish_upper <- function(x, range = c(0, 1), only.finite = TRUE) {
  force(range)
  finite <- if (only.finite) is.finite(x) else TRUE
  x[finite & x > range[2]] <- range[2]
  x
}

oob_squish_lower <- function(x, range = c(0, 1), only.finite = TRUE) {
  force(range)
  finite <- if (only.finite) is.finite(x) else TRUE
  x[finite & x > range[1]] <- range[1]
  x
}




interact_plot_custom <- function(model, vcov, data, pred = "GDP_PPP", resp = NULL, saturation_level = NULL, x.label  = NULL, 
                                 y.label = NULL, physical_boundary = NULL, reverse = FALSE, x_max = NA, y_max = NA, y_min = NA, percentage = FALSE,  
                                 ignore_significance_right = FALSE, verbose = FALSE, insignificant_lines = "normal",
                                 point.size = 1, line.thickness = 1, trans_x = "identity", slack_factor = 0.005, natural_extremum = TRUE){
  pred <- sym(pred)
  if(is.null(resp)) resp <- jtools::get_response_name(model)
  
  ticks_format <- if(percentage) scales::percent_format(scale = 1) else scales::comma_format()
  
  p <- interact_plot(model, 
                     vcov = vcov,
                     pred = !!pred, 
                     resp = resp,
                     saturation_level = saturation_level,
                     reverse = reverse,
                     modx = regime_type, 
                     modx.values = c(#"revolutionary", 
                       "capitalist_periphery",
                       "capitalist_core",  
                       "socialist"), 
                     modx.labels = c(#"revolutionary",
                       "Capitalist periphery",
                       "Capitalist core",  
                       "Socialist"),
                     legend.main = "Regime type",
                     colors = c(#"lightgrey",
                       "#004488", "#DDAA33", "#BB5566"),
                     vary.lty = FALSE,
                     data = data, 
                     centered = "all",
                     x.label = x.label, 
                     y.label = y.label, 
                     robust = TRUE, 
                     plot.points = TRUE, 
                     point.shape = TRUE, 
                     point.size = point.size, 
                     point.alpha = 0.7,
                     line.thickness = line.thickness,
                     interval = TRUE, 
                     rug=FALSE,
                     restrict_lines_data = TRUE,
                     insignificant_lines = insignificant_lines,
                     show_guide_line = FALSE,
                     trans_x = trans_x,
                     verbose = verbose,
                     x_max = x_max,
                     slack_factor = slack_factor) +
    geom_point(mapping = aes(x=!!pred, y=!!sym(resp), colour = regime_type, shape = regime_type),
               show.legend = FALSE,
               inherit.aes = FALSE,
               data = data %>% filter(regime_type == "revolutionary") %>% mutate(regime_type = "Revolutionary"),
               #color = "darkgrey",
               #shape = 17,
               size = 1.5 * point.size,
               alpha = 0.7
    ) +
    scale_discrete_manual(aesthetics = c("colour", "fill"),
                          breaks = c("Capitalist core", "Capitalist periphery", "Socialist", "Revolutionary"),
                          values = c("Capitalist core" = "#004488", "Capitalist periphery" = "#DDAA33", "Socialist" = "#BB5566", "Revolutionary" = "#777777"),
                          name = "Regime type") +
    scale_shape_manual(breaks = c("Capitalist core", "Capitalist periphery", "Socialist", "Revolutionary"),
                       values = c("Capitalist core" = 19, "Capitalist periphery" = 19, "Socialist" = 19, "Revolutionary" = 17),
                       name = "Regime type") +
    scale_y_continuous(limits = c(y_min, y_max), oob= ifelse(natural_extremum, oob_squish, oob_squish_upper), labels = ticks_format) +
    #coord_cartesian(xlim=c(0,NA), ylim=c(0,NA), expand = FALSE) +
    {if(!is.null(physical_boundary)) geom_vline(aes(xintercept=physical_boundary), color = "#555555")} +
    theme(text = element_text(size=9, family="Computer Modern Sans"), legend.position="top") +
    guides(linetype = FALSE)
  
  
  slack_diff <-  slack_factor*(min(c(max(data[[pred]], na.rm = TRUE), x_max), na.rm = TRUE) - 
                  min(data[[pred]], na.rm = TRUE))
  
  slack_ratio <-  (min(c(max(data[[pred]], na.rm = TRUE), x_max), na.rm = TRUE) / 
                                 min(data[[pred]], na.rm = TRUE))^slack_factor
  
  # print("slack_ratio grey boxes")
  # print(slack_ratio)
  
  max_socialist <- max(data[data[["regime_type"]] == "socialist",][[pred]], na.rm = TRUE) 
  min_socialist <- min(data[data[["regime_type"]] == "socialist",][[pred]], na.rm = TRUE) 
  #this makes sure that the grey boxes 
  #only show up in the range the socialist countries actually occupy, with a little bit of slack
  
  x_min_grey <- if_else(!is.na(trans_x) & trans_x == "log10",
                        max(c(layer_scales(p)$x$range$range[1], min_socialist / slack_ratio), na.rm = TRUE),
                        max(c(layer_scales(p)$x$range$range[1], min_socialist - slack_diff), na.rm = TRUE))
  
  x_max_grey <- 
    min(c(
      x_max,
      layer_scales(p)$x$range$range[2],
      ifelse(!is.na(trans_x) & trans_x == "log10",
             max_socialist * slack_ratio,
             max_socialist + slack_diff)), na.rm = TRUE)
  
  # print("minimal and maximal value:")
  # print(c(x_min_grey, x_max_grey))
                 
  
  insignificant_ranges <-  significant_ranges(obj = model, 
                                              vcov = vcov,
                                              varnames = c("regime_typesocialist", 
                                                           paste0("log(", as_name(pred), ")")),
                                               find_insignificant = TRUE) %>% lapply(exp) %>%
    lapply(function(range){
      lower <- max(range[1], x_min_grey)
      upper <- min(range[2], x_max_grey) #ifelse(lower < x_max_grey, min(range[2], x_max_grey), range[2])
      return(c(lower, upper))
    })
  
  # print("insignificant ranges raw:")
  # print(significant_ranges(obj = model, 
  #                          vcov = vcov,
  #                          varnames = c("regime_typesocialist", 
  #                                       paste0("log(", as_name(pred), ")")),
  #                          find_insignificant = TRUE) %>% lapply(exp))
  # 
  # print("insignificant ranges:")
  # print(insignificant_ranges)
  # print("result of sapply:")
  # print(sapply(insignificant_ranges, function(range) range[1]< x_max_grey))
  # print("type of result:")
  # print(typeof(sapply(insignificant_ranges, function(range) range[1]< x_max_grey)))
  
  if(length(insignificant_ranges)>0){
    insignificant_box_layers <- 
      lapply(insignificant_ranges[sapply(insignificant_ranges, function(range) range[1]< x_max_grey & range[2] > x_min_grey)], 
             function(range){
                # print("range")
                # print(range)
                geom_rect(aes(xmin = range[1], 
                              xmax = range[2], 
                              ymin = -Inf, ymax = Inf), 
                          inherit.aes = FALSE, fill = "lightgrey", alpha = 0.014)
    })
    
  
    p$layers <- append(p$layers, insignificant_box_layers, after=0)
  }
  

  
  p +
    scale_x_continuous(limits = c(ifelse(trans_x == "identity", 0, NA) ,x_max), 
                       labels = scales::comma_format(), 
                       trans = trans_x) + #, expand = c(0,0.05)) +
    {if(layer_scales(p)$y$range$range[1] < 0) scale_y_continuous(limits=c(0,y_max),oob= ifelse(natural_extremum, oob_squish, oob_squish_upper), labels = ticks_format)} +
    {if(reverse) list(scale_y_reverse(limits=c(y_max,y_min),oob= ifelse(natural_extremum, oob_squish, oob_squish_upper), labels = ticks_format))}
}




countries_to_tex <- function(data, file_path){
  
  library(glue)
  
  # Define the order of regime types and their corresponding LaTeX labels
  regime_order <- c("socialist", "capitalist_periphery", "capitalist_core", "revolutionary")
  regime_labels <- c(
    "socialist" = "Socialist:",
    "capitalist_periphery" = "Capitalist periphery:",
    "capitalist_core" = "Capitalist core:",
    "revolutionary" = "Revolutionary:"
  )
  
  # Create a data frame to ensure all regime types are included
  all_regimes <- tibble(regime_type = regime_order)
  
  # Generate the LaTeX fragment, ensuring all regime types are present
  latex_output <- data %>%
    filter(regime_type %in% regime_order) %>%  # Only keep the regime types in the order list
    mutate(regime_type = factor(regime_type, levels = regime_order)) %>%  # Ensure the order
    group_by(regime_type) %>%
    summarize(countries = paste(sort(country_name), collapse = ", ")) %>%
    right_join(all_regimes, by = "regime_type") %>%  # Ensure all regime types are included
    mutate(countries = ifelse(is.na(countries), "/", countries)) %>%  # Set "/" for missing regime types
    mutate(regime_type = factor(regime_type, levels = regime_order)) %>%  # Re-apply the correct order
    arrange(regime_type) %>%  # Sort by the factor levels to preserve the desired order
    mutate(latex_entry = glue("\\item[{regime_labels[regime_type]}] {countries}")) %>%
    pull(latex_entry)
  
  
  # Combine into the full LaTeX fragment
  latex_fragment <- glue("
  \\begin{{description}}
  {paste(latex_output, collapse = '\\n')}
  \\end{{description}}
  ")
  
  # Save the result
  writeLines(latex_fragment, con = file_path)
  
}


generate_supplementary_information <- function(model, vcov, model_sat = NULL, vcov_sat = NULL, data, pred, 
                                               resp, y.label, percentage, reverse, sat, y_max = NA, y_min = NA,
                                               natural_extremum = TRUE){
  
  # create subfolder
  subfolder <- paste0("../figures/supplementary_information/", pred, "_", resp)
  if (!dir.exists(subfolder)) {
    dir.create(subfolder)
  }
  
  # determine input labels
  x.name <- switch(pred,
                   "GDP_PPP" = "GDP per capita",
                   "material_footprint" = "Material footprint per capita",
                   "co2" = "CO$_2$ emissions per capita")
  x.label <- switch(pred,
                    "GDP_PPP" = "GDP ($ per capita)",
                    "material_footprint" = "Material footprint (tons per capita)",
                    "co2" = expression(paste(CO[2], " emissions (tons per capita)"))
  )
  
  
  
  # Logarithm plot
  log_plot <-
    interact_plot_custom(model, vcov, data, pred = pred, x.label = x.label, x_max = NA,
                         point.size = 0.8, line.thickness = 0.8, insignificant_lines = "dashed", y.label = y.label,
                         reverse = reverse, percentage = percentage, resp = resp,
                         trans_x = "log10",
                         y_max = y_max, y_min = y_min,
                         natural_extremum = natural_extremum)
  
  ggsave(file = paste0("../figures/supplementary_information/", pred, "_", resp, "/big_plot.pdf"), log_plot, width = 160, height = 115, units = "mm", device=cairo_pdf)
  
  
  
  # Saturation plot
  if(!is.null(model_sat)){
    sat_plot <-
      interact_plot_custom(model_sat, vcov_sat, data, pred = pred, x.label = x.label, x_max = NA,
                           point.size = 0.8, line.thickness = 0.8, insignificant_lines = "dashed", y.label = y.label,
                           reverse = reverse, percentage = percentage, resp = resp, y_max = y_max, y_min = y_min, 
                           trans_x = "log10", saturation_level = sat,
                           natural_extremum = natural_extremum)
    
    ggsave(file = paste0("../figures/supplementary_information/", pred, "_", resp, "/big_plot_sat.pdf"),
           sat_plot, width = 160, height = 115, units = "mm", device=cairo_pdf)
  }
  
  
  
  
  # stargazer table
  if(!is.null(model_sat)){
    models <- list(model,model_sat)
    standard_errors <- list(sqrt(diag(vcov)), sqrt(diag(vcov_sat)))
  } else {
    models <- list(model)
    standard_errors <- list(sqrt(diag(vcov)))
  }
  
  stargazer(models,
            #type = "text",
            out = paste0("../figures/supplementary_information/", pred, "_", resp, "/stargazer.tex"),
            report = 'vc*',
            se = standard_errors,
            intercept.bottom = FALSE,
            #column.labels = c("Model 1", "Model 2", "Model 3"),
            model.numbers = FALSE,
            align = TRUE,
            dep.var.caption = y.label,
            dep.var.labels.include = TRUE,
            no.space = TRUE,
            title = paste0("Model coefficients"),
            dep.var.labels = c("Logarithmic model", "Saturation model"),
            covariate.labels =  c("Intercept", "Capitalist core", "Revolutionary", "Socialist", paste0("Log(", x.name, ")"),
                                  paste0("Log(", x.name, ") $\\times$ capitalist core"), paste0("Log(", x.name, ") $\\times$ revolutionary"),
                                  paste0("Log(", x.name, ") $\\times$ socialist")),
            keep.stat = c("n", "adj.rsq"),
            add.lines = switch(is.null(model_sat)+1, list(c("Saturation constant",
                                                            "",
                                                            signif(sat,4)
            )), NULL),
            star.char = c("+","*","**","***"),
            star.cutoffs =c(0.1, 0.05, 0.01, 0.001),
            notes = c("$+$ p$<$0.1; $^*$ p$<$0.05; $^{**}$ p$<$0.01; $^{***}$ p$<$0.001"),
            notes.append = F,
            label = "tab:main")
  
  
  # statement which model is better
  if(!is.null(model_sat)){
    better_model <- ifelse(summary(model)$adj.r.squared > summary(model_sat)$adj.r.squared, 
                           "logarithmic model",
                           "saturation model")
    
    
    writeLines(paste0("\\noindent The ", better_model, " provides a better fit to the data according to the adjusted $R^2$ values."), 
               con = paste0("../figures/supplementary_information/", pred, "_", resp, "/better_model.tex"))
  }
  
  # Marginal effect plots
  plot_arguments <- list(data = data,
                         category_name = "regime_type",
                         title = FALSE,
                         rug = TRUE,
                         color_rug = TRUE,
                         twoways = FALSE,
                         rugwidth = 0.5,
                         unlog_x = TRUE,
                         log_scale_x_axis = TRUE,
                         jitter_factor = 0.2,
                         font_family = "Computer Modern Sans",
                         font_size = 9)
  
  # Log plots
  # socialist - capitalist periphery
  do.call(ggintplot, c(list(
    obj = model,
    vcov = vcov,
    varnames = c("regime_typesocialist", paste0("log(", pred, ")")
    ),
    varlabs=c("socialist regime type", x.label),
    ref = "capitalist_periphery",
    ref_lab = "capitalist periphery regime type",
    color_rug_ref = "#DDAA33",
    color_rug_v1 = "#BB5566"), plot_arguments)) +
    {if(reverse) scale_y_reverse()}
  ggsave(file = paste0("../figures/supplementary_information/", pred, "_", resp, "/marginal_effect_socialist_periphery.pdf"),
         width = 160, height  = 115, units = "mm", device=cairo_pdf) 
  
  
  # capitalist core - capitalist periphery
  do.call(ggintplot, c(list(
    obj = model,
    vcov = vcov,
    varnames = c("regime_typecapitalist_core", paste0("log(", pred, ")")
    ),
    varlabs=c("capitalist core regime type", x.label),
    ref = "capitalist_periphery",
    ref_lab = "capitalist periphery regime type",
    color_rug_ref = "#DDAA33",
    color_rug_v1 = "#004488"), plot_arguments)) +
    {if(reverse) scale_y_reverse()}
  ggsave(file = paste0("../figures/supplementary_information/", pred, "_", resp, "/marginal_effect_core_periphery.pdf"),
         width = 160, height  = 115, units = "mm", device=cairo_pdf)
  
  # socialist - capitalist core
  do.call(ggintplot, c(list(
    obj = model,
    vcov = vcov,
    varnames = c("regime_typesocialist", paste0("log(", pred, ")")
    ),
    varlabs=c("socialist regime type", x.label),
    alt_ref = "regime_typecapitalist_core",
    ref_lab = "capitalist core regime type",
    color_rug_ref = "#004488",
    color_rug_v1 = "#BB5566"), plot_arguments)) +
    {if(reverse) scale_y_reverse()}
  ggsave(file = paste0("../figures/supplementary_information/", pred, "_", resp, "/marginal_effect_socialist_core.pdf"),
         width = 160, height  = 115, units = "mm", device=cairo_pdf)
  
  # Sat plots
  if(!is.null(model_sat)){
    # socialist - capitalist periphery
    do.call(ggintplot, c(list(
      obj = model_sat,
      vcov = vcov_sat,
      varnames = c("regime_typesocialist", paste0("log(", pred, ")")
      ),
      varlabs=c("socialist regime type", x.label),
      ref = "capitalist_periphery",
      ref_lab = "capitalist periphery regime type",
      color_rug_ref = "#DDAA33",
      color_rug_v1 = "#BB5566",
      exp_y = FALSE), plot_arguments)) +
      scale_y_reverse()
    ggsave(file = paste0("../figures/supplementary_information/", pred, "_", resp, "/marginal_effect_socialist_periphery_sat.pdf"),
           width = 160, height  = 115, units = "mm", device=cairo_pdf)
    
    
    # capitalist core - capitalist periphery
    do.call(ggintplot, c(list(
      obj = model_sat,
      vcov = vcov_sat,
      varnames = c("regime_typecapitalist_core", paste0("log(", pred, ")")
      ),
      varlabs=c("capitalist core regime type", x.label),
      ref = "capitalist_periphery",
      ref_lab = "capitalist periphery regime type",
      color_rug_ref = "#DDAA33",
      color_rug_v1 = "#004488",
      exp_y = FALSE), plot_arguments)) +
      scale_y_reverse()
    ggsave(file = paste0("../figures/supplementary_information/", pred, "_", resp, "/marginal_effect_core_periphery_sat.pdf"),
           width = 160, height  = 115, units = "mm", device=cairo_pdf)
    
    # socialist - capitalist core
    do.call(ggintplot, c(list(
      obj = model_sat,
      vcov = vcov_sat,
      varnames = c("regime_typesocialist", paste0("log(", pred, ")")
      ),
      varlabs=c("socialist regime type", x.label),
      alt_ref = "regime_typecapitalist_core",
      ref_lab = "capitalist core regime type",
      color_rug_ref = "#004488",
      color_rug_v1 = "#BB5566",
      exp_y = FALSE), plot_arguments)) +
      scale_y_reverse()
    ggsave(file = paste0("../figures/supplementary_information/", pred, "_", resp, "/marginal_effect_socialist_core_sat.pdf"),
           width = 160, height  = 115, units = "mm", device=cairo_pdf)
  }
  
  # Countries excluded
  countries_used <- data %>% pull(country_name)
  countries_excluded <- all_countries_considered %>% 
    filter(country_name %in% outersect(all_countries_considered %>% pull(country_name), countries_used))
  countries_to_tex(countries_excluded, paste0("../figures/supplementary_information/", pred, "_", resp, "/countries_excluded.tex"))
  
  
  
}
