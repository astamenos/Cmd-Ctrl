#----- JSD Function -----
jsd_im <- function(im1, im2, base = 2) {
  
  # Extract pixel values
  v1 <- as.vector(im1$v)
  v2 <- as.vector(im2$v)
  
  # Remove NA values (optional but safer, especially for density maps)
  valid <- !is.na(v1) & !is.na(v2)
  v1 <- v1[valid]
  v2 <- v2[valid]
  
  # Normalize so they are probability distributions (sum to 1)
  v1 <- v1 / sum(v1)
  v2 <- v2 / sum(v2)
  
  # Avoid zero entries to prevent log(0) problems
  v1[v1 == 0] <- 1e-12
  v2[v2 == 0] <- 1e-12
  
  # Mixture
  M <- 0.5 * (v1 + v2)
  
  # Compute KLD parts
  KLD1 <- sum(v1 * (log(v1, base = base) - log(M, base = base)))
  KLD2 <- sum(v2 * (log(v2, base = base) - log(M, base = base)))
  
  # JSD
  JSD <- sqrt(0.5 * (KLD1 + KLD2))
  
  return(JSD)
}

#----- Exponential correlation -----
exp_corr <- function(dists, phi) {
  temp <- exp(-dists * phi)
  attr(temp, 'dimnames') <- NULL
  temp
}

#----- Covariate map -----
plot_covariate_map <- function(fill_var,
                               data_sf,
                               title = "Map",
                               legend_title = "Value",
                               fill_option = "inferno",
                               facet_period = FALSE) {

  p <- ggplot() +
    # Main fill layer
    geom_sf(data = data_sf, aes_string(fill = fill_var)) +

    # Optional facet wrap by period
    { if (facet_period) facet_wrap(~period) else NULL } +
    
    # Bridge overlap fill
    geom_sf(data = filter(data_sf, touches_bridge == TRUE), aes_string(fill = fill_var)) +
    
    # Contextual buffers
    geom_sf(data = st_buffer(water_sf, dist = 0.07), linewidth = 0.5, color = 'black',
            fill = 'steelblue', alpha = 0.9) +
    geom_sf(data = st_buffer(roads_sf, dist = 0.02), linewidth = 0.3, color = 'black',
            fill = 'darkgreen', alpha = 1) +
    geom_sf(data = st_buffer(mainroads_sf, dist = 0.07), linewidth = 0.5, color = 'black',
            fill = 'darkgreen', alpha = 1) +
    geom_sf(data = police_pols, color = 'white', stroke = 0.6, fill = 'darkblue',
            inherit.aes = FALSE) +
    
    # Color scale
    scale_fill_viridis_c(name = legend_title, option = fill_option) +
    
    # Theme and styling
    theme_void() +
    theme(
      legend.position = if (facet_period) 'right' else 'bottom',
      legend.box = if (facet_period) 'horizontal' else 'vertical',
      legend.direction = if (facet_period) 'vertical' else 'horizontal',
      plot.title = element_text(hjust = 0.5, size = 22, face = 'bold'),
      plot.subtitle = element_text(hjust = 0.5, size = 22, face = 'bold'),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size = 8, face = 'bold'),
      legend.text = element_text(size = 16)
    ) +
    labs(title = title)
  
  return(p)
}

#----- Check Problem Knots -----
check_gp_covariate_collinearity <- function(B, X, knots, window_sf,
                                            threshold = 0.7,
                                            title = "GP Basis: Problematic Knots Highlighted") {
  stopifnot(nrow(B) == nrow(X))  # dimensions must match
  
  # Center matrices
  B_centered <- scale(B, center = TRUE, scale = FALSE)
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  
  # Correlation: basis functions vs covariates
  cor_mat <- cor(B_centered, X_centered)
  
  # Identify problematic basis functions (columns of B) based on max correlation
  prob_z1_idx <- which(apply(abs(cor_mat), 1, function(x) any(x > threshold)))
  
  if (length(prob_z1_idx) == 0) {
    message("No problematic GP basis elements found with correlation > ", threshold)
  } else {
    cat(round(100 * length(prob_z1_idx) / ncol(B), 2), "% of knots are problematic\n",
        "Problematic basis indices:", prob_z1_idx, "\n")
  }
  
  # Extract coordinates
  if (inherits(knots, "sf")) {
    knot_coords <- st_coordinates(knots)
    knot_df <- knots
    knot_df$problematic <- FALSE
    knot_df$problematic[prob_z1_idx] <- TRUE
  } else {
    knot_coords <- as.data.frame(knots)
    knot_coords$problematic <- FALSE
    knot_coords$problematic[prob_z1_idx] <- TRUE
  }
  
  # Base plot
  p <- ggplot() + 
    geom_sf(data = city_poly, fill = 'white', linewidth = 2) +
    geom_sf(data = safe_zone, fill = 'darkgreen', alpha = 0.4) +
    
    geom_sf(data = covariate_buffer, fill = 'darkred', alpha = 0.4) +
    geom_sf(data = water_sf, fill = 'steelblue', alpha = 0.9, linewidth = 0.7) +
    geom_point(data = knot_coords, aes(x = x, y = y), 
               color = 'black', size = 1.5, shape = 17) +
    geom_point(data = knot_coords[prob_z1_idx, , drop = FALSE], 
               aes(x = x, y = y), 
               color = 'darkred', size = 3.5, shape = 17) +
    theme_minimal() +
    labs(title = title,
         subtitle = bquote(.(length(prob_z1_idx)) ~ "out of" ~ .(n_K) ~ "basis functions exceed" ~ "|" * rho * "|" > .(threshold))) +
    theme(
      plot.title = element_text(hjust = 0.5, face = 'bold', size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  print(p)
  
  return(list(
    correlation_matrix = cor_mat,
    problematic_indices = prob_z1_idx,
    knot_coords = knot_coords
  ))
}

#----- Initial Value Generator -----
inits_fn <- function() {
  list(
    beta = rnorm(p_int, beta_MLE, sd=0.1),
    gamma = rnorm(p_mark, gamma_MLE, sd=0.1),
    u1   = rnorm(n_periods, 0, 1),
    u2   = rnorm(n_periods, 0, 1),
    z1 = rnorm(n_K, 0, 1),
    log_sigma = rnorm(1, constants_list$mu_sigma, 0.1),
    log_tau = rnorm(2, constants_list$mu_tau, 0.1)
  )
}

#----- Test Run -----
test_run <- function(model_code, data_list, constants_list, cov_beta, cov_gamma, n_iter = 5000, n_burnin = 0, n_chains = 1) {
  
  inits <- inits_fn()
  
  # Build model 
  model <- nimbleModel(model_code,
                       data = data_list,
                       constants = constants_list,
                       inits = inits)
  
  config <- configureMCMC(model)
  
  # Remove defaults
  config$removeSampler(c('z1' #, 'z2'
  )) 
  config$removeSampler(c('u1', 'u2')) 
  config$removeSampler('beta')
  config$removeSampler('gamma')
  config$removeSampler('log_sigma')
  config$removeSampler('log_tau')

  config$addSampler(target = c("u1"),
                    type   = 'AF_slice',
                    control = list(
                      adaptive      = TRUE,
                      adaptInterval =  500
                    ))
  
  config$addSampler(target = c("u2"),
                    type   = 'AF_slice',
                    control = list(
                      adaptive      = TRUE,
                      adaptInterval =  500
                    ))
  
  config$addSampler(target = c('log_sigma', 'z1'),
                    type   = 'RW_block',
                    control = list(propCov = diag(1+constants_list$n_K),
                                   adaptive      = TRUE,
                                   adaptInterval =  500,
                                   scale         = 0.01
                    ))
  
  config$addSampler(
    target = c('beta'),
    type   = 'RW_block',
    control = list(propCov  = cov_beta,
                   adaptive      = TRUE,
                   adaptInterval =  500,
                   scale         = 0.05
    )
  )
  
  config$addSampler(
    target = c('gamma'), 
    type = 'RW_block', 
    control = list(propCov       = cov_gamma, 
                   adaptive      = TRUE,
                   adaptInterval = 500,
                   scale         = 0.05
    )
  )
  
  config$addSampler(
    target = c('log_tau[1]', 'log_tau[2]'),
    type   = 'RW_block',
    control = list(adaptive      = TRUE,
                   adaptInterval = 500)
  )
  
  # Adding monitors for derived quantities of interest
  config$addMonitors(c('lambda_sV_dom', 'lambda_sV_nondom', 
                       'delta_1', 
                       'delta_2',
                       'W1star[1]', 'W1star[5]', 'W1star[13]',
                       'sigma', 'tau'
  ))
  
  mcmc <- buildMCMC(config)
  compiled <- compileNimble(model, mcmc, showCompilerOutput = TRUE)
  
  system.time({
    samples_list <- runMCMC(compiled$mcmc,
                            inits = inits_fn(),
                            nburnin = n_burnin,
                            niter = n_iter,
                            thin = 1,
                            nchains = n_chains,
                            samplesAsCodaMCMC = TRUE,
                            progressBar = TRUE)
  })
  
  samples <- coda::mcmc.list(samples_list)
  rm(samples_list)
  
  params_beta  <- paste0('beta[', 1:p_int, ']')
  params_gamma <- paste0('gamma[', 1:p_mark, ']')
  params_z1    <- paste0('z1[', 1:n_K, ']')
  params_varcomp <- c('sigma', 'tau[1]', 'tau[2]')
  
  ##--- Updated proposal covariances ---##
  samples_comb <- do.call(rbind, samples)
  cov_beta_new <- cov(samples_comb[, params_beta])
  cov_gamma_new <- cov(samples_comb[, params_gamma])
  cov_GP_new   <- cov(samples_comb[, c('log_sigma', params_z1)])
  
  ##--- Posterior checks ---##
  if(n_chains > 1) {
    gelman_rubin <- coda::gelman.diag(samples[, c(params_beta, params_gamma)], multivariate = FALSE)
    print(gelman_rubin)
  }
  
  plot(samples[, c(params_beta, params_gamma)], density = FALSE)
  plot(samples[, params_varcomp], density = FALSE)
  
  model_summary <- tibble(
    metric = c('Knots', 'Integration Points', 'GP1 Effective Range', 'Burnin', 'Total Iterations'),
    value = c(n_K, nrow(ips_locs), 3 / phi_1, n_burnin, n_iter),
    type = c('int', 'int', 'num', 'int', 'int')
  )
  print(model_summary[,-3])
  
  print(summary(samples[, c(params_beta, params_gamma, params_varcomp)]))
  
  posterior_intensity(samples_comb, model_summary)
  
  return(list(
    compiled   = compiled,
    samples    = samples,
    cov_beta   = cov_beta_new,
    cov_gamma   = cov_gamma_new,
    cov_GP     = cov_GP_new
  ))
}

#----- Posterior Intensity Map -----
posterior_intensity <- function(samples_comb, model_summary) {
  # Summary label for plotting purposes
  label_df <- model_summary %>%
    mutate(label = if_else(type == 'num',
                           sprintf(paste0(metric, ': %.2f'), value),
                           sprintf(paste0(metric, ': %d'), as.integer(value))))
  
  # Storing posterior estimates for the intensities and mark probabilities
  post_mean_lambda_dom <- colMeans(samples_comb[,grep('lambda_sV_dom', colnames(samples_comb))])
  post_mean_lambda_nondom <- colMeans(samples_comb[,grep('lambda_sV_nondom', colnames(samples_comb))])
  
  to_plot <- st_drop_geometry(ips_sf) %>%
    mutate(Domestic = post_mean_lambda_dom,
           `Non-Domestic` = post_mean_lambda_nondom)
  
  ggplot(to_plot) + 
    geom_point(aes(x=log(`Non-Domestic`), y=log(Domestic), color = period)) +
    labs(x = expression(lambda['0']),  y = expression(lambda['1']))
  
  to_plot <- to_plot %>%
    pivot_longer(names_to = 'Offense', cols = c('Domestic', 'Non-Domestic'), values_to = 'offense_lambda')
  
  # Add a new variable to indicate 'Observed' vs 'Model'
  obs_pp$source <- 'Observed'
  to_plot$source <- 'Model'
  to_plot$N <- NA
  
  
  # Combine the observed points and ips
  combined_data <- bind_rows(
    obs_pp %>% select(x, y, period, Offense, offense_lambda, source),
    to_plot %>% select(x, y, period, Offense, offense_lambda, source)
  )
  
  # Ensure `to_plot` is split into Domestic and Non-Domestic subsets
  
  to_plot_domestic <- to_plot %>%
    filter(Offense == 'Domestic', period %in% periods) %>%
    #mutate(geometry_str = st_as_text(geometry)) %>%
    #st_drop_geometry() %>%
    rename(lambda_domestic = offense_lambda) %>%
    left_join(clipped) %>%
    st_as_sf() %>%
    st_set_geometry('geometry')
  
  #mutate(centroid = st_centroid(geometry))
  
  to_plot_nondomestic <- to_plot %>% 
    filter(Offense == 'Non-Domestic', period %in% periods) %>%
    #st_drop_geometry() %>%
    rename(lambda_nondomestic = offense_lambda) %>%
    left_join(clipped) %>%
    st_as_sf() %>%
    st_set_geometry('geometry')
  #mutate(centroid = st_centroid(geometry))
  
  obs_pp_filtered <- obs_pp %>% 
    filter(period %in% periods) %>%
    slice_sample(by = Offense, prop = 0.7) %>%
    mutate(plt.alpha = if_else(Offense == 'Domestic', 0.5, 0.4),
           plt.size = if_else(Offense == 'Domestic', 0.5, 0.4))
  rm(to_plot)
  
  # Mean posterior intensity surface by offense type and year
  plt <- ggplot() +
    # Non-Domestic
    geom_sf(data = to_plot_nondomestic, 
            aes(fill = 1*(lambda_nondomestic))) +
    scale_fill_viridis_c(
      name = 'Non-Domestic λ(s, t)', option = 'inferno',
      #breaks = seq(0, max(to_plot_nondomestic$lambda_nondomestic, na.rm = TRUE), by = 300),
      aesthetics = 'fill',
      guide = guide_colorbar(order = 2, title.position = 'top', barwidth = 20, barheight = 2)
    ) +
    new_scale_fill() +
    
    # Domestic
    geom_sf(data = to_plot_domestic, 
            aes(fill = 1*(lambda_domestic))) +
    scale_fill_viridis_c(
      name = 'Domestic λ(s, t)', option = 'inferno',
      #breaks = seq(0, max(to_plot_domestic$lambda_domestic, na.rm = TRUE), by = 30),
      aesthetics = 'fill',
      guide = guide_colorbar(order = 1, title.position = 'top', barwidth = 20, barheight = 2)
    ) +
    
    geom_point(data = obs_pp_filtered, aes(x = x, y = y, alpha = plt.alpha, size = plt.size),
               color = 'white') +
    
    # Landmarks
    geom_sf(data = st_buffer(water_sf, dist = 0.07), linewidth = 0.5, color = 'black', 
            fill = 'steelblue', alpha = 1) +
    geom_sf(data = filter(ips_sf, touches_bridge==T), aes(fill = N), color = NA, size = 0) +
    geom_sf(data=police_pols,
            color = 'white', stroke = 0.6, fill = 'darkblue', inherit.aes = FALSE) +
    geom_sf(data = st_buffer(roads_sf, dist = 0.02), linewidth = 0.3, color = 'black', 
            fill = 'darkgreen', alpha = 1) +
    geom_sf(data = st_buffer(mainroads_sf, dist = 0.07), linewidth = 0.5, color = 'black', 
            fill = 'darkgreen', alpha = 1) +

    scale_alpha_identity(guide = 'none') +
    scale_size_identity(guide = 'none') +
    scale_size_continuous(
      name = 'Assaults per Location (N)',
      range = c(0.1, 0.5),
      breaks = pretty(obs_pp_filtered$N, n = 5),
      guide = guide_legend(title.position = 'top', title.hjust = 0.5)
    ) +
    
    #geom_sf(data = my_acs, linewidth = 0.5, fill = adjustcolor('white', alpha = 0)) +
    
    facet_grid(Offense ~ period, switch = 'y') +
    
    theme_minimal() +
    theme(
      legend.position = 'right',
      legend.box = 'vertical',  # Required to get side-by-side colorbars
      legend.direction = 'horizontal',  # Makes both bars horizontal
      plot.title = element_text(hjust = 0.5, size = 22, face = 'bold'),
      plot.subtitle = element_text(hjust = 0.5, size = 22, face = 'bold'),
      strip.text.x = element_text(size = 20, face = 'bold', color = 'white'),
      strip.text.y = element_text(size = 20, face = 'bold', color = 'white'),
      strip.background = element_rect(fill = 'black', color = NA),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size = 18, face = 'bold'),
      legend.text = element_text(size = 16)
    ) +
    labs(title = 'Offense-Specific Intensity λ(s, t)', x = NULL, y = NULL,
         subtitle = bquote(phi[1] == .(round(phi_1, 3))
                        )
         )
  
  
  plt <- ggdraw(plt) +
    draw_label(
      label_df$label[1], x = plt_style$x, y = 0.80, hjust = plt_style$hjust, vjust = plt_style$vjust, 
      size = plt_style$size, fontface = plt_style$fontface
    ) +
    draw_label(
      label_df$label[2], x = plt_style$x, y = 0.75, 
      hjust = plt_style$hjust, vjust = plt_style$vjust, 
      size = plt_style$size, fontface = plt_style$fontface
    ) +
    draw_label(
      label_df$label[3], x = plt_style$x, y = 0.70, 
      hjust = plt_style$hjust, vjust = plt_style$vjust, 
      size = plt_style$size, fontface = plt_style$fontface
    ) +
    draw_label(
      label_df$label[4], x = plt_style$x, y = 0.30, 
      hjust = plt_style$hjust, vjust = plt_style$vjust, 
      size = plt_style$size, fontface = plt_style$fontface
    ) +
    draw_label(
      label_df$label[5], x = plt_style$x, y = 0.25, 
      hjust = plt_style$hjust, vjust = plt_style$vjust, 
      size = plt_style$size, fontface = plt_style$fontface
    ) +
    draw_label(
      label_df$label[6], x = plt_style$x, y = 0.20, 
      hjust = plt_style$hjust, vjust = plt_style$vjust, 
      size = plt_style$size, fontface = plt_style$fontface
    )
  print(plt)  
}

#----
visualize_problematic_knots <- function(samples_comb, knots, window_sf,
                                        threshold = 0.7) {
  
  # Only keep relevant parameters
  fixed_effects <- c('log_sigma', paste0('beta[', 1:p_int, ']'))
  z1_names <- paste0('z1[', 1:n_K, ']')
  all_names <- c(fixed_effects, z1_names)
  samples_sub <- samples_comb[, all_names]
  
  # Compute posterior correlation matrix
  cor_mat <- cor(samples_sub)
  
  # Get problematic z1 indices (correlation > threshold with any fixed effect)
  prob_z1_idx <- which(apply(cor_mat[fixed_effects, z1_names, drop = FALSE], 2, function(x) any(abs(x) > threshold)))
  
  if (length(prob_z1_idx) == 0) {
    message("No problematic GP basis elements found with correlation > ", threshold)
  } else {
    cat(round(100*length(prob_z1_idx)/n_K, 2), "% of knots are problematic\n", "Problematic basis indices (z1):", prob_z1_idx, "\n")
  }
  
  # Extract coordinates if sf
  if (inherits(knots, "sf")) {
    knot_coords <- st_coordinates(knots)
  } else {
    knot_coords <- knots
  }
  
  # Plotting
  p <- ggplot() + 
    geom_sf(data = city_poly, fill = 'white', linewidth = 2) +
    geom_sf(data = safe_zone, fill = 'darkgreen', alpha = 0.3) +
    
    geom_sf(data = covariate_buffer, fill = 'darkred', alpha = 0.4) +
    geom_sf(data = water_sf, fill = 'steelblue', alpha = 0.9, linewidth = 0.7) +
    
    geom_sf(data=police_pols,
            color = 'white', stroke = 0.6, fill = 'darkblue', inherit.aes = FALSE, alpha = 0.8) +
    geom_point(data = knots, aes(x=x, y=y), color = 'black', size = 1.5, shape = 17) +
    geom_point(data = knots[prob_z1_idx,], aes(x=x, y=y), color = 'darkred', size = 3.5, shape = 17) +
    theme_minimal() +
    theme(
      legend.position = 'bottom',
      legend.box = 'vertical',  # Required to get side-by-side colorbars
      legend.direction = 'horizontal',  # Makes both bars horizontal
      plot.title = element_text(hjust = 0.5, size = 22, face = 'bold'),
      plot.subtitle = element_text(hjust = 0.5, size = 22, face = 'bold'),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size = 18, face = 'bold'),
      legend.text = element_text(size = 16)
    ) +
    labs(title = "GP Basis: Problematic Knots Highlighted", 
         subtitle = bquote(.(length(prob_z1_idx)) ~ "out of" ~ .(n_K) ~ "basis functions exceed" ~ "|" * rho * "|" > .(threshold) ~' ,'~ 
                           n[K] == .(round(n_K, 0)) ~' ,'~ 
                           r == .(round(r, 2)) ~' ,'~ 
                           d(knots) == .(round(med_knot_dist, 2)) ~' km')
         )
  
  print(p | p2 + geom_point(data = knots[prob_z1_idx,], aes(x=x, y=y), color = 'darkred', size = 5.5, shape = 17) 
        | p3 + geom_point(data = knots[prob_z1_idx,], aes(x=x, y=y), color = 'darkred', size = 5.5, shape = 17) 
        | p4 + geom_point(data = knots[prob_z1_idx,], aes(x=x, y=y), color = 'darkred', size = 5.5, shape = 17))
  
  # Return useful objects
  return(list(correlation_matrix = cor_mat,
              problematic_indices = prob_z1_idx,
              knot_coords = knot_coords))
}

# LGCP log-likelihood for Poisson count with binomial marks
d_lgcp <- nimbleFunction(
  run = function(x         = double(1),
                 lambda_sV = double(1),
                 lambda_s  = double(1),
                 E_V         = double(1),
                 #volume    = double(0),
                 #n_V       = double(0),
                 log       = integer(0, default = 1)) {
    returnType(double())
    lambda_DT <- sum(E_V*lambda_sV)
    val       <- sum(log(lambda_s)) - lambda_DT
    if(log)    return(val)
    else       return(exp(val))
  }
)

# LGCP log-likelihood for Poisson count with binomial marks
r_lgcp <- nimbleFunction(
  run = function(n         = integer(0),
                 lambda_sV = double(1),
                 lambda_s  = double(1),
                 E_V         = double(1)
                 #volume    = double(0),
                 #n_V       = double(0)
                 ) {
    returnType(double(1))
    ## here you must return a length‐n_points integer vector of simulated s's 
    ## (for WAIC you only need the log‐lik, so you can simply return the observed s
    ##  or a dummy, e.g. all zeros, if you don’t care about the *replicate*)
    out <- rep(0.0, length(lambda_s)) 
    return(out)
  }
)

#registerDistributions(list(
#  d_lgcp = list(
#    BUGSdist = "d_lgcp(x, lambda_sV, lambda_s, N, volume, n_V)",
#    types = c("x = double(1)", 
#              "lambda_sV = double(1)", 
#              "lambda_s = double(1)", 
#              #"pi = double(1)", 
#              #"Y = double(1)", 
#              "N = double(1)", 
#              "volume = double(0)", 
#              "n_V = double(0)"),
#    discrete = FALSE
#  )
#))
