#----- Libraries -----
require(scales)
require(patchwork)
require(cowplot)
require(spdep)
require(ggnewscale)
require(parallel)
suppressPackageStartupMessages(require(spatstat))
require(sf)
require(tidyverse)
require(nimble)
require(lme4)
require(FNN)

#----- MCMC Arguments -----
niter <- 150000
nburnin <- 0.5*niter
thin <- 10
nviable <- (niter-nburnin)/thin
nchains <- 3
x_rng <- window$xrange
y_rng <- window$yrange
intpoint_knot_ratio <- 17

plt_style <- list(x = 0.65, size = 25, hjust = 0, vjust = 1, fontface = 'bold')

#----- Functions -----
# Exponential correlation function
exp_corr <- function(dists, phi) {
  temp <- exp(-dists * phi)
  attr(temp, 'dimnames') <- NULL
  temp
}

#----- Model -----
model_code <- nimbleCode({
  # LGCP likelihood
  s[1:n_points] ~ d_lgcp(
    lambda_sV = lambda_sV[1:n_V],
    lambda_s  = lambda[1:n_points],
    E_V       = area[1:n_V]
  )
  
  ## Draw non-centered GP normals
  for(k in 1:n_K) {
    z1[k] ~ dnorm(0, 1)
    
  }
  for(k in 1:n_K) {
    W1star[k] <- inprod(L1[k, 1:n_K], z1[1:n_K])
  }
  
   meanW1star <- mean(W1star[1:n_K])

   # Center W1 star
   for(k in 1:n_K) {
     W1star_c[k] <- W1star[k] - meanW1star
   }
  
  ## *Observed* unique locations
  for(k in 1:n_locations) {
    ## spatial GP tilts
    W1_loc[k]   <- sigma * inprod(C_cross_C_inv_obs_1[k,1:n_K], W1star_c[1:n_K])
    
    ## fixed-effects linear predictor
    Xbeta_loc[k]  <- inprod(X_int[k, 1:p_int], beta[1:p_int])
    
    ## combined predictor at observed loc
    eta_loc[k]   <- Xbeta_loc[k] + W1_loc[k] + log_pop_dens[k]
  }
  
  ## Iintegration point unique locations
  for(k in 1:n_A) {
    W1_loc_ips[k]   <- sigma * inprod(C_cross_C_inv_int_1[k,1:n_K], W1star_c[1:n_K])
    
    Xbeta_loc_ips[k]  <- inprod(X_int_V[k, 1:p_int], beta[1:p_int])
    eta_loc_ips[k]   <- Xbeta_loc_ips[k] + W1_loc_ips[k] + log_pop_dens_V[k]
  }
  
  # Observed-point likelihood
  for(i in 1:n_points) {
    # Intensity model
    log(lambda[i]) <- eta_loc[loc[i]] + delta_1[year[i]]
    
    # Mark model
    Y[i] ~ dbin(prob = pi[i], size = 1)
    logit(pi[i]) <- inprod(X_mark[i, 1:p_mark], gamma[1:p_mark]) + delta_2[year[i]]
  }
  
  # Integration points
  for(j in 1:n_V) {
    # Intensity model
    log(lambda_sV[j]) <- eta_loc_ips[loc_V[j]] + delta_1[year_V[j]]
    
    # Mark model
    logit(pi_V[j])    <- inprod(X_mark_V[j, 1:p_mark], gamma[1:p_mark]) + delta_2[year_V[j]]
    
    # Posterior surfaces
    lambda_sV_dom[j]    <- lambda_sV[j] * pi_V[j]
    lambda_sV_nondom[j] <- lambda_sV[j] * (1-pi_V[j])
  }

  # Temporal random effects
  for(t in 1:n_years) {
    u1[t] ~ dnorm(0, sd = 1)
    u2[t] ~ dnorm(0, sd = 1)
  }
  mean_u1 <- mean(u1[1:n_years])
  mean_u2 <- mean(u2[1:n_years])
  for(t in 1:n_years) {
    u1_c[t] <- u1[t] - mean_u1
    u2_c[t] <- u2[t] - mean_u2
    
    delta_1[t] <- tau[1]*u1_c[t]
    delta_2[t] <- tau[2]*u2_c[t]
  }
  
  # Priors for covariate coefficients
  for(k in 1:p_int) beta[k]  ~ dnorm(mu_beta[k], sd = nu_beta[k])
  for(k in 1:p_mark)  gamma[k] ~ dnorm(mu_gamma[k], sd = nu_gamma[k])

  # GP variance
  log_sigma ~ dnorm(mu_sigma, sd = nu_sigma)
  sigma     <- exp(log_sigma)
  sigma2    <- sigma^2

  # Temporal RE variance
  for(l in 1:2) {
    log_tau[l]   ~ dnorm(mu_tau[l],   sd = nu_tau[l])
    tau[l]       <- exp(log_tau[l])
  }
})

## ----- Loading Data -----
MPLS_crime <- read_csv('data/MPLS_crime.csv')
MPLS_crime_sf <- st_as_sf(MPLS_crime, 
                          coords = c('x','y'), 
                          crs=crs_use) %>%
  mutate(year = factor(year),
         GEOID = factor(GEOID)) 

#----- Point Processes -----
# Combined point process
assaults_ppp <- as.ppp(st_filter(MPLS_crime_sf, window_sf), W = window)

marks(assaults_ppp) <- with(st_filter(MPLS_crime_sf, window_sf), 
                            data.frame(Y = Y, year = year))

# Point process for domestic
dom_ppp <- as.ppp(subset(assaults_ppp, marks(assaults_ppp)$Y == 1),
                  W = window)

# Point process for non-domestic
nondom_ppp <- as.ppp(subset(assaults_ppp, marks(assaults_ppp)$Y == 0),
                     W = window)

# Densities
dens_dom <- density.ppp(subset(dom_ppp, marks(dom_ppp)$year == '2021'))
dens_nondom <- density.ppp(subset(nondom_ppp, marks(nondom_ppp)$year == '2021')) # Fix
jsd_im(dens_dom, dens_nondom)


#L_dom <- envelope(subset(dom_ppp, marks(dom_ppp) == '2021'), Lest, nsim = 10)
#plot(L_dom)

#L_nondom <- envelope(subset(nondom_ppp, marks(nondom_ppp) == '2021'), Lest, nsim = 10)
#plot(L_nondom)

#----- GP -----
##----- Knots -----
# Build a coarse, regular grid at spacing ≈ r
phi0 <- 170/100            
r       <- 3 / phi0 
dx <- r
dy <- r

coarse_grid <- expand.grid(
  x = seq(x_rng[1], x_rng[2], by=0.5),
  y = seq(y_rng[1], y_rng[2], by=0.5)
)

coarse_grid <- coarse_grid[
 inside.owin(x = coarse_grid$x,
             y = coarse_grid$y,
             w = window),] 


# Count events in each coarse cell (within radius r/2)
densities <- vapply(seq_len(nrow(st_drop_geometry(coarse_grid))), function(i) {
  cx <- coarse_grid$x[i]
  cy <- coarse_grid$y[i]
  
  # Squared distances to all observed points
  d2 <- (assaults_ppp$x - cx)^2 + (assaults_ppp$y - cy)^2
  sum(d2 <= (r/2)^2)
}, integer(1))

# Variogram
coarse_grid <- coarse_grid %>%
  st_as_sf(coords = c('x', 'y'), crs=crs_use)
h_max <- max(dist(st_coordinates(coarse_grid$geometry)))
bin_width <- h_max / 20             

glm_mod <- glm(densities~1, family = 'poisson', data=st_drop_geometry(coarse_grid))
coarse_grid$pearson_resids <- resid(glm_mod, type = 'pearson')
emp_vg <- variogram(pearson_resids ~ 1, data = coarse_grid, cutoff = 6, width = bin_width)

vgm_fit <- fit.variogram(emp_vg, model =  vgm("Gau", range = 2.5, psill = 800, nugget = 50))
vgm_fit
plot(emp_vg, vgm_fit)
phi_G0 <- 1/vgm_fit$range[2]
phi_E0 <- phi_G0*sqrt(-log(0.05))
r <- 3/phi_E0

# Create main knots
spacing_factor <- 31/100
h_knot <- spacing_factor*r # spacing based on rule-of-thumb
x_seq <- seq(x_rng[1], x_rng[2], by = h_knot)
y_seq <- seq(y_rng[1], y_rng[2], by = h_knot)

knots <- expand.grid(x = x_seq, y = y_seq)

knots <- knots[inside.owin(x = knots$x,
                           y = knots$y,
               w = window),]

dx <- h_knot
dy <- h_knot
thresh_prob <- 0.3

# Count events in each coarse cell (within radius r/2)
densities <- vapply(seq_len(nrow(knots)), function(i) {
  cx <- knots$x[i]
  cy <- knots$y[i]
  
  # Squared distances to all observed points
  d2 <- (assaults_ppp$x - cx)^2 + (assaults_ppp$y - cy)^2
  sum(d2 <= (r/3)^2)
}, integer(1))

# Pick the “hottest” cells
thresh_dens    <- quantile(densities, thresh_prob)
to_refine <- which(densities >= thresh_dens)

# Sprinkle a small subgrid around each hot cell
extra_knots <- do.call(rbind, lapply(to_refine, function(i) {
  cx <- knots$x[i]
  cy <- knots$y[i]
  expand.grid(
    x = seq(cx - dx/2, cx + dx/2, length.out = 1),
    y = seq(cy - dy/2, cy + dy/2, length.out = 1)
  )
}))

knots <- unique(rbind(knots, extra_knots))

# Knot buffer
knots_sf <- st_as_sf(knots, coords = c("x", "y"), crs = crs_use)
city_inner <- st_buffer(city_poly, dist = -r*(20/100))
water_buffer <- st_buffer(water_sf, dist = r*(12/100)) # Expand water slightly to avoid knots near lake edges
mainroad_buffer <- st_buffer(mainroads_sf, dist = r*(8/100))
police_buffer <- st_buffer(st_union(police_pols), dist = r*(2/100))
covariate_buffer <- st_union(mainroad_buffer, police_buffer)
covariate_buffer <- st_union(water_buffer, covariate_buffer)
safe_zone <- st_difference(city_inner, covariate_buffer) # Create a "safe zone" = inside shrunken city AND not near buffered water

# Keep knots that fall inside the safe zone
knots <- knots[st_within(knots_sf, safe_zone, sparse = FALSE), ]

# Distances between knots
n_K <- nrow(knots)
dists_knots     <- as.matrix(dist(knots))
d_nn <- apply(dists_knots, 1, function(r) min(r[r>0]))
med_knot_dist <- median(d_nn)

# Knot placement
ggplot() + 
  geom_sf(data = city_poly, fill = 'white', linewidth = 2) +
  geom_sf(data = safe_zone, fill = 'darkgreen', alpha = 0.2) +

  geom_sf(data = covariate_buffer, fill = 'darkred', alpha = 0.4) +
  geom_sf(data = water_sf, fill = 'steelblue', alpha = 0.9, linewidth = 0.7) +
  #geom_sf(data = st_buffer(police_pols, dist=-0.1), fill = 'darkblue', alpha = 0.9, linewidth = 0.7) +
  geom_sf(data = roads_sf, linewidth = 0.7, color = 'black', 
          fill = 'black', alpha = 1) +
  geom_sf(data = mainroads_sf, linewidth = 1, color = 'black', 
          fill = 'black', alpha = 1) +
  geom_point(data = knots, aes(x=x, y=y)) +
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
  labs(title = 'Knot Placement', 
       subtitle = bquote(n[K] == .(round(nrow(knots), 0)) ~' ,'~ 
                         r == .(round(r, 2)) ~' ,'~ 
                         d(knots) == .(round(med_knot_dist, 2)) ~' km'
       )
  )

##----- Integration Points -----
A <- area(window)
cell_size_1 <- r/2                                      # Resolve the GP
cell_size_2 <- sqrt(A/(n_K*intpoint_knot_ratio))
cell_size <- min(cell_size_1, cell_size_2)

n_x   <- ceiling(diff(x_rng) / cell_size)
n_y   <- ceiling(diff(y_rng) / cell_size)
ips <- expand.grid(
  x = seq(x_rng[1], x_rng[2], length.out = n_x),
  y = seq(y_rng[1], y_rng[2], length.out = n_y)
)

ips <- ips[inside.owin(x = ips$x, 
                       y = ips$y, 
                       w = window), ]

# Distinct integration point locations
ips_locs <- ips %>% distinct(x, y)
nrow(ips_locs)

# A bunch of spatial operations to remove integration points from water but still connect East + West Bank through bridges
land_minus_water <- st_difference(window_sf, water_sf)
land_patches <- st_cast(land_minus_water, "POLYGON") # Break into individual polygons
patch_areas <- st_area(land_patches) # Find the “main” patch (the one with the largest area)
areas <- as.numeric(st_area(land_patches)) # Compute their areas

# Pick the top 2 biggest patches (East & West Bank)
keep_idxs <- order(areas, decreasing = TRUE)[1:2]
main_land  <- land_patches[keep_idxs, ]

# Build raw polygon grid over the study window
grid_polys <- st_make_grid(
  main_land,
  cellsize = c(cell_size, cell_size),
  what     = "polygons",
  square   = TRUE
)
grid_sf <- st_sf(
  grid_id  = seq_along(grid_polys),
  geometry = grid_polys
) %>% 
  mutate(full_area = as.numeric(st_area(geometry)))

# Clip to the land‐only window (no bridges yet)
clipped <- st_intersection(grid_sf, window_sf) %>%
  mutate(clipped_area = as.numeric(st_area(geometry)))

# Pull out exactly those bits of roads that lie over water
bridge_polys <- st_intersection(roads_sf, water_sf)
bridges_sf <- bridge_polys %>% st_cast("MULTILINESTRING")

# Rescue any clipped grid‐cell whose patch overlaps one of these
rescue_ids <- clipped$grid_id[
  lengths(st_intersects(clipped, bridges_sf)) > 0
]

# 5a) Rescue those cells by giving them back their full grid square
rescue_sf <- grid_sf %>% 
  filter(grid_id %in% rescue_ids) %>%
  mutate(
    clipped_area   = full_area,
    touches_bridge = TRUE,
    was_clipped    = FALSE
  )
keep_sf <- clipped %>% 
  filter(!grid_id %in% rescue_ids) %>%
  mutate(
    touches_bridge = FALSE,
    was_clipped    = clipped_area < full_area
  )


clipped2 <- bind_rows(rescue_sf, keep_sf) %>%
  mutate(centroid = st_centroid(geometry)) %>%
  filter(
    lengths(st_intersects(centroid, main_land)) > 0 | touches_bridge == T
  )

## Get centroids and pull out x,y,E_V into a plain tibble
ips_locs <- clipped2 %>%
  st_centroid() %>%
  transmute(
    grid_id,
    x   = st_coordinates(.)[,1],
    y   = st_coordinates(.)[,2],
    E_V = clipped_area
  ) %>%
  st_drop_geometry()

# Expand integration point grid by year
ips <- expand_grid(ips_locs, year = years_fac) %>%
  mutate(loc_V = row_number())

## Join *every* crime‐point to its containing grid cell
events_with_cell <- st_join(
  MPLS_crime_sf,            # your sf of crime‐points
  clipped2 %>% select(grid_id),
  join = st_intersects,     # so points on a boundary still match
  left = FALSE              # drop any that truly fall outside
) %>%
  st_drop_geometry() %>%
  mutate(
    # find the row in `ips` for this (grid_id, year)
    loc = match(
      paste0(grid_id, "_", year),
      paste0(ips$grid_id, "_", ips$year)
    )
  )

## Aggregate at each cell×year
counts_df <- events_with_cell %>%
  group_by(grid_id, year) %>%
  summarise(
    N = sum(N),
    Y = sum(Y),   
    domestic_flag = mean(domestic_flag),
    .groups = "drop"
  )

## Attach those counts back to `ips` (zero if none)
ips <- ips %>%
  left_join(counts_df, by = c("grid_id","year")) %>%
  mutate(
    N = replace_na(N, 0),
    Y = replace_na(Y, 0),
    domestic_flag = replace_na(domestic_flag, 0)
  )

ips_sf <- st_as_sf(ips, coords = c('x','y'), crs = st_crs(MPLS_crime_sf))

# Spatially join each centroid to the one (and only one) patch polygon
# Then join the summary‐columns *into* clipped2, but with the GIS operation:
ips_sf <- clipped2 %>% 
  # Cross each patch‐polygon with each year
  expand_grid(year = years_fac) %>%
  
  # Join counts by grid_id & year
  left_join(counts_df, by = c("grid_id","year")) %>%
  mutate(
    N             = replace_na(N, 0),
    Y             = replace_na(Y, 0),
    domestic_flag = replace_na(domestic_flag, 0)
  ) %>%
  
  # Now convert to sf once, naming the geometry column
  st_as_sf(sf_column_name = "geometry") %>%
  
  # Compute centroids & coords & local index
  mutate(
    centroid = st_centroid(geometry),
    x        = st_coordinates(centroid)[,1],
    y        = st_coordinates(centroid)[,2],
    loc_V    = row_number()
  )

ips <- st_drop_geometry(ips_sf)
ips_coords <- st_coordinates(ips_sf$centroid)
ips_locs <- ips %>% distinct(centroid) %>% 
  mutate(x = st_coordinates(centroid)[,'X'], 
         y = st_coordinates(centroid)[,'Y'])

# Some summaries
n_A <- nrow(ips_locs)
n_V <- nrow(ips)
A_diff <- 100*(sum(ips$clipped_area)/3-A)/A
sprintf('# of integration points (per year): n_A = %d', n_A)
sprintf('# of integration points (across years): n_V = %d', n_V)
sprintf('R = n_A:n = %f', n_A/n_K)
sprintf('%% Diff between A=|D| and area of integration points: %.2f %%', 
        A_diff)

##----- d-Nearest Neighbors ----- 
# To smooth point-level covariates and create a % domestic flag spatial covariate
flag_nn_vec <- ips$domestic_flag

# Search radius
dx <- mean(diff(sort(unique(ips$x))))
dy <- mean(diff(sort(unique(ips$y))))
d  <- med_knot_dist/2

for(t in years) {
  # Subset the integration points in year t
  ips_t   <- ips_sf %>% filter(year == t)
  ids_t   <- ips_t$loc_V
  coords  <- st_coordinates(ips_t$centroid)
  
  # Precompute pairwise distances among the ips for year t
  Dmat <- as.matrix(dist(coords))
  
  for(ii in seq_along(ids_t)) {
    full_i <- ids_t[ii]
    
    di <- Dmat[ii, ]
    di[ii] <- Inf   # exclude self
    
    # Find all other integration‐points within radius d
    neigh <- which(di <= d)
    neigh
    if(length(neigh) > 0) {
      # average over those neighbours
      full_neigh_ids     <- ids_t[neigh]
      flag_nn_vec[full_i] <- mean(flag_nn_vec[full_neigh_ids], na.rm = TRUE)
    } else {
      # fallback to the single nearest neighbour
      j0                 <- which.min(di)
      full_j0            <- ids_t[j0]
      
      flag_nn_vec[full_i] <- flag_nn_vec[full_j0]
    }
  }
}

sum(is.na(flag_nn_vec))

# Attach point-level covariates to the integration points
ips_sf <- ips_sf %>%
  mutate(
    domestic_flag = flag_nn_vec,
    other = 1 - domestic_flag
  )

ips <- st_drop_geometry(ips_sf)

# Calculate tract-level averages
ips_sf <- st_join(
  st_set_geometry(ips_sf, 'centroid'), 
  my_acs %>% select(pop_dens, pct_pov, GEOID), 
  join = st_within,   # only keeps rows where the point falls inside a tract
  left = FALSE        # drop any integration points that don’t lie in a tract
) %>%
  st_set_geometry('geometry') %>%
  mutate(pop = pop_dens*clipped_area)

ips <- st_drop_geometry(ips_sf)

MPLS_crime_sf <- MPLS_crime_sf %>%
  mutate(temp_year = year) %>%
  select(-year, -domestic_flag) %>%
  st_join(
    ips_sf %>% select(grid_id, year, domestic_flag, pop),
    join = st_within,    
    left = FALSE,         # drop crimes outside all cells
    by = 'year'
  ) %>%
  filter(year == temp_year) %>%
  select(Precinct, Ward, GEOID, grid_id, year, Case_Number, Type, 
         Offense_Category, Offense, Y, pop, pop_dens, pct_pov, domestic_flag)


# Compute the distance from each crime point to its nearest police station:
nearest_idx <- st_nearest_feature(MPLS_crime_sf, police_pols)
dist_police_obs <- st_distance(
  MPLS_crime_sf,
  police_pols[nearest_idx, ],
  by_element = TRUE
)

nearest_idx <- st_nearest_feature(MPLS_crime_sf, water_sf)
dist_water_obs <- st_distance(
  MPLS_crime_sf,
  water_sf[nearest_idx, ],
  by_element = TRUE
)

nearest_idx <- st_nearest_feature(MPLS_crime_sf, mainroads_sf)
dist_roads_obs <- st_distance(
  MPLS_crime_sf,
  mainroads_sf[nearest_idx, ],
  by_element = TRUE
)

# Create covariate
MPLS_crime_sf <- MPLS_crime_sf %>%
  mutate(dist_police = as.numeric(dist_police_obs),
         dist_water = as.numeric(dist_water_obs),
         dist_roads = as.numeric(dist_roads_obs))

# Repeat for integration points (ips):
nearest_idx <- st_nearest_feature(ips_sf, police_pols)
dist_police_int <- st_distance(
  ips_sf,
  police_pols[nearest_idx, ],
  by_element = TRUE
)

nearest_idx <- st_nearest_feature(ips_sf, water_sf)
dist_water_int <- st_distance(
  ips_sf,
  water_sf[nearest_idx, ],
  by_element = TRUE
)

nearest_idx <- st_nearest_feature(ips_sf, mainroads_sf)
dist_roads_int <- st_distance(
  ips_sf,
  mainroads_sf[nearest_idx, ],
  by_element = TRUE
)

ips_sf <- ips_sf %>%
  mutate(dist_police = as.numeric(dist_police_int),
         dist_water = as.numeric(dist_water_int),
         dist_roads = as.numeric(dist_roads_int))
ips <- st_drop_geometry(ips_sf)

n_A <- length(unique(ips_sf$centroid))
n_V <- length(ips_sf$centroid)

#-----Visualizations -----
##----- Map of observed assaults aggragated by integration point -----
plt <- ggplot() +
  facet_wrap(~year) +
  geom_sf(data = ips_sf, aes(fill = N), color = NA, size = 0) +
  #geom_sf(data = my_acs, linewidth = 0.5, fill = adjustcolor('white', alpha = 0)) +
  geom_sf(data = st_buffer(water_sf, dist = 0.07), linewidth = 0.5, color = 'black', 
          fill = 'steelblue', alpha = 1) +
  geom_sf(data = filter(ips_sf, touches_bridge==T), aes(fill = N), color = NA, size = 0) +
  geom_sf(data = st_buffer(roads_sf, dist = 0.02), linewidth = 0.3, color = 'black', 
          fill = 'darkgreen', alpha = 1) +
  geom_sf(data = st_buffer(mainroads_sf, dist = 0.07), linewidth = 0.5, color = 'black', 
          fill = 'darkgreen', alpha = 1) +
  scale_fill_viridis_c(name = 'Number of Assaults',
                       option = 'inferno') +
  
  # the knot locations
  geom_point(data = knots,
             shape = 23,
             aes(x = x, y = y),
             fill = 'darkgreen', color = 'black', stroke = 0.5,
             size = 0.1, alpha = 0.9) +

  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.box = 'vertical',  # Required to get side-by-side colorbars
    legend.direction = 'horizontal',  # Makes both bars horizontal
    plot.title = element_text(hjust = 0.5, size = 22, face = 'bold'),
    plot.subtitle = element_text(hjust = 0.5, size = 22, face = 'bold'),
    strip.text.x = element_text(size = 20, face = 'bold', color = 'white'),
    strip.text.y = element_text(size = 20, face = 'bold', color = 'white'),
    strip.background = element_rect(fill = 'darkred', color = NA),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 18, face = 'bold'),
    legend.text = element_text(size = 16)
  ) +
  labs(title = sprintf('Number of Assaults (at %s Integration Points of Size %.2f km x %.2f km )', n_A, cell_size, cell_size),
       subtitle = bquote(phi[1] == .(round(phi0, 2)) ~','~ 
                         phi[2] == .(round(phi0, 2))  ~','~ 
                         n[A]:n[K] == .(round(n_A, 0))~':'~.(round(n_K, 0)) == .(round(n_A/n_K, 0)) ~' ,'~ 
                         r == .(round(r, 2)) ~' ,'~ 
                           Delta[A] == .(round(A_diff, 2)) ~ '%'
       ))
plt

# Covariate maps
p1 <- plot_covariate_map("domestic_flag", ips_sf,
                         title = "% of Calls that Indicate Domestic",
                         legend_title = "% of Calls",
                         fill_option = "inferno",
                         facet_year = TRUE)

p2 <- plot_covariate_map("dist_police", ips_sf,
                         title = "Distance from Nearest Police Station",
                         legend_title = "Kilometers")

p3 <- plot_covariate_map("dist_roads", ips_sf,
                         title = "Distance from Nearest Main Road",
                         legend_title = "Kilometers")

p4 <- plot_covariate_map("dist_water", ips_sf,
                         title = "Distance from Nearest Body of Water",
                         legend_title = "Kilometers")

p5 <- plot_covariate_map("pop", ips_sf,
                         title = "Population by Integration Point",
                         legend_title = "People")

plt <- p1 / (p2 | p3 | p4 | p5)
print(plt)

# Thin points for plotting purposes
obs_pp <- MPLS_crime_sf %>%
  group_by(geometry, year, Offense) %>%
  summarise(N = n()) %>%
  ungroup() %>%
  group_by(year, Offense) %>%
  slice_sample(prop = 1) %>% 
  ungroup()

# For observed points, assign offense_lambda = NA
coords <- st_coordinates(obs_pp)
obs_pp$offense_lambda <- NA
obs_pp$x <- coords[, 'X']
obs_pp$y <- coords[, 'Y']

# Integration grid
plt <- ggplot() +
  facet_wrap(~Offense) +
  geom_sf(aes(), data = my_acs, color = 'black', linewidth = 1, fill = 'white', inherit.aes = FALSE) + 
  geom_sf(data = st_buffer(water_sf, dist = 0.05), linewidth = 1.3, fill = 'steelblue', alpha = 1) +
  # observed points
  geom_point(data = filter (obs_pp, year %in% c('2020', '2021', '2022')),
             aes(x = x, y = y),
             color = 'black', size = 2, alpha = 1) +
  
  # the integration grid locations
  geom_point(data = ips_locs,
             aes(x = x, y = y),
             color = 'darkred', size = 0.7, alpha = 0.7) +
  
  # the knot locations
  geom_point(data = knots,
             shape = 23,
             aes(x = x, y = y),
             fill = 'darkgreen', color = 'black', stroke = 0.5,
             size = 5, alpha = 0.9) +
  
  # Police stations
  geom_sf(data=police_pols,
          color = 'white', stroke = 0.7, fill = 'darkblue', inherit.aes = FALSE) +
  
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 30, face = 'bold'),
        plot.subtitle = element_text(hjust = 0.5, size = 36, face = 'bold'),
        strip.text.x = element_text(size = 20, face = 'bold', angle = 0, hjust = 0.5, vjust = 0.5, color = 'white'),
        strip.text.y = element_text(size = 20, face = 'bold', angle = 0, hjust = 1, color = 'white'),
        strip.background = element_rect(fill = 'darkred', color = NA),
        axis.ticks = element_blank(),
        axis.text = element_blank())  +
  labs(title = 'Integration Grid', x = NULL, y = NULL,
       subtitle = bquote(n[A]:n[K] == .(round(n_A, 0))~':'~.(round(n_K, 0)) == .(round(n_A/n_K, 0)) ~' ,'
                         ~ r == .(round(r, 2)) ~' ,'
                         ~ Delta[A] == .(round(A_diff, 2)) ~ '%'))

print(plt)

#---- Hyperparameter Estimation -----
##---- Design Matrices -----
covariates_int <- c('dist_water')
covariates_mark <- c('dist_water','domestic_flag')
intercept <- T
MPLS_crime <- st_drop_geometry(MPLS_crime_sf)

# design matrix for covariates for observed data and integration points
X_mark <- as.matrix(MPLS_crime[,covariates_mark])
X_mark_V <- as.matrix(ips[,covariates_mark])

X_int <- as.matrix(MPLS_crime[,covariates_int])
X_int_V <- as.matrix(ips[,covariates_int])

X_mark <- scale(X_mark, center = T, scale = T)
X_mark_V <- scale(X_mark_V, center = T, scale = T)

ips$X_int_c <- scale(X_int_V, center = T, scale = T)

if(intercept == T) {
  X_int <- cbind(1, X_int)
  X_int_V <- cbind(1, X_int_V)
  X_mark <- cbind(1, X_mark)
  X_mark_V <- cbind(1, X_mark_V)
  
  f_int <- formula(N ~ X_int_c + (1|year) + offset(1*(log(pop_dens/mean(pop_dens)))))
  f_mark <- formula(cbind(Y, N-Y) ~ (1|year) + X_mark)
} else {
  f_int <- formula(N ~ 0 + X_int_c + (1|year) + offset(1*(log(pop_dens/mean(pop_dens)))))
  f_mark <- formula(cbind(Y, N-Y) ~ 0 + (1|year) + X_mark)
}

MPLS_crime$X_int <- X_int
MPLS_crime$X_mark <- X_mark
ips$X_int <- X_int_V
ips$X_mark <- X_mark_V

p_int <- ncol(X_int)
p_mark <- ncol(X_mark)


##---- Frequentist GLMMs -----
# Poisson GLMM to get starting values for beta
h_max <- max(dist(st_coordinates(ips_sf$centroid)))
bin_width <- h_max / 25        
int_mod <- glmer(f_int, data=ips, family = 'poisson',
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(npt = 5)))
ips_sf$pearson_resid_int <- residuals(int_mod, type = "pearson")

# Empirical variogram
vg_emp <- variogram(pearson_resid_int ~ 1, data = ips_sf, width = bin_width, cutoff = 6)
vgm0   <- vgm(psill = max(vg_emp$gamma)-15, 
              model = "Exp", range = median(vg_emp$dist)-2, nugget = 15)
vg_fit <- fit.variogram(vg_emp[-1,], model = vgm0)
plot(vg_emp, vg_fit, main = "Empirical & Fitted Variogram")
vg_fit

r <- vg_fit[vg_fit$model=="Exp","range"] * 3
print(r)
phi_1 <- 3/r
psill_spatial <- vg_fit[vg_fit$model=="Exp","psill"]

sigma1_hat <- sqrt(psill_spatial)

glmm_summ <- summary(int_mod)
glmm_summ

if(intercept == T) {
  if(p_int==2) {
    beta_names <- c("(Intercept)", "X_int_c")
  } else {
    beta_names <- c("(Intercept)", paste0("X_int_c", covariates_int))
  }
} else {
  if(p_int==1) {
    beta_names <- c("X_int_c")
  } else {
    beta_names <- paste0("X_int_c", covariates_int)
  }
}

# Covariance matrix for RW block
params_beta <- paste0('beta[', 1:p_int, ']')
cov_beta <- as.matrix( vcov(glmm_summ)[ beta_names, beta_names ] )
dimnames(cov_beta) <- list(
  params_beta, params_beta
)

# Frequentist estimates
beta_MLE <- glmm_summ$coefficients[,1]
SE_beta_MLE <- glmm_summ$coefficients[,2]
tau1_hat <- as.numeric(sqrt(VarCorr(int_mod)$year))
delta_1_MLE <- ranef(int_mod)$year$`(Intercept)`


# Hyperparameter Empirical Bayes
m_beta   <- mean(beta_MLE)            # should be ~ mu_beta
s2_beta  <- var(beta_MLE)             # sample variance of the MLEs
v_beta   <- mean(SE_beta_MLE^2)           # average sampling variance

nu2_beta <- pmax(s2_beta - v_beta, 1e-8)  # ensure non‐negative
nu_beta   <- sqrt(nu2_beta)               # this is your EB prior SD

# Binomial GLMM to get starting values for gamma
mark_mod <- lme4::glmer(f_mark, data=filter(ips, N>0), family = 'binomial',
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(npt = 5)))
ips_sf$pearson_resid_mark <- NA
ips_sf$pearson_resid_mark[ips_sf$N > 0] <- residuals(mark_mod, type = "pearson")
vg_emp <- variogram(pearson_resid_mark ~ 1, data = filter(ips_sf, N > 0), width = bin_width, cutoff = 5)

vgm0   <- vgm(psill = max(vg_emp$gamma)-1, 
              model = "Exp", range = 1/6, nugget = 1)
vg_fit <- fit.variogram(vg_emp[-1,], model = vgm0)
plot(vg_emp, vg_fit, main = "Empirical & Fitted Variogram")
glmm_summ <- summary(mark_mod)
glmm_summ

if(intercept == T) {
    gamma_names <- c("(Intercept)", paste0("X_mark", covariates_mark))
} else {
  if(p_mark==1) {
    gamma_names <- c("X_mark")
  } else {
    gamma_names <- paste0("X_mark", covariates_mark)
  }
}

# Covariance matrix for RW block
params_gamma <- paste0('gamma[', 1:p_mark, ']')
cov_gamma <- as.matrix( vcov(glmm_summ)[ gamma_names, gamma_names ] )
dimnames(cov_gamma) <- list(
  params_gamma, params_gamma
)

# Frequentist estimates
tau2_hat <- as.numeric(sqrt(VarCorr(mark_mod)$year))
gamma_MLE <- glmm_summ$coefficients[,1]
SE_gamma_MLE <- glmm_summ$coefficients[,2]
delta_2_MLE <- ranef(mark_mod)$year$`(Intercept)`

# Hyperparameter Empirical Bayes
m_gamma   <- mean(gamma_MLE)            # should be ~ mu_gamma
s2_gamma  <- var(gamma_MLE)             # sample variance of the MLEs
v_gamma   <- mean(SE_gamma_MLE^2)           # average sampling variance

nu2_gamma <- pmax(s2_gamma - v_gamma, 1e-8)  # ensure non‐negative
nu_gamma   <- sqrt(nu2_gamma)               # EB prior SD


##----  Hyperpamaters nu -----
sigma_hat <- c(sigma1_hat)
v_sigma <- 0.05 * sigma_hat^2  # assume 5% relative error
s2_sigma <- sigma_hat^2

var_log_sigma <- v_sigma / s2_sigma
nu_log_sigma <- sqrt(pmax(var_log_sigma, 1e-8))

tau_hat <- c(tau1_hat, tau2_hat)
s2_tau <- c(var(delta_1_MLE), var(delta_2_MLE))
v_tau  <- 0.05 * s2_tau

var_log_tau <- v_tau / s2_tau
nu_log_tau <- sqrt(pmax(var_log_tau, 1e-8))

#-----NIMBLE Data Structures-----
## ----- GP Covariance -----
phi_1 <- phi_1 

# Observed
MPLS_crime_sf$X_int <- MPLS_crime$X_int
MPLS_crime_sf$X_mark <- MPLS_crime$X_mark
obs_coords <- MPLS_crime_sf %>%
  mutate(x = st_coordinates(geometry)[,'X'],
         y = st_coordinates(geometry)[,'Y'],
         txt = st_as_text(geometry)) %>%
  distinct(txt, .keep_all = TRUE) %>%
  mutate(loc = row_number()) %>%
  select(x, y, pop, pop_dens, X_int, geometry, txt)

loc_idx <- match(st_as_text(MPLS_crime_sf$geometry), obs_coords$txt)


ips_loc <- ips %>% 
  distinct(x, y, X_int, pop, pop_dens) %>% 
  mutate(loc_A = row_number())

loc_V_idx  <- match(
  paste0(ips$x, '-', ips$y),
  paste0(ips_loc$x, '-', ips_loc$y)
)

if(intercept == T) {
  obs_coords$X_int <- cbind(1, scale(obs_coords$X_int[,-1], scale = T))
  ips_loc$X_int <- cbind(1, scale(ips_loc$X_int[,-1], scale = T))
} else {
  obs_coords$X_int <- scale(obs_coords$X_int, scale = T)
  ips_loc$X_int <- scale(ips_loc$X_int, scale = T)
}

# Distance calculations
dists_obs_knots <- as.matrix(proxy::dist(cbind(obs_coords$x, obs_coords$y), knots))
dists_int_knots <- as.matrix(proxy::dist(cbind(ips_loc$x, ips_loc$y), knots))

kern <- 'Exp'
if(kern == 'Gaus') {
  # pre-compute matrices for faster computation since phi is fixed
  W1star_gauss_cov <- gauss_corr(dists_knots, phi = phi_1)
  inv_W1star_gauss_cov <- solve(W1star_gauss_cov)
  C_cross_obs_gauss_cov_1 <- gauss_corr(dists_obs_knots, phi= phi_1)
  C_cross_int_gauss_cov_1 <- gauss_corr(dists_int_knots, phi= phi_1)

  # sigma cancels out of this part since one is inverse
  C_cross_C_inv_obs_1 <- C_cross_obs_gauss_cov_1 %*% inv_W1star_gauss_cov
  C_cross_C_inv_int_1 <- C_cross_int_gauss_cov_1 %*% inv_W1star_gauss_cov

  
  L1 <- chol(W1star_gauss_cov)
} else if (kern == 'Exp') {
  # pre-compute matrices for faster computation since phi is fixed
  W1star_exp_cov <- exp_corr(dists_knots, phi = phi_1)
  inv_W1star_exp_cov <- solve(W1star_exp_cov)
  C_cross_obs_exp_cov_1 <- exp_corr(dists_obs_knots, phi= phi_1)
  C_cross_int_exp_cov_1 <- exp_corr(dists_int_knots, phi= phi_1)

  # sigma cancels out of this part since one is inverse
  C_cross_C_inv_obs_1 <- C_cross_obs_exp_cov_1 %*% inv_W1star_exp_cov
  C_cross_C_inv_int_1 <- C_cross_int_exp_cov_1 %*% inv_W1star_exp_cov
  
  L1 <- chol(W1star_exp_cov)
}

# Check problematic knnots based on covariates
GP_check <- check_gp_covariate_collinearity(C_cross_C_inv_obs_1, obs_coords$X_int, 
                                            knots, 
                                            window_sf, 
                                            threshold = 0.4, 
                                            title = "GP Basis: Problematic Knots Highlighted")

year_idx <- as.integer(as.factor(MPLS_crime_sf$year))
n_years <- length(unique(year_idx)) 
loc_idx <- as.integer(as.factor(st_as_text(MPLS_crime_sf$geometry)))



##-----Data-----
data_list <- list(
  s = as.numeric(loc_idx),
  Y       = MPLS_crime_sf$Y,
  area = ips$clipped_area,
  log_pop_dens     =  as.numeric(1*(log(obs_coords$pop_dens))) - mean(log(obs_coords$pop_dens)),
  log_pop_dens_V   = as.numeric(1*(log(ips_loc$pop_dens))) - mean(log(ips_loc$pop_dens)),
  X_int   = obs_coords$X_int,  # n_locations × p_int,
  X_int_V = ips_loc$X_int,
  X_mark  = MPLS_crime_sf$X_mark,
  X_mark_V = ips$X_mark
)

##-----Constants-----
constants_list <- list(
  loc         = loc_idx,
  loc_V       = loc_V_idx,
  year = year_idx,
  year_V = as.integer(as.factor(ips$year)),
  n_points    = nrow(MPLS_crime_sf),
  n_locations = nrow(obs_coords),
  n_A         = n_A,
  n_V         = n_V,
  n_K         = n_K,
  n_years     = n_years,
  L1 = L1,
  C_cross_C_inv_obs_1 = C_cross_C_inv_obs_1,
  C_cross_C_inv_int_1 = C_cross_C_inv_int_1,
  p_int       = p_int,
  p_mark      = p_mark,
  mu_beta     = as.numeric(beta_MLE),
  mu_gamma    = rep(0, p_mark),
  nu_beta     = rep(nu_beta, p_int), 
  nu_gamma    = rep(nu_gamma, p_mark),
  mu_sigma = as.numeric(log(1.5)),
  mu_tau = as.numeric(log(tau_hat)),
  nu_sigma = as.numeric(0.3),
  nu_tau = as.numeric(nu_log_tau)
)

#----- Sampler -----
##----- Test Run -----
# Test run (to get RW_block covariance matrices & to make your life easier :] )
system.time({
  testrun <- test_run(model_code, data_list, constants_list, cov_beta, cov_gamma,
                      n_burnin = 6000, n_iter = 8000, n_chains=1)
})

# Grab proposal covariances matrices for RW_block
samples_comb <- do.call(rbind, testrun$samples)
cov_GP <- testrun$cov_GP
cov_beta <- testrun$cov_beta
cov_gamma <- testrun$cov_gamma

# Check for problem knots based on test run
knot_problem <- visualize_problematic_knots(samples_comb = samples_comb, knots = knots, 
                            window_sf = window_sf,
                            threshold = 0.8)

# Create parallel cluster
cl <- makeCluster(nchains)

# Export necessVry objects to each worker
clusterExport(cl, c('model_code', 'data_list', 'constants_list', 'niter', 'n_years',
                    'n_K', 'p_int', 'p_mark', 
                    'nburnin', 'thin', 'knots', 'beta_MLE', 'cov_GP', 'cov_beta', 'cov_gamma', 'delta_2_MLE', 'delta_2_MLE',
                    'gamma_MLE'))

# Load nimble on workers
clusterEvalQ(cl, {
  library(nimble)
  library(beepr)
  source('utils.R')
  nimbleOptions(verbose = F)
})

#----- Model Run -----
  samples_list <- pbapply::pblapply(1:nchains,function(chain_id) {
    set.seed(chain_id)
    
    inits <- inits_fn()
    
    # Save to file1star
    saveRDS(inits, file = paste0('inits_chain_', chain_id, '.rds'))
    
    if(!exists('d_lgcp'))     stop('d_lgcp is missing on this worker!')
    
  
    # Build model 
    model <- nimbleModel(model_code,
                         data = data_list,
                         constants = constants_list,
                         inits = inits)
    
    config <- configureMCMC(model)
    
    # Remove defaults
    config$removeSampler(c('z1')) 
    config$removeSampler(c('u1', 'u2')) 
    config$removeSampler('beta')
    config$removeSampler('gamma')
    config$removeSampler('log_sigma')
    config$removeSampler('log_tau')
    
    # Add optimized sVmplers
    # Block temporal RE
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
    
    # Block GP bases
    config$addSampler(target = c('log_sigma', 'z1'),
                      type   = 'RW_block',
                      control = list(propCov = cov_GP,
                                     adaptive      = TRUE,
                                     adaptInterval =  500,
                                     scale         = 0.01
                      ))
    
    # Block beta
    config$addSampler(
      target = c('beta'),
      type   = 'RW_block',
      control = list(propCov  = cov_beta,
                     adaptive      = TRUE,
                     adaptInterval =  500,
                     scale         = 0.05
        )
      )
    
    # Block gamma
    config$addSampler(
      target = c('gamma'), 
      type = 'RW_block', 
      control = list(propCov       = cov_gamma, 
                     adaptive      = TRUE,
                     adaptInterval = 500,
                     scale         = 0.05
                     )
      )
  
    # Block tau
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
    
    # Running one chain
    runMCMC(compiled$mcmc,
              nburnin = nburnin,
              niter = niter,
              thin = thin,
              nchains = 1,
              samplesAsCodaMCMC = TRUE,
              progressBar = T)
    
  }, cl = cl)

stopCluster(cl)


#----- MCMC Summaries -----
# GP summary
model_summary <- tibble(
  metric = c('Knots', 'Integration Points', 'GP Effective Range', 'Burnin', 'Total Iterations', 
             'Thinning'#, 'Multivariate PSRF'
  ),
  value = c(n_K, nrow(ips_locs), 3/phi_1, nburnin, niter, thin#, gelman_rubin$mpsrf
  ),
  type = c('int', 'int', 'num', 'int', 'int', 'int'#, 'num'
  )
)
print(model_summary[,-3])

# Combine into mcmc.list for diagnostics
samples <- coda::mcmc.list(samples_list)
rm(samples_list)

# Parameter names
params_beta <- paste0('beta[', 1:p_int, ']')
params_gamma <- paste0('gamma[', 1:p_mark, ']')
params_delta <- c(paste0( 'delta_1[', 1:n_years, ']'), 
                  paste0('delta_2[', 1:n_years, ']')
                  )
params_varcomp <- c('sigma', # 'sigma[2]', 
                    'tau[1]', 'tau[2]')

##----- Convergence -----
# Check Gelman-Rubin
gelman_rubin <- coda::gelman.diag(samples[,c(params_beta, params_gamma,
                                             params_delta, params_varcomp)], 
                                  multivariate = F)
print(gelman_rubin)

# Plot traceplots
plot(samples[,params_varcomp], density = T)
plot(samples[,c(params_beta, params_gamma)], density = F) 
plot(samples[,params_delta], density = F) 

#----- Posterior Summaries -----
##----- Posterior Inference -----
summary(samples[,c(params_beta, params_gamma,
                   params_delta, 
                   params_varcomp
)])

# Combine the samples
samples_comb <- do.call(rbind, samples)

##----- Posterior Intensity Maps -----
posterior_intensity(samples_comb = samples_comb, model_summary = model_summary)

ggsave('plots/posterior_mean.png', width = 12, height = 10, dpi = 300)


