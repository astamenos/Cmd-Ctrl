# ----- Crime -----
# Loading the Minneapolis crime data
url <- "https://services.arcgis.com/afSMGVsC7QlRK1kZ/arcgis/rest/services/Crime_Data/FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson"
crime_data <- st_read(url)

# ----- Landscapes -----
##----- Police Stations -----
# Get the bounding box for Minneapolis, MN
bbox_minneapolis <- getbb("Minneapolis, MN")

# Build an Overpass query for all OSM features tagged amenity=police
query_police <- opq(bbox_minneapolis) %>%
  add_osm_feature(key = "amenity", value = "police")

# Download the results as sf
police_osm <- osmdata_sf(query_police)

# Extract point‐mapped stations, plus centroids of any polygons
police_pols <- police_osm$osm_polygons %>%
  st_transform(crs = crs_use) %>%
  mutate(geometry = geometry/1000)

# Get centroids for polygonal station footprints
police_cents <- st_centroid(police_pols$geometry)
police_pols <- police_pols %>%
  st_sf() %>% 
  select(osm_id, name) %>%
  st_set_geometry(police_cents) %>%
  st_buffer(
    dist         = 0.2,
    joinStyle    = "MITRE",
    endCapStyle  = "SQUARE"
  ) %>%
  st_set_crs(crs_use)


##----- MPLS -----
# Loading Minneapolis neighborhood shapefile
MPLS <- st_read('data/Minneapolis_Neighborhoods/Minneapolis_Neighborhoods.shp')


# 
MPLS <- MPLS %>%
  st_transform(crs = crs_use) %>%
  mutate(geometry = geometry/1000)
city_poly <- st_union(MPLS)
st_crs(city_poly) <- crs_use

#xmin <- min(st_drop_geoometry(city_poly[,1]))
#ymin <- min(city_poly[,2])

#----- ACS -----
my_acs <- get_acs(
  geography = "tract",
  variables = c(pop = "B01003_001", pov = "B17001_002"),
  output    = 'wide',
  state     = "MN",
  county    = "Hennepin", 
  year      = 2019,
  geometry  = TRUE
) %>%
  st_transform(crs = crs_use) %>%
  mutate(
    area_km2 = as.numeric(st_area(geometry)) / 1e6,
    pop_dens = popE / area_km2,
    pct_pov  = 100*povE/popE
  ) %>%
  mutate(geometry = geometry/1000) %>%
  st_set_crs(crs_use) %>%  
  st_intersection(city_poly, .predicate = st_within) 

##----- Water -----
water_sf <- tigris::area_water(state = '27', county = '053') %>%
  st_transform(crs = crs_use) %>%
  mutate(geometry = geometry/1000) %>%
  st_set_crs(crs_use) %>%
  st_intersection(city_poly, .predicate = st_within) %>%
  mutate(cts = st_centroid(geometry)) %>%
  st_union()    

##----- Roads -----
roads_sf <- tigris::primary_secondary_roads(state = '27') %>%
  st_transform(crs = crs_use) %>%
  mutate(geometry = geometry/1000) %>%
  st_set_crs(crs_use) %>%
  st_intersection(city_poly, .predicate = st_within) %>%
  mutate(cts = st_centroid(geometry)) %>%
  st_union()

mainroads_sf <- tigris::primary_roads() %>%
  st_transform(crs = crs_use) %>%
  mutate(geometry = geometry/1000) %>%
  st_set_crs(crs_use) %>%
  st_intersection(city_poly, .predicate = st_within) %>%
  mutate(cts = st_centroid(geometry)) %>%
  st_union()

# ----- Crime -----
crime_data <- crime_data %>%
  mutate(
    occurred_datetime = as.POSIXct(Occurred_Date / 1000, origin = "1970-01-01", tz = "UTC"),
    occurred_date = format(occurred_datetime, "%Y-%m-%d"),
    period = paste0(format(occurred_datetime, "%Y"), "Q",
                  ceiling(as.integer(format(occurred_datetime, "%m")) / 3))) %>%
  filter(period %in% periods)

# The kind of crimes people commit...
offense_cat <- st_drop_geometry(crime_data) %>%
  select(Offense_Category, Offense) %>%
  distinct() %>%
  arrange(Offense_Category, Offense)

#
MPLS_crime_sf <- crime_data %>%
  drop_na(Problem_Initial | Problem_Final) %>%
  mutate(
    period = factor(period, levels = levels(periods_fac)),
    call_type = trimws(gsub("~PD", "", gsub("\\s*\\([^\\)]+\\)", "", Problem_Initial))),
    Offense_Category = factor(Offense_Category),
    Offense = if_else(grepl('Domestic', Offense), 'Domestic', 'Non-Domestic'),
    Offense = factor(Offense),
    Crime_Count_cat = case_when(
      Crime_Count == 1 ~ '1',
      Crime_Count > 1 & Crime_Count <= 4 ~ '1 - 4',
      Crime_Count > 4 ~ '4 +'
    )
  ) %>%
  filter(
    Offense_Category %in% offense_cat$Offense_Category[grep('Aggravated Assault', offense_cat$Offense)],
    period %in% periods,
    Latitude != 0,
    Longitude != 0,
  ) %>%
  drop_na(Precinct, Ward) %>%
  select(Precinct, Ward, geometry, occurred_datetime, occurred_date, period, 
         Case_Number, Type, Offense_Category, Offense,
         Crime_Count, Crime_Count_cat, Problem_Initial, Problem_Final, call_type) %>%
  distinct() %>%
  st_transform(crs = crs_use) %>%
  mutate(geometry = geometry/1000,
         geometry_str = st_as_text(geometry)) %>%
  st_set_crs(crs_use)

# Deduplicate based on datetime and geometry
dup_datetime_geom <- MPLS_crime_sf %>%
  st_drop_geometry() %>%
  group_by(occurred_date, geometry_str) %>%
  filter(n_distinct(Offense) > 1) %>%
  ungroup() %>%
  distinct(occurred_datetime, occurred_date, geometry_str)

# Keep only Domestic where duplicated, and all others as-is
MPLS_crime_sf <- MPLS_crime_sf %>%
  anti_join(dup_datetime_geom, by = c("occurred_date", "geometry_str")) %>%
  bind_rows(
    MPLS_crime_sf %>%
      semi_join(dup_datetime_geom, by = c("occurred_date", "geometry_str")) %>%
      filter(Offense == "Domestic")
  ) %>%
  arrange(Precinct, Ward, geometry_str, period)


# Derive domestic_flag based on relevant call types
domestic_pattern <- str_c(domestic_terms, collapse = "|")

MPLS_crime_sf <- MPLS_crime_sf %>%
  rename(N=Crime_Count) %>%
  mutate(Y = if_else(Offense == 'Domestic', 1, 0),
         domestic_flag = if_else(
           str_detect(call_type, regex(domestic_pattern, ignore_case = TRUE)),
           1L, 0L
         ),
         other = 1 - domestic_flag,
         across(c(other, domestic_flag),
                ~ if_else(is.na(.x), median(.x, na.rm=TRUE), .x))
  ) %>%
  filter(N==1)

# Jitter points slightly
duplicated_pts <- MPLS_crime_sf %>%
  group_by(Precinct, Ward, geometry, period) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(Precinct, Ward, geometry, period) %>%
  pivot_wider(names_from = period, values_from = n, values_fill = 0) %>%
  filter(if_any(all_of(periods), ~ . > 1)) %>%
  mutate(geometry_str = st_as_text(geometry)) %>%
  select(geometry, geometry_str) %>%
  mutate(coords = st_coordinates(geometry),
         coords_j = jitter(coords, d_jitt),
         x = coords[,'X'],
         y = coords[,'Y'],
         x_j = coords_j[,'X'],
         y_j = coords_j[,'Y']) %>%
  select(-coords, -coords_j) %>%
  st_drop_geometry()

# Only jitter non-unique locations
MPLS_crime_sf <- MPLS_crime_sf %>%
  mutate(coords = st_coordinates(MPLS_crime_sf),
         x = coords[,'X'],
         y = coords[,'Y']) %>%
  st_drop_geometry() %>%
  left_join(duplicated_pts) %>%
  mutate(x = if_else(is.na(x_j), x, x_j),
         y = if_else(is.na(y_j), y, y_j)) %>%
  st_as_sf(coords = c('x', 'y'), crs = st_crs(MPLS_crime_sf)) %>%
  select(-x_j, -y_j) %>%
  filter(N==1) %>%
  st_set_crs(crs_use)

# Yearly data summary of assaults
crime_summ <- st_drop_geometry(MPLS_crime_sf) %>% 
  group_by(period, Offense) %>% 
  summarise(N = sum(N),
            across(c(domestic_flag, other), ~mean(.x, na.rm = TRUE)))

#----- Window -----
# Create window using shapefile boundaries
window <- as.owin(city_poly)
window_sf <- st_as_sf(window)
st_crs(window_sf) <- st_crs(my_acs)

# Remove bodies of water from window
window_sf <- st_difference(window_sf, water_sf)
window <- as.owin(window_sf)

# Remove points outside of window
MPLS_crime_sf <- MPLS_crime_sf %>%
  st_filter(window_sf) 
MPLS_crime_sf <- st_join(
    MPLS_crime_sf,
    my_acs,
    join = st_within,
    left = TRUE
  ) %>%
  drop_na(pop_dens)

# Remove geometry to create CSV
MPLS_crime <- MPLS_crime_sf %>%
  mutate(x = coords[,'X'],
         y = coords[,'Y']) %>%
  st_drop_geometry() %>%
  select(-coords) 

# Write to CSV
write_csv(MPLS_crime, 'data/MPLS_crime.csv')
