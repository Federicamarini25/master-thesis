rm(list = ls(all = TRUE))     
cat("\014")                    

# Load necessary libraries
library(httr)
library(jsonlite)
library(ggplot2)
library(caret)  
library(dplyr)  
library(mgcv)
library(tseries)
library(forecast)
library(imputeTS)
library(zoo)  
library(tidyr)
library(future.apply)
library(animation)
library(furrr)
library(geodata)
library(terra)
library(splines)
library(cglasso)
library(stats)
library(gridExtra)
library(grid)
library(ggpubr)

# ---- Functions ----
# Fetch locations data
fetch_locations_data <- function(domain) {
  base_url <- "http://www.oasi.ti.ch/web/rest/locations"
  params <- list(domain = domain)
  response <- GET(base_url, query = params)
  
  if (response$status_code == 200) {
    data <- fromJSON(content(response, as = "text", encoding = "UTF-8"), flatten = TRUE)
    if (!"name" %in% names(data) || !"code" %in% names(data)) {
      warning("Invalid data structure received from locations API.")
      return(NULL)
    }
    return(data)
  } else {
    warning("Failed to get locations data. Status code:", response$status_code)
    return(NULL)
  }
}

# Fetch time series data
fetch_time_series_data <- function(domain, location_code, parameter, resolution, from_date, to_date) {
  base_url <- "http://www.oasi.ti.ch/web/rest/measure/csv"
  params <- list(
    domain = domain,
    location = location_code,
    parameter = parameter,
    resolution = resolution,
    from = format(as.Date(from_date), "%Y-%m-%d"),
    to = format(as.Date(to_date), "%Y-%m-%d")
  )
  
  response <- GET(base_url, query = params)
  
  if (response$status_code == 200) {
    return(content(response, as = "text", encoding = "UTF-8"))
  } else {
    warning("Failed to get data. Status code:", response$status_code)
    return(NULL)
  }
}

# Process and append data
process_and_append_data <- function(response_data, parameter, location_name, data_frame) {
  data_lines <- strsplit(response_data, "\n")[[1]]
  data_lines <- data_lines[!grepl("^#", data_lines)]
  
  if (length(data_lines) > 0) {
    data_clean <- paste(data_lines, collapse = "\n")
    data_df <- read.csv(text = data_clean, sep = ";", header = TRUE)
    
    if ("data" %in% names(data_df) && parameter %in% names(data_df)) {
      data_df$data <- as.POSIXct(data_df$data, format="%d.%m.%Y %H:%M", tz = "UTC")
      if (all(is.na(data_df$data))) {
        warning("Date conversion failed. Check the date format.")
      }
      data_df[[parameter]] <- as.numeric(data_df[[parameter]])
      data_df$Location <- location_name
      data_frame <- rbind(data_frame, data_df)
    }
  }
  return(data_frame)
}

# Preprocess data (convert negative values to NA)
preprocess_data <- function(data_frame) {
  numeric_cols <- sapply(data_frame, is.numeric)
  data_frame[numeric_cols] <- lapply(data_frame[numeric_cols], function(x) ifelse(x < 0, NA, x))
  return(data_frame)
}

# Setting working directory
working_dir = "/Users/federicamarini/Desktop/Master Thesis"            
setwd(working_dir)     

# --- Define Parameters ---
domain <- "surfacewater"
parameter <- "Q"
from_date <- as.Date("2020-01-01")
to_date <- as.Date("2023-12-31")

# Parameters for precipitations 
domain_meteo <- "meteo"
parameter_prec <- "Prec"

# Fetch location data to obtain coordinates
locations_df_q <- fetch_locations_data(domain)
if (is.null(locations_df_q) || !"name" %in% names(locations_df_q) || !"coordinates.x" %in% names(locations_df_q) || !"coordinates.y" %in% names(locations_df_q)) {
  stop("Failed to fetch location data for river flow.")
}

location_names <- as.vector(locations_df_q$name)

# Initialize an empty dataframe for storing all locations' data
all_q_data <- data.frame()

# Loop through each location to fetch and process data
for (location_name in location_names) {
  
  location_row <- locations_df_q[locations_df_q$name == location_name, ]
  
  # Skip if location not found
  if (nrow(location_row) == 0) next
  
  location_code <- location_row$code
  coord_x <- location_row$coordinates.x
  coord_y <- location_row$coordinates.y
  
  # Fetch time series data
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  
  # Skip if data fetch fails
  if (is.null(q_data_response)) next
  
  # Process data
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  
  # Skip if no data available
  if (nrow(data) == 0) next
  
  # Preprocess data
  data <- preprocess_data(data)
  
  # Keep only relevant columns and add coordinate & location information
  data <- data %>%
    select(Q, data) %>%
    mutate(
      location_name = location_name,
      coordinates.x = coord_x,
      coordinates.y = coord_y
    )
  
  # Append data to final dataset
  all_q_data <- rbind(all_q_data, data)
}


# --- Same Process for Precipitation Data ---

# Fetch location data for precipitation
locations_df_prec <- fetch_locations_data(domain_meteo)
if (is.null(locations_df_prec) || !"name" %in% names(locations_df_prec) || !"coordinates.x" %in% names(locations_df_prec) || !"coordinates.y" %in% names(locations_df_prec)) {
  stop("Failed to fetch location data for precipitation.")
}


# If row 43 is incorrect, filter by condition rather than row index
locations_df_prec <- locations_df_prec %>% filter(code != "auto_57") 

prec_locations <- as.vector(locations_df_prec$name)

# Initialize an empty dataframe for storing all locations' precipitation data
all_prec_data <- data.frame()

# Loop through each location to fetch and process precipitation data
for (prec_location in prec_locations) {
  
  location_row <- locations_df_prec[locations_df_prec$name == prec_location, ]
  
  # Skip if location not found
  if (nrow(location_row) == 0) next
  
  location_code <- location_row$code
  coord_x <- location_row$coordinates.x
  coord_y <- location_row$coordinates.y
  
  # Fetch time series data for precipitation
  prec_data_response <- fetch_time_series_data(domain_meteo, location_code, parameter_prec, "d", from_date, to_date)
  
  # Skip if data fetch fails
  if (is.null(prec_data_response)) next
  
  # Process data
  data <- process_and_append_data(prec_data_response, parameter_prec, prec_location, data.frame())
  
  # Skip if no data available
  if (nrow(data) == 0) next
  
  # Preprocess data
  data <- preprocess_data(data)
  
  # Keep only relevant columns and add coordinate & location information
  data <- data %>%
    select(prec = !!parameter_prec, data) %>%
    mutate(
      location_name = prec_location,
      coordinates.x = coord_x,
      coordinates.y = coord_y
    )
  
  # Append data to final dataset
  all_prec_data <- rbind(all_prec_data, data)
}

all_q_data[, "data"] <- as.Date(all_q_data[, "data"])
all_prec_data[, "data"] <- as.Date(all_prec_data[, "data"])

# Imputation ----
library(dplyr)
library(tidyr)
library(imputeTS)

# Delete the outliers
all_prec_data$prec[all_prec_data$prec > 500] <- NA
all_q_data$Q[all_q_data$Q > 1000] <- NA

# Impute all_prec_data
all_prec_data_imputed <- all_prec_data %>%
  group_by(location_name) %>%
  arrange(as.Date(data)) %>%
  mutate(
    prec = zoo::na.approx(prec, x = as.Date(data), na.rm = FALSE),
    prec = zoo::na.locf(prec, na.rm = FALSE),
    prec = zoo::na.locf(prec, fromLast = TRUE, na.rm = FALSE)
  ) %>%
  ungroup()

# Impute all_q_data
all_q_data_imputed <- all_q_data %>%
  group_by(location_name) %>%
  arrange(as.Date(data)) %>%
  mutate(
    Q = zoo::na.approx(Q, x = as.Date(data), na.rm = FALSE),
    Q = zoo::na.locf(Q, na.rm = FALSE),
    Q = zoo::na.locf(Q, fromLast = TRUE, na.rm = FALSE)
  ) %>%
  ungroup()

# Check for remaining NAs
sum(is.na(all_q_data_imputed$Q))  # 68667
sum(is.na(all_prec_data_imputed$prec))  # 2922

# --- Remove locations with all NA values in all_q_data ---
all_q_data_imputed <- all_q_data_imputed %>%
  group_by(location_name) %>%
  filter(!all(is.na(Q))) %>%  # Remove locations with only missing values
  ungroup()

# --- Remove locations with all NA values in all_prec_data ---
all_prec_data_imputed <- all_prec_data_imputed %>%
  group_by(location_name) %>%
  filter(!all(is.na(prec))) %>%  # Remove locations with only missing values
  ungroup()

# --- Impute remaining NAs in all_prec_data with location-specific mean ---
all_prec_data_imputed <- all_prec_data_imputed %>%
  group_by(location_name) %>%
  mutate(prec = ifelse(is.na(prec), mean(prec, na.rm = TRUE), prec)) %>%
  ungroup()

# Verify that no missing values remain
colSums(is.na(all_prec_data_imputed))
colSums(is.na(all_q_data_imputed))

# Center the observations of river flow 
all_q_data_imputed <- all_q_data_imputed %>%
  mutate(logQ = log1p(Q)) %>%
  group_by(location_name) %>%
  mutate(Q_centered = logQ - mean(logQ, na.rm = TRUE)) %>%
  ungroup()

# ADD ELEVATION TO THE DATASETS 
library(geodata)
library(terra)
library(dplyr)

# STEP 1: Download both SRTM tiles for Switzerland and Italy
elev_che <- elevation_30s(country = "CHE", path = tempdir())
elev_ita <- elevation_30s(country = "ITA", path = tempdir())

# STEP 2: Reproject both to LV95 (EPSG:2056)
elev_che_lv95 <- project(elev_che, "EPSG:2056")
elev_ita_lv95 <- project(elev_ita, "EPSG:2056")

# STEP: Resample Italian raster to match Swiss raster (same resolution + extent grid)
elev_ita_lv95_resampled <- resample(elev_ita_lv95, elev_che_lv95)

# Now mosaic works!
elev_combined_lv95 <- mosaic(elev_che_lv95, elev_ita_lv95_resampled)


# STEP 4: Create expanded bounding box from your data
coords_mat <- as.matrix(all_prec_data_imputed[, c("coordinates.x", "coordinates.y")])
pts_vect <- vect(coords_mat, type = "points", crs = "EPSG:2056")

bbox <- ext(pts_vect)
bbox_expanded <- ext(bbox$xmin - 20000, bbox$xmax + 20000,
                     bbox$ymin - 20000, bbox$ymax + 20000)

# STEP 5: Crop to expanded bounding box (includes parts of Italy)
elev_expanded <- crop(elev_combined_lv95, bbox_expanded)

# STEP 7: Extract elevation at precipitation points
prec_points <- vect(all_prec_data_imputed, geom = c("coordinates.x", "coordinates.y"), crs = "EPSG:2056")
prec_elev <- extract(elev_expanded, prec_points)

# STEP 8: Merge into original data
all_prec_data_with_elev <- bind_cols(all_prec_data_imputed, elevation = prec_elev[, 2])

# same for river flow ----
# STEP 1: Create bounding box around river flow station coordinates ----
coords_q <- as.matrix(all_q_data_imputed[, c("coordinates.x", "coordinates.y")])
q_pts_vect <- vect(coords_q, type = "points", crs = "EPSG:2056")

# STEP 2: Expand bounding box by 20km (same as for precipitation) ----
bbox_q <- ext(q_pts_vect)
bbox_q_expanded <- ext(bbox_q$xmin - 20000, bbox_q$xmax + 20000,
                       bbox_q$ymin - 20000, bbox_q$ymax + 20000)

# STEP 3: Crop elevation raster to this expanded extent ----
elev_expanded_q <- crop(elev_combined_lv95, bbox_q_expanded)

# STEP 5: Convert Q station dataframe to SpatVector
q_points <- vect(all_q_data_imputed,
                 geom = c("coordinates.x", "coordinates.y"),
                 crs = "EPSG:2056")

# STEP 6: Extract elevation values at each river flow location ----
q_elev <- extract(elev_expanded_q, q_points)

# STEP 7: Merge elevation into the river flow dataset ----
all_q_data_with_elev <- bind_cols(all_q_data_imputed, elevation = q_elev[, 2])

# STEP 8: Inspect
head(all_q_data_with_elev)

# Extract ranges from the datasets
lat_range <- range(all_q_data_imputed$coordinates.y)
long_range <- range(all_q_data_imputed$coordinates.x)

# Create a grid of points ----
dim.grid <- 25
a <- dim.grid 
b <- dim.grid 
N_T <- a * b

new.data <- matrix(NA, a * b, 2, dimnames = list("obs" = paste0("obs", seq(1, a*b)), "coord" = c("coordinates.y", "coordinates.x") ))
new.data[, 1] <- sort(rep(seq(lat_range[1], lat_range[2], length.out = a), b))
new.data[, 2] <- rep(seq(long_range[1], long_range[2], length.out = b), a)
new.data <- as.data.frame(new.data)

# Renaming columns for consistency
colnames(new.data) <- c("coordinates.y", "coordinates.x")
N_T <- dim(new.data)[1]

# NEW YEAR 2024----
# Fetch new data
domain <- "surfacewater"
parameter <- "Q"
from_date_new <- as.Date("2024-01-01")
to_date_new <- as.Date("2024-12-31")

# Parameters for precipitations 
domain_meteo <- "meteo"
parameter_prec <- "Prec"

# Fetch location data to obtain coordinates
locations_df_q <- fetch_locations_data(domain)
if (is.null(locations_df_q) || !"name" %in% names(locations_df_q) || !"coordinates.x" %in% names(locations_df_q) || !"coordinates.y" %in% names(locations_df_q)) {
  stop("Failed to fetch location data for river flow.")
}

location_names <- as.vector(locations_df_q$name)

# Initialize an empty dataframe for storing all locations' data
all_q_data_new <- data.frame()

# Loop through each location to fetch and process data
for (location_name in location_names) {
  
  location_row <- locations_df_q[locations_df_q$name == location_name, ]
  
  # Skip if location not found
  if (nrow(location_row) == 0) next
  
  location_code <- location_row$code
  coord_x <- location_row$coordinates.x
  coord_y <- location_row$coordinates.y
  
  # Fetch time series data
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date_new, to_date_new)
  
  # Skip if data fetch fails
  if (is.null(q_data_response)) next
  
  # Process data
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  
  # Skip if no data available
  if (nrow(data) == 0) next
  
  # Preprocess data
  data <- preprocess_data(data)
  
  # Keep only relevant columns and add coordinate & location information
  data <- data %>%
    select(Q, data) %>%
    mutate(
      location_name = location_name,
      coordinates.x = coord_x,
      coordinates.y = coord_y
    )
  
  # Append data to final dataset
  all_q_data_new <- rbind(all_q_data_new, data)
}

# --- Same Process for Precipitation Data ---

# Fetch location data for precipitation
locations_df_prec <- fetch_locations_data(domain_meteo)
if (is.null(locations_df_prec) || !"name" %in% names(locations_df_prec) || !"coordinates.x" %in% names(locations_df_prec) || !"coordinates.y" %in% names(locations_df_prec)) {
  stop("Failed to fetch location data for precipitation.")
}


# If row 43 is incorrect, filter by condition rather than row index
locations_df_prec <- locations_df_prec %>% filter(code != "auto_57") 

prec_locations <- as.vector(locations_df_prec$name)

# Initialize an empty dataframe for storing all locations' precipitation data
all_prec_data_new <- data.frame()

# Loop through each location to fetch and process precipitation data
for (prec_location in prec_locations) {
  
  location_row <- locations_df_prec[locations_df_prec$name == prec_location, ]
  
  # Skip if location not found
  if (nrow(location_row) == 0) next
  
  location_code <- location_row$code
  coord_x <- location_row$coordinates.x
  coord_y <- location_row$coordinates.y
  
  # Fetch time series data for precipitation
  prec_data_response <- fetch_time_series_data(domain_meteo, location_code, parameter_prec, "d", from_date_new, to_date_new)
  
  # Skip if data fetch fails
  if (is.null(prec_data_response)) next
  
  # Process data
  data <- process_and_append_data(prec_data_response, parameter_prec, prec_location, data.frame())
  
  # Skip if no data available
  if (nrow(data) == 0) next
  
  # Preprocess data
  data <- preprocess_data(data)
  
  # Keep only relevant columns and add coordinate & location information
  data <- data %>%
    select(prec = !!parameter_prec, data) %>%
    mutate(
      location_name = prec_location,
      coordinates.x = coord_x,
      coordinates.y = coord_y
    )
  
  # Append data to final dataset
  all_prec_data_new <- rbind(all_prec_data_new, data)
}

all_q_data_new[, "data"] <- as.Date(all_q_data_new[, "data"])
all_prec_data_new[, "data"] <- as.Date(all_prec_data_new[, "data"])

# Imputation ----
library(dplyr)
library(tidyr)
library(imputeTS)

# Delete the outliers
all_prec_data_new$prec[all_prec_data_new$prec > 500] <- NA
all_q_data_new$Q[all_q_data_new$Q > 1000] <- NA

# Impute all_prec_data_new
all_prec_data_new_imputed <- all_prec_data_new %>%
  group_by(location_name) %>%
  arrange(as.Date(data)) %>%
  mutate(
    prec = na.approx(prec, x = as.Date(data), na.rm = FALSE),
    prec = na.locf(prec, na.rm = FALSE),
    prec = na.locf(prec, fromLast = TRUE, na.rm = FALSE)
  ) %>%
  ungroup()

# Impute all_q_data_new
all_q_data_new_imputed <- all_q_data_new %>%
  group_by(location_name) %>%
  arrange(as.Date(data)) %>%
  mutate(
    Q = na.approx(Q, x = as.Date(data), na.rm = FALSE),
    Q = na.locf(Q, na.rm = FALSE),
    Q = na.locf(Q, fromLast = TRUE, na.rm = FALSE)
  ) %>%
  ungroup()

# Check for remaining NAs
sum(is.na(all_q_data_new_imputed$Q))  # 17202
sum(is.na(all_prec_data_new_imputed$prec))  # 1098

# --- Remove locations with all NA values in all_q_data_new ---
all_q_data_new_imputed <- all_q_data_new_imputed %>%
  group_by(location_name) %>%
  filter(!all(is.na(Q))) %>%  # Remove locations with only missing values
  ungroup()

# --- Remove locations with all NA values in all_prec_data_new ---
all_prec_data_new_imputed <- all_prec_data_new_imputed %>%
  group_by(location_name) %>%
  filter(!all(is.na(prec))) %>%  # Remove locations with only missing values
  ungroup()

# --- Impute remaining NAs in all_prec_data_new with location-specific mean ---
all_prec_data_new_imputed <- all_prec_data_new_imputed %>%
  group_by(location_name) %>%
  mutate(prec = ifelse(is.na(prec), mean(prec, na.rm = TRUE), prec)) %>%
  ungroup()

# Verify that no missing values remain
colSums(is.na(all_prec_data_new_imputed))
colSums(is.na(all_q_data_new_imputed))

# Center observations of river flow 
all_q_data_new_imputed <- all_q_data_new_imputed %>%
  mutate(logQ = log1p(Q)) %>%
  group_by(location_name) %>%
  mutate(Q_centered = logQ - mean(logQ, na.rm = TRUE)) %>%
  ungroup()

# ADD ELEVATION TO THE DATASETS 
library(geodata)
library(terra)
library(dplyr)

# STEP 1: Download both SRTM tiles for Switzerland and Italy
elev_che <- elevation_30s(country = "CHE", path = tempdir())
elev_ita <- elevation_30s(country = "ITA", path = tempdir())

# STEP 2: Reproject both to LV95 (EPSG:2056)
elev_che_lv95 <- project(elev_che, "EPSG:2056")
elev_ita_lv95 <- project(elev_ita, "EPSG:2056")

# STEP: Resample Italian raster to match Swiss raster (same resolution + extent grid)
elev_ita_lv95_resampled <- resample(elev_ita_lv95, elev_che_lv95)

# Now mosaic works!
elev_combined_lv95 <- mosaic(elev_che_lv95, elev_ita_lv95_resampled)


# STEP 4: Create expanded bounding box from your data
coords_mat <- as.matrix(all_prec_data_new_imputed[, c("coordinates.x", "coordinates.y")])
pts_vect <- vect(coords_mat, type = "points", crs = "EPSG:2056")

bbox <- ext(pts_vect)
bbox_expanded <- ext(bbox$xmin - 20000, bbox$xmax + 20000,
                     bbox$ymin - 20000, bbox$ymax + 20000)

# STEP 5: Crop to expanded bounding box (includes parts of Italy)
elev_expanded <- crop(elev_combined_lv95, bbox_expanded)

# STEP 7: Extract elevation at precipitation points
prec_points <- vect(all_prec_data_new_imputed, geom = c("coordinates.x", "coordinates.y"), crs = "EPSG:2056")
prec_elev <- extract(elev_expanded, prec_points)

# STEP 8: Merge into original data
all_prec_data_new_with_elev <- bind_cols(all_prec_data_new_imputed, elevation = prec_elev[, 2])

# same for river flow ----
# STEP 1: Create bounding box around river flow station coordinates ----
coords_q <- as.matrix(all_q_data_new_imputed[, c("coordinates.x", "coordinates.y")])
q_pts_vect <- vect(coords_q, type = "points", crs = "EPSG:2056")

# STEP 2: Expand bounding box by 20km (same as for precipitation) ----
bbox_q <- ext(q_pts_vect)
bbox_q_expanded <- ext(bbox_q$xmin - 20000, bbox_q$xmax + 20000,
                       bbox_q$ymin - 20000, bbox_q$ymax + 20000)

# STEP 3: Crop elevation raster to this expanded extent ----
elev_expanded_q <- crop(elev_combined_lv95, bbox_q_expanded)

# STEP 5: Convert Q station dataframe to SpatVector
q_points <- vect(all_q_data_new_imputed,
                 geom = c("coordinates.x", "coordinates.y"),
                 crs = "EPSG:2056")

# STEP 6: Extract elevation values at each river flow location ----
q_elev <- extract(elev_expanded_q, q_points)

# STEP 7: Merge elevation into the river flow dataset ----
all_q_data_new_with_elev <- bind_cols(all_q_data_new_imputed, elevation = q_elev[, 2])

# BACK TO TRAINING
###############################################################################
# 3) Restructure all_q_data_with_elev by Unique Dates
###############################################################################
y.n <- list()
all_dates <- unique(all_q_data_with_elev$data)

for (i in seq_along(all_dates)) {
  current_date <- all_dates[i]
  subset_data <- all_q_data_with_elev[all_q_data_with_elev$data == current_date, ]
  
  # Include elevation in sub_df
  sub_df <- data.frame(
    location_name = subset_data$location_name,
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    elevation     = subset_data$elevation,      # NEW
    Q             = subset_data$Q,
    data          = subset_data$data, 
    Q_centered    = subset_data$Q_centered
  )
  
  y.n[[i]] <- sub_df
}

str(y.n[1:5])  # Check structure

###############################################################################
# 4) Restructure all_prec_data_with_elev by Unique Dates
###############################################################################
x.n <- list()
all_dates_prec <- unique(all_prec_data_with_elev$data)

for (i in seq_along(all_dates_prec)) {
  current_date <- all_dates_prec[i]
  subset_data <- all_prec_data_with_elev[all_prec_data_with_elev$data == current_date, ]
  
  # Include elevation in sub_df
  sub_df <- data.frame(
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    elevation     = subset_data$elevation,      # NEW
    prec          = subset_data$prec,
    data          = subset_data$data
  )
  
  x.n[[i]] <- sub_df
}

str(x.n[1:5])  # Check structure


###############################################################################
# GAM Model for Q (River Flow) ----
###############################################################################
y.names <- "Q_centered"
p <- length(y.names)
N <- length(y.n)

models.plot <- vector(mode = "list", length = p)
names(models.plot) <- y.names

max_Tpoints <- max(sapply(y.n, nrow))
check.weight.y <- weight.y <- Y.star <- array(
  NA,
  dim = c(p, N, max_Tpoints),
  dimnames = list(y.names, paste0("N_", seq_len(N)), paste0("T_", seq_len(max_Tpoints)))
)

no.unit.y <- c()
h <- 0

for (j in seq_len(p)) {
  models.plot[[y.names[j]]] <- list()
  
  for (i in seq_len(N)) {
    subset_data <- y.n[[i]]
    y  <- as.data.frame(subset_data[, y.names[j]])
    x1 <- as.data.frame(subset_data$coordinates.x)
    x2 <- as.data.frame(subset_data$coordinates.y)
    elev <- as.data.frame(subset_data$elevation)  # NEW
    
    data <- data.frame("y" = y,
                       "x1" = as.numeric(unlist(x1)),
                       "x2" = as.numeric(unlist(x2)),
                       "elev" = as.numeric(unlist(elev)))  # NEW
    colnames(data) <- c("y", "x1", "x2", "elev")
    
    unique_x1 <- length(unique(data$x1))
    unique_x2 <- length(unique(data$x2))
    unique_elev <- length(unique(data$elev))  # NEW
    
    # Proceed only if enough unique points
    if (unique_x1 > 1 && unique_x2 > 1 && unique_elev > 1) {
      k_val  <- 25
      data$y <- log1p(data$y)
      
      # Try fitting the model with an added smooth term for elevation
      prova <- try(
        gam(y ~ s(x1, x2, elev, k = k_val), data = data, family = gaussian())
      )
      
      if (!inherits(prova, "try-error")) {
        model_i <- gam(y ~ s(x1, x2, elev, k = k_val), data = data, family = gaussian())
        
        models.plot[[y.names[j]]][[i]] <- list(model = model_i)
        
        # Predict only at the river flow station coordinates (for this day)
        new_data_pred <- data.frame(
          x1   = data$x1,
          x2   = data$x2,
          elev = data$elev  # NEW
        )
        
        predictions <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        
        w_ji <- 1 / ((predictions$se.fit)^2)
        w_ji[is.na(w_ji) | is.infinite(w_ji)] <- 0
        
        num_points <- length(predictions$fit)
        weight.y[y.names[j], i, 1:num_points] <- w_ji
        Y.star[y.names[j], i, 1:num_points]   <- predictions$fit
        
      } else {
        weight.y[y.names[j], i, ] <- NA
        Y.star[y.names[j], i, ]   <- NA
        h <- h + 1
        no.unit.y[h] <- i
      }
    } else {
      weight.y[y.names[j], i, ] <- NA
      Y.star[y.names[j], i, ]   <- NA
      h <- h + 1
      no.unit.y[h] <- i
    }
  }
}

print(no.unit.y)

###############################################################################
# GAM MODEL for X (Precipitation) with Elevation as an Additional Regressor
###############################################################################
library(mgcv)

# Define the dependent variable name for precipitation
x.names <- "prec"
q <- length(x.names)
N <- length(x.n)  # Total number of dates in the x.n list

# Initialize the list to store models for precipitation
models.plot.x <- vector(mode = "list", length = q)
names(models.plot.x) <- x.names

no.unit.x <- c()
h <- 0

# Ensure the prediction grid 'new.data' has an elevation column.
# If not available, assign the mean elevation from all_prec_data_with_elev.
if (!("elevation" %in% colnames(new.data))) {
  mean_elev <- mean(all_prec_data_with_elev$elevation, na.rm = TRUE)
  new.data$elevation <- mean_elev
}

# Define the number of grid points (N_T)
N_T <- nrow(new.data)

# Initialize arrays to store predictions, weights, and standard errors.
check.weight.x <- weight.x <- X.star <- array(
  NA, 
  dim = c("var" = q, "unit" = N, "Tpoint" = N_T), 
  dimnames = list(x.names, paste0("N_", seq_len(N)), paste0("T_", seq_len(N_T)))
)

# Loop over each date (i.e., each sub-dataset in x.n)
for (j in seq_len(q)) {
  models.plot.x[[x.names[j]]] <- list()  # Initialize model storage for this variable
  
  for (i in seq_len(N)) {
    # Extract the current sub-dataset for the date
    subset_data <- x.n[[i]]
    
    # Prepare the data for model fitting:
    # Convert the relevant columns to numeric vectors.
    y   <- as.data.frame(subset_data[, x.names[j]])
    x1  <- as.data.frame(subset_data$coordinates.x)
    x2  <- as.data.frame(subset_data$coordinates.y)
    elev <- as.data.frame(subset_data$elevation)
    
    data <- data.frame(
      y    = as.numeric(unlist(y)),
      x1   = as.numeric(unlist(x1)),
      x2   = as.numeric(unlist(x2)),
      elev = as.numeric(unlist(elev))
    )
    
    # Check for enough unique values for each regressor
    unique_x1   <- length(unique(data$x1))
    unique_x2   <- length(unique(data$x2))
    unique_elev <- length(unique(data$elev))
    
    # Only proceed if there is variability in all predictors
    if (unique_x1 > 1 && unique_x2 > 1 && unique_elev > 1) { 
      # Set smoothing parameters based on the number of unique values
      k_val  <- 25
      data$y <- log1p(data$y)
      
      # Try fitting the GAM model including elevation as an additional smooth term.
      prova <- try(
        gam(y ~ s(x1, x2, elev, k = k_val), data = data, family = gaussian())
      )
      
      if (!inherits(prova, "try-error")) {
        model_i <- gam(y ~ s(x1, x2, elev, k = k_val), data = data, family = gaussian())
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        
        # Prepare the newdata for predictions.
        # Here, new.data is our grid; ensure it has a numeric elevation column.
        new_data_pred <- data.frame(
          x1   = as.numeric(new.data$coordinates.x),
          x2   = as.numeric(new.data$coordinates.y),
          elev = as.numeric(new.data$elevation)
        )
        
        # Generate predictions over the new grid along with standard errors.
        predictions <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        
        # Compute weights using the inverse squared standard error.
        w_ji <- 1 / ((predictions$se.fit)^2)
        w_ji[is.na(w_ji) | is.infinite(w_ji)] <- 0
        
        # Store the weights and predicted values.
        weight.x[x.names[j], i, ] <- w_ji
        X.star[x.names[j], i, ]   <- predictions$fit
        
      } else {
        # In case of model fitting failure, store NA values.
        weight.x[x.names[j], i, ] <- NA
        X.star[x.names[j], i, ]   <- NA
        h <- h + 1
        no.unit.x[h] <- i
      }
      
      # Clean up temporary variables.
      rm(unique_x1, unique_x2, unique_elev, k_val)
      
    } else {
      # If there is not enough variability, assign NA.
      weight.x[x.names[j], i, ] <- NA
      X.star[x.names[j], i, ]   <- NA
      h <- h + 1
      no.unit.x[h] <- i
    }
  }
}

# Print indices of any failed dates (should be NULL if all models fitted successfully)
print(no.unit.x)

# Summaries
sum(is.na(weight.x))
sum(is.na(weight.y))

# --- Calculate Weighted Mean Predictions for Q (Y.mean) ---
max_Tpoints <- max(sapply(y.n, function(subset_data) nrow(subset_data)))  # Dynamic dimension calculation

Y.mean <- matrix(0, p, max_Tpoints)
rownames(Y.mean) <- y.names
colnames(Y.mean) <- paste0("T_", seq(1, max_Tpoints))

for (j in 1:p) {
  for (i in 1:N) {
    num_points <- sum(!is.na(weight.y[y.names[j], i, ]))  # Number of valid points for this day
    if (num_points > 0) {
      Y.mean[y.names[j], 1:num_points] <- Y.mean[y.names[j], 1:num_points] + 
        weight.y[y.names[j], i, 1:num_points] * Y.star[y.names[j], i, 1:num_points]
    }
  }
  
  # Normalize by the sum of weights
  weight_sum <- apply(weight.y[y.names[j], , ], 2, sum, na.rm = TRUE)
  weight_sum[weight_sum == 0] <- NA  # Avoid division by zero
  Y.mean[y.names[j], ] <- Y.mean[y.names[j], ] / weight_sum
}

# --- Calculate Weighted Mean Predictions for Prec (X.mean) ---
X.mean <- matrix(0, q, N_T)
rownames(X.mean) <- x.names
colnames(X.mean) <- paste0("T_", seq(1, N_T))

for (j in 1:q) {
  for (i in 1:N) {
    # Check if weight.x and X.star have valid (non-NA) values
    if (all(!is.na(weight.x[x.names[j], i, ])) && all(!is.na(X.star[x.names[j], i, ]))) {
      X.mean[x.names[j], ] <- X.mean[x.names[j], ] + weight.x[x.names[j], i, ] * X.star[x.names[j], i, ]
    }
  }
  
  # Normalize by the sum of weights
  weight_sum <- apply(weight.x[x.names[j], , ], 2, sum)
  X.mean[x.names[j], ] <- X.mean[x.names[j], ] / weight_sum
}



###################
### TWO BASIS SYSTEMS
###################

## - SIGMAj MATRIX : Covariance matrix for each variable - dimension N_T x N_T 
## - For Y (River Flow)sigma.Y_j <- vector(mode = "list", length = p)
sigma.Y_j <- vector(mode = "list", length = p)
names(sigma.Y_j) <- y.names

for (j in seq_len(p)) {
  max_Tpoints <- max(sapply(y.n, function(subset_data) nrow(subset_data)))  # Dynamic dimension calculation
  
  sigma.Y_ji_num <- sigma.Y_ji_den <- matrix(0, max_Tpoints, max_Tpoints)
  
  for (i in seq_len(N)) {
    current_weight_y <- weight.y[y.names[j], i, ]
    current_Y_star <- Y.star[y.names[j], i, ]
    
    # Check how many valid points exist for this particular day
    num_points <- sum(!is.na(current_weight_y))
    
    if (num_points > 0) {  # Proceed if weights are valid
      # Only consider the first 'num_points' values
      valid_weights <- current_weight_y[1:num_points]
      valid_Y_star <- current_Y_star[1:num_points]
      valid_Y_mean <- Y.mean[y.names[j], 1:num_points]
      
      A <- matrix(valid_weights * (valid_Y_star - valid_Y_mean), nrow = 1)
      sigma.Y_ji_num <- sigma.Y_ji_num + t(A) %*% A 
      sigma.Y_ji_den <- sigma.Y_ji_den + outer(valid_weights, valid_weights)
    }
  }
  
  sigma.Y_j[[y.names[j]]] <- sigma.Y_ji_num / sigma.Y_ji_den  # Store the covariance matrix for variable j
  print(paste("Completed calculation for Y variable:", y.names[j]))
}

## - For X (Precipitation)
sigma.X_j <- vector(mode = "list", length = q)
names(sigma.X_j) <- x.names

for (j in seq_len(q)) {
  sigma.X_ji_num <- sigma.X_ji_den <- matrix(0, N_T, N_T)
  
  for (i in seq_len(N)) {
    current_weight_x <- weight.x[x.names[j], i, ]
    current_X_star <- X.star[x.names[j], i, ]
    
    if (!all(is.na(current_weight_x))) {  # Proceed if weights are valid
      A <- matrix(current_weight_x * (current_X_star - X.mean[x.names[j], ]), nrow = 1)
      sigma.X_ji_num <- sigma.X_ji_num + t(A) %*% A 
      sigma.X_ji_den <- sigma.X_ji_den + outer(current_weight_x, current_weight_x)
    }
  }
  
  sigma.X_j[[x.names[j]]] <- sigma.X_ji_num / sigma.X_ji_den  # Store the covariance matrix for variable j
  print(paste("Completed calculation for X variable:", x.names[j]))
}

## - Compute H.Y and H.X (Mean covariance matrices)
H.X <- matrix(0, N_T, N_T)

# Sum over all variables to get the mean covariance matrix
# --- Initialize H.Y with the correct dimension ---
max_Tpoints <- max(sapply(y.n, function(subset_data) nrow(subset_data)))
H.Y <- matrix(0, max_Tpoints, max_Tpoints)

# --- Aggregate covariance matrices for Y (River Flow) ---
for (j in seq_len(p)) {
  if (!is.null(sigma.Y_j[[j]])) {  # Check if the matrix was calculated successfully
    H.Y <- H.Y + sigma.Y_j[[j]]
  }
}

# --- Normalize H.Y by the number of variables (p) ---
H.Y <- H.Y / p  # This is now a proper average, based on the actual sizes of the covariance matrices


for (j in seq_len(q)) {
  H.X <- H.X + (1 / q) * sigma.X_j[[j]]
}

# Display ranges to check for correctness
range(H.Y) # -0.02588097  0.24368886
range(H.X) # -0.0007313096  0.0074061257

# ORTHONORMAL BASIS ----
library(splines)

# Define the number of basis functions you want to use
L <- 10  # For Q (river flow)
H <-10  # For Prec (precipitation)

# Generating basis functions for Q (River Flow)
phi.Y <- bs(seq(0, 1, length.out = nrow(H.Y)), df = L, intercept = TRUE)


# Generating basis functions for Prec (Precipitation)
phi.X <- bs(seq(0, 1, length.out = nrow(H.X)), df = H, intercept = TRUE)

# SCORES ----
# --- Initialize Score Arrays ---
gamma <- array(NA, dim = c(p, N, L), 
               dimnames = list(y.names, paste0("N_", seq(1, N)), paste0("L_", seq(1, L))))
chi <- array(NA, dim = c(q, N, H), 
             dimnames = list(x.names, paste0("N_", seq(1, N)), paste0("H_", seq(1, H))))

# --- Calculate Scores for Y (Q) ---
for (j in seq_len(p)) {
  for (i in seq_len(N)) {
    current_weight_y <- weight.y[j, i, ]
    current_Y_star <- Y.star[j, i, ]
    
    # Number of valid points for this day
    num_points <- sum(!is.na(current_weight_y))
    
    if (num_points > 0) {
      # Extract relevant portions of the weight and data
      valid_weights <- current_weight_y[1:num_points]
      valid_Y_star <- current_Y_star[1:num_points]
      valid_Y_mean <- Y.mean[j, 1:num_points]
      
      # Perform the basis function decomposition
      tryCatch({
        gamma[j, i, ] <- solve(t(phi.Y) %*% diag(valid_weights) %*% phi.Y) %*% 
          t(phi.Y) %*% diag(valid_weights) %*% (valid_Y_star - valid_Y_mean)
      }, error = function(e) {
        gamma[j, i, ] <- rep(0, L)  # Fill with zeros if calculation fails
      })
    } else {
      gamma[j, i, ] <- rep(0, L)  # If no valid points, assign zeros
    }
  }
}



# --- Calculate Scores for X (Prec) ---
for (j in seq_len(q)) {
  for (i in seq_len(N)) {
    # 1) Extract the day's data
    subset_data <- x.n[[i]]
    
    # 2) FALLBACK: if there's no spatial variation, give a uniform first‐basis score
    if (var(subset_data$prec, na.rm = TRUE) == 0) {
      # assign all mass to the first basis function
      chi[j, i, ] <- c(1, rep(0, H - 1))
      next  # skip the rest of this loop iteration
    }
    
    # 3) Otherwise proceed with your usual weight/X-star extraction
    current_weight_x <- weight.x[j, i, ]
    current_X_star   <- X.star[j, i, ]
    
    # 4) Now do the solve() … code you already have
    tryCatch({
      chi[j, i, ] <- solve(
        t(phi.X) %*% diag(current_weight_x) %*% phi.X
      ) %*% t(phi.X) %*% diag(current_weight_x) %*% (current_X_star - X.mean[j, ])
    }, error = function(e) {
      chi[j, i, ] <- rep(0, H)  # optional: keep your original error‐fallback
    })
  }
}

# --- Check Scores ---
which(apply(chi[1,,], 1, function(h) sum(h)==0)) # 0
which(apply(gamma[1,,], 1, function(h) sum(h)==0)) # 0

# REGRESSION ----

# STEP 1: Construct lagged chi scores (chi_{t-1})
x_lag <- matrix(NA, N, q * H)
colnames(x_lag) <- paste0("lag.h.", sort(rep(seq_len(H), q)), "_", rep(x.names, H))
for (h in seq_len(H)) {
  x_lag[, h] <- chi[1, , h]
}
x_lag <- rbind(rep(NA, q * H), x_lag[-N, ])  # shift down by one row

# STEP 2: Construct lagged gamma scores (gamma_{t-1})
gamma_lag <- matrix(NA, N, p * L)
colnames(gamma_lag) <- paste0("lag.l.", sort(rep(seq_len(L), p)), "_", rep(y.names, L))
for (l in seq_len(L)) {
  gamma_lag[, l] <- gamma[1, , l]
}
gamma_lag <- rbind(rep(NA, p * L), gamma_lag[-N, ])  # shift down by one row

# STEP 3: Build the full design matrix x (lagged precipitation + lagged river flow)
x <- cbind(x_lag, gamma_lag)
x <- x[-1, ]  # drop first row due to NA lags
x <- cbind(1, x)  # add intercept

# STEP 4: Build response matrix y (today's gamma)
y <- matrix(NA, N, p * L)
colnames(y) <- paste0("l.", sort(rep(seq_len(L), p)), "_", rep(y.names, L))
for (l in seq_len(L)) {
  y[, l] <- gamma[1, , l]
}
y <- y[-1, ]  # align with x (remove first row)

# STEP 5: Fit the penalized regression model
library(cglasso)
data <- datacggm(Y = y, X = x)
model <- cglasso(. ~ ., data = data, nlambda = 2, lambda.min.ratio = 0, nrho = 1)

# STEP 6: Extract and reshape the coefficient matrix
B_vector <- model$B[-1, , 2, 1]  # drop intercept row from model$B
B_matrix <- matrix(B_vector, nrow = ncol(x), ncol = ncol(y))

# STEP 7: Predict gamma and compute residuals
gamma_pred <- x %*% B_matrix
residuals <- y - gamma_pred

# STEP 8: Compute Mean Squared Error
MSE <- mean(residuals^2, na.rm = TRUE)
cat("Mean Squared Error (MSE):", MSE, "\n") # 0.03464865 

# STEP 9: Plotting predicted vs actual gamma
plot(y, gamma_pred,
     xlab = "Actual Gamma Values", ylab = "Predicted Gamma Values",
     main = "Predicted vs Actual Gamma Values",
     col = "blue", pch = 19, 
     cex.main = 1)
abline(0, 1, col = "red")

# MODEL RESIDUALS FROM MODEL 3 ----
# FIT AR(1) TO THE RESIDUALS 
# to model 3, ie model with elevation option 1 
# Step: Fit AR(1) to residuals and extract innovation noise
library(stats)

# Initialize containers for AR(1) coefficient (phi) and innovation SD (tau)
phi <- numeric(ncol(residuals))       # AR(1) coefficients for each basis function
tau <- numeric(ncol(residuals))       # Innovation (white noise) standard deviations
eta_matrix <- matrix(NA, nrow = N - 1, ncol = ncol(residuals))  # Matrix of innovations η_t

# Loop over each basis function
for (l in seq_len(ncol(residuals))) {
  res_l <- residuals[, l]
  
  # Fit AR(1) model without intercept (residuals already centered)
  fit <- try(arima(res_l, order = c(1, 0, 0), include.mean = FALSE), silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    phi[l] <- fit$coef["ar1"]
    tau[l] <- sqrt(fit$sigma2)
    
    # Compute η_t = ε_t - φ ε_{t-1}
    eps_t   <- res_l[2:N]
    eps_tm1 <- res_l[1:(N - 1)]
    eta_matrix[, l] <- eps_t - phi[l] * eps_tm1
  } else {
    # Fallback: set to zero if AR(1) fit fails
    phi[l] <- 0
    tau[l] <- NA
    eta_matrix[, l] <- NA
  }
}

# Compute innovation-based RMSE
RMSE_AR1 <- mean(eta_matrix^2, na.rm = TRUE)
cat("RMSE from AR(1) innovations:", round(RMSE_AR1, 5), "\n") #  0.03443 

# Optional: summary of AR(1) coefficients
summary(phi)

# PREDICTIONS ----

# --- STEP 1: Retrieve latest day index (assumes gamma and chi are aligned)
last_day_index <- dim(gamma)[2]  # Last day of training

# --- STEP 2: Get today's precipitation scores (chi_t)
chi_t <- chi[1, last_day_index, ]  # [length H]

# --- STEP 3: Get today's river flow scores (gamma_t)
gamma_t <- gamma[1, last_day_index, ]  # [length L]

# --- STEP 4: Construct design vector for tomorrow
X_forecast <- c(
  1,         # Intercept
  chi_t,     # Lagged precipitation
  gamma_t    # Lagged river flow
)

# --- STEP 5: Predict tomorrow's gamma scores
gamma_pred_next <- X_forecast %*% B_matrix  # [1 x L]

# --- STEP 6: Apply AR(1) residual correction
resid_corr <- numeric(L)
for (l in seq_len(L)) {
  gamma_t_minus_1 <- gamma[1, last_day_index - 1, l]
  gamma_t_current <- gamma[1, last_day_index, l]
  resid_corr[l] <- -phi[l] * (gamma_t_current - gamma_t_minus_1)
}
gamma_pred_next <- gamma_pred_next + resid_corr  # [1 x L]

# --- STEP 7: Reconstruct predicted log river flow
Q_pred_log <- gamma_pred_next %*% t(phi.Y) + Y.mean  # [1 x Tpoints]

# --- STEP 7.1: Compute prediction standard error in log-scale
Q_var_log <- rowSums((phi.Y^2) * matrix(tau^2, nrow = nrow(phi.Y), ncol = length(tau), byrow = TRUE))  # [Tpoints]
Q_se_log <- sqrt(Q_var_log)

# --- STEP 7.2: Compute 95% confidence intervals in log scale
z_score <- qnorm(0.975)  # 1.96 for 95% CI
Q_pred_log_upper <- Q_pred_log + z_score * Q_se_log
Q_pred_log_lower <- Q_pred_log - z_score * Q_se_log

# --- STEP 8: Add station-specific mean and back-transform to Q
station_means <- all_q_data_imputed %>%
  mutate(logQ = log1p(Q)) %>%
  group_by(location_name) %>%
  summarize(mean_logQ = mean(logQ, na.rm = TRUE)) %>%
  arrange(location_name)

mean_logQ_vec <- station_means$mean_logQ  # [Tpoints]

# Back-transform predictions
Q_forecast <- exp(Q_pred_log + mean_logQ_vec) - 1
Q_lower <- exp(Q_pred_log_lower + mean_logQ_vec) - 1
Q_upper <- exp(Q_pred_log_upper + mean_logQ_vec) - 1

# --- STEP 9: Build forecast dataframe

# 1. Get forecast date from test data
forecast_date <- min(all_q_data_new_with_elev$data)

# 2. Get station names
stations <- unique(all_q_data_new_with_elev$location_name)
stations <- stations[stations != "Reuss Briggboden"]

# 3. Extract actual Q for today
Q_today_df <- all_q_data_new_with_elev %>%
  filter(data == forecast_date) %>%
  arrange(location_name) %>%
  select(location_name, Q_true = Q)

# 4. Extract station coordinates
coords_x <- all_q_data_with_elev %>%
  group_by(location_name) %>%
  summarize(coordinates.x = first(coordinates.x), .groups = "drop")

coords_y <- all_q_data_with_elev %>%
  group_by(location_name) %>%
  summarize(coordinates.y = first(coordinates.y), .groups = "drop")

# 5. Assemble forecast table
forecast_df <- data.frame(
  date = forecast_date,
  location_name = stations,
  Q_forecast = as.numeric(Q_forecast),
  Q_lower = as.numeric(Q_lower),
  Q_upper = as.numeric(Q_upper)
) %>%
  left_join(Q_today_df, by = "location_name") %>%
  left_join(coords_x, by = "location_name") %>%
  left_join(coords_y, by = "location_name") %>%
  relocate(date, location_name, coordinates.x, coordinates.y, Q_true, Q_forecast, Q_lower, Q_upper)

# 6. Print the forecast table with 95% CI
print(forecast_df)


# PLOT 2 ----
library(ggplot2)
library(ggrepel)  # for better label placement

ggplot(forecast_df, aes(x = Q_true, y = Q_forecast)) +
  geom_point(color = "black", size = 1.5) +   # small black dots
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +  # 45-degree line
  geom_text_repel(aes(label = location_name), size = 3, max.overlaps = 10) +  # labels near points
  labs(
    title = "True vs Forecasted River Flow",
    x = "True River Flow",
    y = "Forecasted River Flow"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


# plot 
ggplot(forecast_df) + geom_point(aes(x = Q_true, y = Q_forecast)) + geom_abline()

# PLOTS OF 6 MONTHS OF TRAINING AND PREDICTION FOT 2024-01-01 ----
# 1. Prepare training data (last 6 months)
training_df <- all_q_data_with_elev %>%
  filter(data >= as.Date("2023-07-01") & data <= as.Date("2023-12-31")) %>%
  rename(date = data, Q_training = Q) %>%
  mutate(Q_forecast = NA, Q_lower = NA, Q_upper = NA, Q_true = NA)

# 2. Prepare forecast data (2024-01-01)
forecast_df_long <- forecast_df %>%
  rename(date = date) %>%
  mutate(Q_training = NA)  # No training data for forecast date

# 3. Combine datasets
combined_df <- bind_rows(
  training_df %>%
    select(location_name, date, Q_training, Q_forecast, Q_lower, Q_upper, Q_true),
  forecast_df_long %>%
    select(location_name, date, Q_training, Q_forecast = Q_forecast, Q_lower, Q_upper, Q_true)
)

# 5. Initialize empty list to collect plots
plot_list <- list()

for (station in stations) {
  df_station <- combined_df %>% filter(location_name == station)
  
  # Reshape forecast points for consistent legends
  forecast_points <- df_station %>% filter(!is.na(Q_forecast))
  true_points <- df_station %>% filter(!is.na(Q_true))
  
  p <- ggplot(df_station, aes(x = date)) +
    # Training line
    geom_line(aes(y = Q_training, color = "Training"), linewidth = 0.5, na.rm = TRUE) +
    
    # Forecast point
    geom_point(
      data = forecast_points,
      aes(y = Q_forecast, color = "Forecast"),
      size = 3
    ) +
    
    # Forecast error bar
    geom_errorbar(
      data = forecast_points,
      aes(ymin = Q_lower, ymax = Q_upper, color = "Forecast"),
      linewidth = 0.9
    ) +
    
    # True observed point
    geom_point(
      data = true_points,
      aes(y = Q_true, color = "True"),
      size = 3
    ) +
    
    # Labels
    labs(
      title = paste("River Flow Time Series with Forecast and CI -", station),
      x = "Date",
      y = "River Flow (Q)",
      color = NULL
    ) +
    
    # Custom colors for legend
    scale_color_manual(
      values = c(
        "Training" = "black",
        "Forecast" = "blue",
        "True" = "red"
      )
    ) +
    
    # Styling
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 18,
        hjust = 0
      ),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.position = "bottom",
      legend.text = element_text(size = 14)
    )
  
  # Store plot
  plot_list[[station]] <- p
}


# 7. Display all plots (as a list, for example)
# This prints them one after the other
for (p in plot_list) print(p)



