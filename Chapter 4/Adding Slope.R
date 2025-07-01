rm(list = ls(all = TRUE))     
cat("\014")   

# Load required libraries ---
library(dplyr)
library(tidyr)
library(imputeTS)
library(ggplot2)
library(mgcv)
library(splines)
library(cglasso)
library(forecast)     
library(gridExtra)
library(grid)
library(ggpubr)
library(geodata)
library(terra)
library(stats)
library(FNN)
library(dplyr)
library(httr)            
library(jsonlite)           
library(zoo)

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
# Delete the outliers
all_prec_data$prec[all_prec_data$prec > 500] <- NA
all_q_data$Q[all_q_data$Q > 1000] <- NA

# Impute all_prec_data
all_prec_data_imputed <- all_prec_data %>%
  group_by(location_name) %>%
  arrange(as.Date(data)) %>%
  mutate(
    prec = na.approx(prec, x = as.Date(data), na.rm = FALSE),
    prec = na.locf(prec, na.rm = FALSE),
    prec = na.locf(prec, fromLast = TRUE, na.rm = FALSE)
  ) %>%
  ungroup()

# Impute all_q_data
all_q_data_imputed <- all_q_data %>%
  group_by(location_name) %>%
  arrange(as.Date(data)) %>%
  mutate(
    Q = na.approx(Q, x = as.Date(data), na.rm = FALSE),
    Q = na.locf(Q, na.rm = FALSE),
    Q = na.locf(Q, fromLast = TRUE, na.rm = FALSE)
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
# STEP 1: Download both SRTM tiles for Switzerland and Italy
elev_che <- elevation_30s(country = "CHE", path = tempdir())
elev_ita <- elevation_30s(country = "ITA", path = tempdir())

# STEP 2: Reproject both to LV95 (EPSG:2056)
elev_che_lv95 <- project(elev_che, "EPSG:2056")
elev_ita_lv95 <- project(elev_ita, "EPSG:2056")

# STEP: Resample Italian raster to match Swiss raster (same resolution + extent grid)
elev_ita_lv95_resampled <- resample(elev_ita_lv95, elev_che_lv95)

elev_combined_lv95 <- mosaic(elev_che_lv95, elev_ita_lv95_resampled)


# STEP 4: Create expanded bounding box from data
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

# adding slope for precipitation ----
# Unique precip stations
prec_stations <- all_prec_data_with_elev %>%
  distinct(location_name, coordinates.x, coordinates.y, elevation)

# Extract coordinates and elevation
coords_prec <- prec_stations %>% select(coordinates.x, coordinates.y)
elev_prec <- prec_stations$elevation
n_prec <- nrow(prec_stations)

# Nearest neighbors (k = 7 including itself)
nn_prec <- get.knn(coords_prec, k = 7)
neighbor_indices_prec <- nn_prec$nn.index

# Compute slope for each precipitation station
slope_prec <- numeric(n_prec)

for (i in 1:n_prec) {
  neigh_ids <- neighbor_indices_prec[i, 2:7]
  pairs <- combn(neigh_ids, 2, simplify = FALSE)
  
  slopes <- sapply(pairs, function(pair) {
    h <- pair[1]; l <- pair[2]
    
    if (is.na(elev_prec[h]) || is.na(elev_prec[l])) return(NA_real_)
    
    dx <- coords_prec[h, 1] - coords_prec[l, 1]
    dy <- coords_prec[h, 2] - coords_prec[l, 2]
    dist_hl <- sqrt(dx^2 + dy^2)
    
    if (dist_hl == 0) return(NA_real_)
    
    return(abs(elev_prec[h] - elev_prec[l]) / dist_hl)
  })
  
  # Remove NAs before averaging
  slopes_numeric <- unlist(slopes)
  slopes_numeric <- slopes_numeric[!is.na(slopes_numeric)]
  
  if (length(slopes_numeric) > 0) {
    slope_prec[i] <- mean(slopes_numeric)
  } else {
    slope_prec[i] <- NA_real_
    cat("⚠️ Station", i, "has no valid slopes\n")
  }
}

# Merge back
prec_stations$slope <- slope_prec

# Join into full data
all_prec_data_with_slope <- all_prec_data_with_elev %>%
  left_join(prec_stations %>% select(location_name, slope), by = "location_name")


# adding slope for river flow ----
library(FNN)
library(dplyr)

# Step 1: Extract unique flow stations with coordinates and elevation
q_stations <- all_q_data_with_elev %>%
  distinct(location_name, coordinates.x, coordinates.y, elevation)

coords_q <- q_stations %>% select(coordinates.x, coordinates.y)
elev_q <- q_stations$elevation
n_q <- nrow(q_stations)

# Step 2: Find 7 nearest neighbors (including self)
nn_q <- get.knn(coords_q, k = 7)
neighbor_indices_q <- nn_q$nn.index

# Step 3: Compute slope for each flow station
slope_q <- numeric(n_q)

for (i in 1:n_q) {
  neigh_ids <- neighbor_indices_q[i, 2:7]  # Exclude the station itself
  pairs <- combn(neigh_ids, 2, simplify = FALSE)
  
  slopes <- sapply(pairs, function(pair) {
    h <- pair[1]; l <- pair[2]
    
    if (is.na(elev_q[h]) || is.na(elev_q[l])) return(NA_real_)
    
    dx <- coords_q[h, 1] - coords_q[l, 1]
    dy <- coords_q[h, 2] - coords_q[l, 2]
    dist_hl <- sqrt(dx^2 + dy^2)
    
    if (dist_hl == 0) return(NA_real_)
    
    return(abs(elev_q[h] - elev_q[l]) / dist_hl)
  })
  
  slopes_numeric <- unlist(slopes)
  slopes_numeric <- slopes_numeric[!is.na(slopes_numeric)]
  
  if (length(slopes_numeric) > 0) {
    slope_q[i] <- mean(slopes_numeric)
  } else {
    slope_q[i] <- NA_real_
    cat("⚠️ Station", i, "has no valid slope value\n")
  }
}

# Step 4: Merge slope into station metadata
q_stations$slope <- slope_q

# Step 5: Join slope back into full Q dataset
all_q_data_with_slope <- all_q_data_with_elev %>%
  left_join(q_stations %>% select(location_name, slope), by = "location_name")


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

###############################################################################
# 3) Restructure all_q_data_with_slope by Unique Dates
###############################################################################
y.n <- list()
all_dates <- unique(all_q_data_with_slope$data)

for (i in seq_along(all_dates)) {
  current_date <- all_dates[i]
  subset_data <- all_q_data_with_slope[all_q_data_with_slope$data == current_date, ]
  
  # Include slope in sub_df
  sub_df <- data.frame(
    location_name = subset_data$location_name,
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    slope     = subset_data$slope,      # NEW
    Q             = subset_data$Q,
    data          = subset_data$data, 
    Q_centered    = subset_data$Q_centered
  )
  
  y.n[[i]] <- sub_df
}

str(y.n[1:5])  # Check structure

###############################################################################
# 4) Restructure all_prec_data_with_slope by Unique Dates
###############################################################################
x.n <- list()
all_dates_prec <- unique(all_prec_data_with_slope$data)

for (i in seq_along(all_dates_prec)) {
  current_date <- all_dates_prec[i]
  subset_data <- all_prec_data_with_slope[all_prec_data_with_slope$data == current_date, ]
  
  # Include slope in sub_df
  sub_df <- data.frame(
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    slope     = subset_data$slope,      # NEW
    prec          = subset_data$prec,
    data          = subset_data$data
  )
  
  x.n[[i]] <- sub_df
}

str(x.n[1:5])  # Check structure

# OPTION 1: Integrate slope into the 3D Smooth  ----
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
    slope <- as.data.frame(subset_data$slope)  # NEW
    
    data <- data.frame("y" = y,
                       "x1" = as.numeric(unlist(x1)),
                       "x2" = as.numeric(unlist(x2)),
                       "slope" = as.numeric(unlist(slope)))  # NEW
    colnames(data) <- c("y", "x1", "x2", "slope")
    
    unique_x1 <- length(unique(data$x1))
    unique_x2 <- length(unique(data$x2))
    unique_slope <- length(unique(data$slope))  # NEW
    
    # Proceed only if enough unique points
    if (unique_x1 > 1 && unique_x2 > 1 && unique_slope > 1) {
      k_val  <- 25
      data$y <- log1p(data$y)
      
      # Try fitting the model with an added smooth term for slope
      prova <- try(
        gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
      )
      
      if (!inherits(prova, "try-error")) {
        model_i <- gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
        
        models.plot[[y.names[j]]][[i]] <- list(model = model_i)
        
        # Predict only at the river flow station coordinates (for this day)
        new_data_pred <- data.frame(
          x1   = data$x1,
          x2   = data$x2,
          slope = data$slope  # NEW
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
# GAM MODEL for X (Precipitation) with Slope as an Additional Regressor
###############################################################################
# Define the dependent variable name for precipitation
x.names <- "prec"
q <- length(x.names)
N <- length(x.n)  # Total number of dates in the x.n list

# Initialize the list to store models for precipitation
models.plot.x <- vector(mode = "list", length = q)
names(models.plot.x) <- x.names

no.unit.x <- c()
h <- 0

# Ensure the prediction grid 'new.data' has a slope column.
# If not available, assign the mean slope from all_prec_data_with_slope.
if (!("slope" %in% colnames(new.data))) {
  mean_slope <- mean(all_prec_data_with_slope$slope, na.rm = TRUE)
  new.data$slope <- mean_slope
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
    
    # Prepare the data for model fitting
    y     <- as.data.frame(subset_data[, x.names[j]])
    x1    <- as.data.frame(subset_data$coordinates.x)
    x2    <- as.data.frame(subset_data$coordinates.y)
    slope <- as.data.frame(subset_data$slope)
    
    data <- data.frame(
      y     = as.numeric(unlist(y)),
      x1    = as.numeric(unlist(x1)),
      x2    = as.numeric(unlist(x2)),
      slope = as.numeric(unlist(slope))
    )
    
    # Check for enough unique values for each regressor
    unique_x1    <- length(unique(data$x1))
    unique_x2    <- length(unique(data$x2))
    unique_slope <- length(unique(data$slope))
    
    if (unique_x1 > 1 && unique_x2 > 1 && unique_slope > 1) { 
      k_val  <- 25
      data$y <- log1p(data$y)
      
      # Try fitting the GAM model with slope
      prova <- try(
        gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
      )
      
      if (!inherits(prova, "try-error")) {
        model_i <- gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        
        # Prepare the newdata for prediction
        new_data_pred <- data.frame(
          x1    = as.numeric(new.data$coordinates.x),
          x2    = as.numeric(new.data$coordinates.y),
          slope = as.numeric(new.data$slope)
        )
        
        predictions <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        
        # Compute weights using the inverse squared standard error
        w_ji <- 1 / ((predictions$se.fit)^2)
        w_ji[is.na(w_ji) | is.infinite(w_ji)] <- 0
        
        weight.x[x.names[j], i, ] <- w_ji
        X.star[x.names[j], i, ]   <- predictions$fit
        
      } else {
        weight.x[x.names[j], i, ] <- NA
        X.star[x.names[j], i, ]   <- NA
        h <- h + 1
        no.unit.x[h] <- i
      }
      
      rm(unique_x1, unique_x2, unique_slope, k_val)
      
    } else {
      weight.x[x.names[j], i, ] <- NA
      X.star[x.names[j], i, ]   <- NA
      h <- h + 1
      no.unit.x[h] <- i
    }
  }
}

# Print indices of any failed fits
print(no.unit.x)

# GAM MODEL FOR SLOPE ----
# Setup
models.plot.slope <- list()
no.unit.slope <- c()
h <- 0

# Ensure the prediction grid 'new.data' has a slope column
if (!("slope" %in% colnames(new.data))) {
  mean_slope <- mean(all_prec_data_with_slope$slope, na.rm = TRUE)
  new.data$slope <- mean_slope
}

# Define number of grid points
N_T <- nrow(new.data)

# Initialize arrays
check.weight.slope <- weight.slope <- slope.star <- array(
  NA, 
  dim = c("unit" = N, "Tpoint" = N_T), 
  dimnames = list(paste0("N_", seq_len(N)), paste0("T_", seq_len(N_T)))
)

# Loop over dates (or set i = 1 if slope is static)
for (i in seq_len(N)) {
  # Use slope as response variable
  subset_data <- x.n[[i]]
  
  y  <- as.numeric(subset_data$slope)
  x1 <- as.numeric(subset_data$coordinates.x)
  x2 <- as.numeric(subset_data$coordinates.y)
  
  data <- data.frame(y = y, x1 = x1, x2 = x2)
  
  # Check for enough variation
  if (length(unique(x1)) > 1 && length(unique(x2)) > 1 && length(unique(y)) > 1) {
    k_val <- 25
    
    # Fit model: slope ~ s(x1, x2)
    fit_try <- try(
      gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
    )
    
    if (!inherits(fit_try, "try-error")) {
      model_i <- gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
      models.plot.slope[[i]] <- list(model = model_i)
      
      # Predict on the grid
      new_data_pred <- data.frame(
        x1 = new.data$coordinates.x,
        x2 = new.data$coordinates.y
      )
      
      preds <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
      
      w_ji <- 1 / (preds$se.fit^2)
      w_ji[is.na(w_ji) | is.infinite(w_ji)] <- 0
      
      weight.slope[i, ] <- w_ji
      slope.star[i, ]   <- preds$fit
      
    } else {
      weight.slope[i, ] <- NA
      slope.star[i, ]   <- NA
      h <- h + 1
      no.unit.slope[h] <- i
    }
    
    rm(k_val)
    
  } else {
    weight.slope[i, ] <- NA
    slope.star[i, ]   <- NA
    h <- h + 1
    no.unit.slope[h] <- i
  }
}

# Print indices of any failed fits
print(no.unit.slope)

# Summary stats
sum(is.na(weight.x))
sum(is.na(weight.y))
sum(is.na(weight.slope))


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
range(H.Y) # -0.001658578  0.100102018
range(H.X) # -0.0002294762  0.0098118685

# ORTHONORMAL BASIS ----
# Define the number of basis functions 
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
      # Extract rslopeant portions of the weight and data
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
    
    # 3) Otherwise proceed with usual weight/X-star extraction
    current_weight_x <- weight.x[j, i, ]
    current_X_star   <- X.star[j, i, ]
    
    tryCatch({
      chi[j, i, ] <- solve(
        t(phi.X) %*% diag(current_weight_x) %*% phi.X
      ) %*% t(phi.X) %*% diag(current_weight_x) %*% (current_X_star - X.mean[j, ])
    }, error = function(e) {
      chi[j, i, ] <- rep(0, H)  
    })
  }
}

# --- Check Scores ---
which(apply(chi[1,,], 1, function(h) sum(h)==0)) # 0
which(apply(gamma[1,,], 1, function(h) sum(h)==0)) # 0

# --- Project SLOPE onto phi.X once (static) ---
# Use the first non-NA slope.star row (slope is same over time)
row_idx <- which(rowSums(!is.na(slope.star)) > 0)[1]
slope_surface <- slope.star[row_idx, ]

# Center slope surface
slope_centered <- slope_surface - mean(slope_surface, na.rm = TRUE)
slope_centered[is.na(slope_centered)] <- 0  # fallback

# Project onto phi.X
slope_scores <- as.vector(t(phi.X) %*% slope_centered)  # length H

# # REGRESSION ----
# --- STEP 1: Construct current-day chi scores (chi_t) ---
x_now <- matrix(NA, N, q * H)
colnames(x_now) <- paste0("h.", sort(rep(seq_len(H), q)), "_", rep(x.names, H))
for (h in seq_len(H)) {
  x_now[, h] <- chi[1, , h]
}

# --- STEP 2: Construct lagged chi scores (chi_{t-1}) ---
x_lag <- rbind(rep(NA, q * H), x_now[-N, ])
colnames(x_lag) <- paste0("lag.h.", sort(rep(seq_len(H), q)), "_", rep(x.names, H))

# --- STEP 3: Construct lagged gamma scores (gamma_{t-1}) ---
gamma_lag <- matrix(NA, N, p * L)
colnames(gamma_lag) <- paste0("lag.l.", sort(rep(seq_len(L), p)), "_", rep(y.names, L))
for (l in seq_len(L)) {
  gamma_lag[, l] <- gamma[1, , l]
}
gamma_lag <- rbind(rep(NA, p * L), gamma_lag[-N, ])

# --- STEP 4: Construct slope matrix (replicated for all time points) ---
slope_mat <- matrix(rep(slope_scores, N), nrow = N, byrow = TRUE)
slope_mat <- slope_mat[-1, ]  # drop first row to align with lagged predictors

# --- STEP 5: Build full design matrix x ---
x <- cbind(x_lag, x_now, gamma_lag)
x <- x[-1, ]  # drop first row due to lagged NA
x <- cbind(1, x, slope_mat)  # add intercept and slope scores

# --- STEP 6: Build response matrix y ---
y <- matrix(NA, N, p * L)
colnames(y) <- paste0("l.", sort(rep(seq_len(L), p)), "_", rep(y.names, L))
for (l in seq_len(L)) {
  y[, l] <- gamma[1, , l]
}
y <- y[-1, ]  # align with design matrix x

# --- STEP 7: Fit the penalized regression model using cglasso ---
data <- datacggm(Y = y, X = x)
model <- cglasso(. ~ ., data = data, nlambda = 2, lambda.min.ratio = 0, nrho = 1)

# --- STEP 8: Extract and reshape the coefficient matrix ---
B_vector <- model$B[-1, , 2, 1]  # exclude global intercept
B_matrix <- matrix(B_vector, nrow = ncol(x), ncol = ncol(y))

# --- STEP 9: Predict gamma and compute residuals ---
gamma_pred <- x %*% B_matrix
residuals <- y - gamma_pred

# --- STEP 10: Compute MSE ---
MSE <- mean(residuals^2, na.rm = TRUE)
cat("Mean Squared Error (MSE):", MSE, "\n") # 0.02261893

# Plotting Predicted vs Actual gamma values
plot(y, gamma_pred,
     xlab = "Actual Gamma Values", ylab = "Predicted Gamma Values",
     main = "Predicted vs Actual Gamma Values",
     col = "blue", pch = 19, 
     cex.main = 1)
abline(0, 1, col = "red")

# MODEL RESIDUALS FROM MODEL 3 ----
# FIT AR(1) TO THE RESIDUALS 
# to model 3, ie model with slope option 1 
# Step: Fit AR(1) to residuals and extract innovation noise

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
MSE_AR1 <- mean(eta_matrix^2, na.rm = TRUE)
cat("MSE from AR(1) innovations:", round(MSE_AR1, 5), "\n") # 0.02213 

# Summary of AR(1) coefficients
summary(phi)


# NEW ----
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


# STEP 4: Create expanded bounding box 
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

# DO EVERYTHING AGAIN FOR MY NEW DATA ----
lat_range <- range(all_q_data_new_imputed$coordinates.y)
long_range <- range(all_q_data_new_imputed$coordinates.x)

library(FNN)
library(dplyr)

# for precipitation ----
# Unique precip stations
prec_stations <- all_prec_data_new_with_elev %>%
  distinct(location_name, coordinates.x, coordinates.y, elevation)

# Extract coordinates and elevation
coords_prec <- prec_stations %>% select(coordinates.x, coordinates.y)
elev_prec <- prec_stations$elevation
n_prec <- nrow(prec_stations)

# Nearest neighbors (k = 7 including itself)
nn_prec <- get.knn(coords_prec, k = 7)
neighbor_indices_prec <- nn_prec$nn.index

# Compute slope for each precipitation station
slope_prec <- numeric(n_prec)

for (i in 1:n_prec) {
  neigh_ids <- neighbor_indices_prec[i, 2:7]
  pairs <- combn(neigh_ids, 2, simplify = FALSE)
  
  slopes <- sapply(pairs, function(pair) {
    h <- pair[1]; l <- pair[2]
    
    if (is.na(elev_prec[h]) || is.na(elev_prec[l])) return(NA_real_)
    
    dx <- coords_prec[h, 1] - coords_prec[l, 1]
    dy <- coords_prec[h, 2] - coords_prec[l, 2]
    dist_hl <- sqrt(dx^2 + dy^2)
    
    if (dist_hl == 0) return(NA_real_)
    
    return(abs(elev_prec[h] - elev_prec[l]) / dist_hl)
  })
  
  # Remove NAs before averaging
  slopes_numeric <- unlist(slopes)
  slopes_numeric <- slopes_numeric[!is.na(slopes_numeric)]
  
  if (length(slopes_numeric) > 0) {
    slope_prec[i] <- mean(slopes_numeric)
  } else {
    slope_prec[i] <- NA_real_
    cat("⚠️ Station", i, "has no valid slopes\n")
  }
}

# Merge back
prec_stations$slope <- slope_prec

# Join into full data
all_prec_data_new_with_slope <- all_prec_data_new_with_elev %>%
  left_join(prec_stations %>% select(location_name, slope), by = "location_name")


# For river flow ----
library(FNN)
library(dplyr)

# Step 1: Extract unique flow stations with coordinates and elevation
q_stations <- all_q_data_new_with_elev %>%
  distinct(location_name, coordinates.x, coordinates.y, elevation)

coords_q <- q_stations %>% select(coordinates.x, coordinates.y)
elev_q <- q_stations$elevation
n_q <- nrow(q_stations)

# Step 2: Find 7 nearest neighbors (including self)
nn_q <- get.knn(coords_q, k = 7)
neighbor_indices_q <- nn_q$nn.index

# Step 3: Compute slope for each flow station
slope_q <- numeric(n_q)

for (i in 1:n_q) {
  neigh_ids <- neighbor_indices_q[i, 2:7]  # Exclude the station itself
  pairs <- combn(neigh_ids, 2, simplify = FALSE)
  
  slopes <- sapply(pairs, function(pair) {
    h <- pair[1]; l <- pair[2]
    
    if (is.na(elev_q[h]) || is.na(elev_q[l])) return(NA_real_)
    
    dx <- coords_q[h, 1] - coords_q[l, 1]
    dy <- coords_q[h, 2] - coords_q[l, 2]
    dist_hl <- sqrt(dx^2 + dy^2)
    
    if (dist_hl == 0) return(NA_real_)
    
    return(abs(elev_q[h] - elev_q[l]) / dist_hl)
  })
  
  slopes_numeric <- unlist(slopes)
  slopes_numeric <- slopes_numeric[!is.na(slopes_numeric)]
  
  if (length(slopes_numeric) > 0) {
    slope_q[i] <- mean(slopes_numeric)
  } else {
    slope_q[i] <- NA_real_
    cat("⚠️ Station", i, "has no valid slope value\n")
  }
}

# Step 4: Merge slope into station metadata
q_stations$slope <- slope_q

# Step 5: Join slope back into full Q dataset
all_q_data_new_with_slope <- all_q_data_new_with_elev %>%
  left_join(q_stations %>% select(location_name, slope), by = "location_name")


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

###############################################################################
# 3) Restructure all_q_data_with_slope by Unique Dates
###############################################################################
y.n.2 <- list()
all_dates <- unique(all_q_data_new_with_slope$data)

for (i in seq_along(all_dates)) {
  current_date <- all_dates[i]
  subset_data <- all_q_data_new_with_slope[all_q_data_new_with_slope$data == current_date, ]
  
  # Include slope in sub_df
  sub_df <- data.frame(
    location_name = subset_data$location_name,
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    slope     = subset_data$slope,      # NEW
    Q             = subset_data$Q,
    data          = subset_data$data, 
    Q_centered    = subset_data$Q_centered
  )
  
  y.n.2[[i]] <- sub_df
}

str(y.n.2[1:5])  # Check structure

###############################################################################
# 4) Restructure all_prec_data_with_slope by Unique Dates
###############################################################################
x.n.2 <- list()
all_dates_prec <- unique(all_prec_data_new_with_slope$data)

for (i in seq_along(all_dates_prec)) {
  current_date <- all_dates_prec[i]
  subset_data <- all_prec_data_new_with_slope[all_prec_data_new_with_slope$data == current_date, ]
  
  # Include slope in sub_df
  sub_df <- data.frame(
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    slope     = subset_data$slope,      # NEW
    prec          = subset_data$prec,
    data          = subset_data$data
  )
  
  x.n.2[[i]] <- sub_df
}

str(x.n.2[1:5])  # Check structure

###############################################################################
# GAM MODEL for X (Precipitation) with slope as an Additional Regressor
###############################################################################
# Define the dependent variable name for precipitation
x.names <- "prec"
q <- length(x.names)
N <- length(x.n.2)  # Total number of dates in the x.n list

# Initialize the list to store models for precipitation
models.plot.x <- vector(mode = "list", length = q)
names(models.plot.x) <- x.names

no.unit.x <- c()
h <- 0

# Ensure the prediction grid 'new.data' has an slope column.
# If not available, assign the mean slope from all_prec_data_with_slope.
if (!("slope" %in% colnames(new.data))) {
  mean_slope <- mean(all_prec_data_with_slope$slope, na.rm = TRUE)
  new.data$slope <- mean_slope
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
    subset_data <- x.n.2[[i]]
    
    # Prepare the data for model fitting:
    # Convert the rslopeant columns to numeric vectors.
    y   <- as.data.frame(subset_data[, x.names[j]])
    x1  <- as.data.frame(subset_data$coordinates.x)
    x2  <- as.data.frame(subset_data$coordinates.y)
    slope <- as.data.frame(subset_data$slope)
    
    data <- data.frame(
      y    = as.numeric(unlist(y)),
      x1   = as.numeric(unlist(x1)),
      x2   = as.numeric(unlist(x2)),
      slope = as.numeric(unlist(slope))
    )
    
    # Check for enough unique values for each regressor
    unique_x1   <- length(unique(data$x1))
    unique_x2   <- length(unique(data$x2))
    unique_slope <- length(unique(data$slope))
    
    # Only proceed if there is variability in all predictors
    if (unique_x1 > 1 && unique_x2 > 1 && unique_slope > 1) { 
      # Set smoothing parameters based on the number of unique values
      k_val  <- 25
      data$y <- log1p(data$y)
      
      # Try fitting the GAM model including slope as an additional smooth term.
      prova <- try(
        gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
      )
      
      if (!inherits(prova, "try-error")) {
        model_i <- gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        
        # Prepare the newdata for predictions.
        # Here, new.data is our grid; ensure it has a numeric slope column.
        new_data_pred <- data.frame(
          x1   = as.numeric(new.data$coordinates.x),
          x2   = as.numeric(new.data$coordinates.y),
          slope = as.numeric(new.data$slope)
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
      rm(unique_x1, unique_x2, unique_slope, k_val)
      
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
  max_Tpoints <- max(sapply(y.n.2, function(subset_data) nrow(subset_data)))  # Dynamic dimension calculation
  
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
max_Tpoints <- max(sapply(y.n.2, function(subset_data) nrow(subset_data)))
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
range(H.Y) # -0.002455248  0.138696242
range(H.X) # 0.0004514201 0.0069412541

# ORTHONORMAL BASIS ----
# Define the number of basis functions
L <- 10  # For Q (river flow)
H <-10  # For Prec (precipitation)

# Generating basis functions for Q (River Flow)
phi.Y <- bs(seq(0, 1, length.out = nrow(H.Y)), df = L, intercept = TRUE)


# Generating basis functions for Prec (Precipitation)
phi.X <- bs(seq(0, 1, length.out = nrow(H.X)), df = H, intercept = TRUE)

# SCORES ----
# --- Initialize Score Arrays ---
gamma_new <- array(NA, dim = c(p, N, L), 
                   dimnames = list(y.names, paste0("N_", seq(1, N)), paste0("L_", seq(1, L))))
chi_new <- array(NA, dim = c(q, N, H), 
                 dimnames = list(x.names, paste0("N_", seq(1, N)), paste0("H_", seq(1, H))))

# --- Calculate Scores for Y (Q) ---
for (j in seq_len(p)) {
  for (i in seq_len(N)) {
    current_weight_y <- weight.y[j, i, ]
    current_Y_star <- Y.star[j, i, ]
    
    # Number of valid points for this day
    num_points <- sum(!is.na(current_weight_y))
    
    if (num_points > 0) {
      # Extract rslopeant portions of the weight and data
      valid_weights <- current_weight_y[1:num_points]
      valid_Y_star <- current_Y_star[1:num_points]
      valid_Y_mean <- Y.mean[j, 1:num_points]
      
      # Perform the basis function decomposition
      tryCatch({
        gamma_new[j, i, ] <- solve(t(phi.Y) %*% diag(valid_weights) %*% phi.Y) %*% 
          t(phi.Y) %*% diag(valid_weights) %*% (valid_Y_star - valid_Y_mean)
      }, error = function(e) {
        gamma_new[j, i, ] <- rep(0, L)  # Fill with zeros if calculation fails
      })
    } else {
      gamma_new[j, i, ] <- rep(0, L)  # If no valid points, assign zeros
    }
  }
}



# --- Calculate Scores for X (Prec) ---
for (j in seq_len(q)) {
  for (i in seq_len(N)) {
    # 1) Extract the day's data
    subset_data <- x.n.2[[i]]
    
    # 2) FALLBACK: if there's no spatial variation, give a uniform first‐basis score
    if (var(subset_data$prec, na.rm = TRUE) == 0) {
      # assign all mass to the first basis function
      chi_new[j, i, ] <- c(1, rep(0, H - 1))
      next  # skip the rest of this loop iteration
    }
    
    # 3) Otherwise proceed with usual weight/X-star extraction
    current_weight_x <- weight.x[j, i, ]
    current_X_star   <- X.star[j, i, ]
    
    tryCatch({
      chi_new[j, i, ] <- solve(
        t(phi.X) %*% diag(current_weight_x) %*% phi.X
      ) %*% t(phi.X) %*% diag(current_weight_x) %*% (current_X_star - X.mean[j, ])
    }, error = function(e) {
      chi_new[j, i, ] <- rep(0, H) 
    })
  }
}


# --- Check Scores ---
which(apply(chi_new[1,,], 1, function(h) sum(h)==0)) # 0
which(apply(gamma_new[1,,], 1, function(h) sum(h)==0)) # 0

# REGRESSION ----
# --- STEP 1: Retrieve dimensions ---
H <- dim(chi_new)[3]
L <- dim(gamma_new)[3]
N <- dim(gamma_new)[2]

# --- STEP 2: Initialize predictions ---
gamma_pred_new <- matrix(NA, nrow = N, ncol = L)
gamma_lag <- gamma[1, dim(gamma)[2], ]  # last gamma from training

# --- STEP 3: Recursive prediction loop ---
for (t in 1:N) {
  chi_t <- chi_new[1, t, ]
  chi_t_minus_1 <- if (t > 1) chi_new[1, t - 1, ] else chi_t
  
  # Replicate slope scores
  slope_t <- slope_scores  # static
  
  # Construct regression input
  X_reg_row <- c(1, chi_t, chi_t_minus_1, gamma_lag, slope_t)
  
  # Predict gamma
  gamma_t_pred <- X_reg_row %*% B_matrix
  gamma_pred_new[t, ] <- gamma_t_pred
  gamma_lag <- gamma_t_pred
}

# --- STEP 4: AR(1) correction ---
residuals_2024 <- matrix(0, nrow = N, ncol = L)
for (l in seq_len(L)) {
  for (t in 3:N) {
    residuals_2024[t, l] <- -phi[l] * (gamma_pred_new[t - 1, l] - gamma_pred_new[t - 2, l])
  }
}
gamma_pred_corrected <- gamma_pred_new + residuals_2024

# --- STEP 5: Reconstruct Q (still log1p scale) ---
Q_pred_log <- gamma_pred_corrected %*% t(phi.Y) +
  matrix(rep(Y.mean, each = N), nrow = N, byrow = FALSE)

# --- STEP 6: Re-add station means and back-transform ---
station_means <- all_q_data_imputed %>%
  mutate(logQ = log1p(Q)) %>%
  group_by(location_name) %>%
  summarize(mean_logQ = mean(logQ, na.rm = TRUE)) %>%
  arrange(location_name)

mean_logQ_vec <- station_means$mean_logQ

Q_pred_real_ar1 <- exp(Q_pred_log + matrix(mean_logQ_vec, nrow = nrow(Q_pred_log), ncol = ncol(Q_pred_log), byrow = TRUE)) - 1

# --- STEP 7: Construct list format (one per date) ---
# Get the corresponding dates (assuming x.n.2 is the forecast input list)
forecast_dates <- sapply(x.n.2, function(df) unique(df$data))

# Get corresponding station coordinates
locations_q <- all_q_data_imputed %>%
  distinct(location_name, coordinates.x, coordinates.y) %>%
  arrange(location_name)

# Rebuild final prediction list
results_df_ar1 <- vector("list", length = N)
for (i in seq_len(N)) {
  df_i <- data.frame(
    date = forecast_dates[i],
    location_name = locations_q$location_name,
    coordinates.x = locations_q$coordinates.x,
    coordinates.y = locations_q$coordinates.y,
    Q_pred = Q_pred_real_ar1[i, ]
  )
  results_df_ar1[[i]] <- df_i
}
results_df_ar1 <- bind_rows(results_df_ar1)


# 1. Reconstruct the final predicted Q with AR(1) correction (already done)
# `Q_pred_real_ar1` is [N x Tpoints] matrix

# 2. Create a list to hold all prediction results per day
results_list_ar1 <- list()

# Loop over each day in y.n.2 (which matches all_q_data_new_imputed$data)
for (i in seq_len(length(y.n.2))) {
  sub_data <- y.n.2[[i]]
  n_obs <- nrow(sub_data)
  
  df_day <- data.frame(
    date           = as.Date(sub_data$data),
    location_name  = sub_data$location_name,
    coordinates.x  = sub_data$coordinates.x,
    coordinates.y  = sub_data$coordinates.y,
    Q_true         = sub_data$Q,  # Real Q from the original dataset
    Q_pred         = Q_pred_real_ar1[i, 1:n_obs]  # AR(1)-corrected prediction
  )
  
  results_list_ar1[[i]] <- df_day
}

# 3. Combine into a single dataframe
results_df_ar1 <- do.call(rbind, results_list_ar1)

# 4. Check structure
head(results_df_ar1)

results_df_ar1 %>% 
  group_by(location_name) %>% 
  summarize(RMSE = sqrt(mean((Q_true - Q_pred)^2, na.rm = TRUE))) %>% 
  print(n=38)


# PLOTS ----
# Get unique locations
location_names <- unique(results_df_ar1$location_name)

# Create an empty list to store all plots
plot_list <- list()

# Loop over each location and create ggplot
for (loc in location_names) {
  df_loc <- results_df_ar1 %>% filter(location_name == loc)
  
  p <- ggplot(df_loc, aes(x = date)) +
    geom_line(aes(y = Q_true, color = "True Q")) +
    geom_line(aes(y = Q_pred, color = "Predicted Q"))+
    labs(title = loc, x = "Date", y = "River Flow (Q)") +
    scale_color_manual(values = c("True Q" = "black", "Predicted Q" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  plot_list[[loc]] <- p
}

print(plot_list)



ggplot(results_df_ar1) + geom_point(aes(x = Q_true, y = Q_pred)) + geom_abline()


# CROSS-VALIDATION ----
# CROSS VALIDATION ----
#############################
# TRAINING PIPELINE FUNCTION
#############################
train_pipeline <- function(q_train_data, prec_train_data) {
  
  ### 1. Restructure Q training data by date (to create list y.n)
  sorted_dates_q <- sort(unique(q_train_data$data))
  y.n <- lapply(sorted_dates_q, function(dt) {
    q_train_data %>% 
      filter(data == dt) %>% 
      select(location_name, Q, data, coordinates.x, coordinates.y, slope)
  })
  N_q <- length(y.n)
  
  ### 2. Restructure Precipitation training data by date (to create list x.n)
  sorted_dates_prec <- sort(unique(prec_train_data$data))
  x.n <- lapply(sorted_dates_prec, function(dt) {
    prec_train_data %>% 
      filter(data == dt) %>% 
      select(coordinates.x, coordinates.y, prec, data, slope)
  })
  N_prec <- length(x.n)
  
  ### 3. Create a spatial grid (new.data) from Q training coordinates 
  dim.grid <- 25
  a <- dim.grid; b <- dim.grid
  lat_range <- range(q_train_data$coordinates.y)
  long_range <- range(q_train_data$coordinates.x)
  new.data <- matrix(NA, a * b, 2)
  new.data[,1] <- sort(rep(seq(lat_range[1], lat_range[2], length.out = a), b))
  new.data[,2] <- rep(seq(long_range[1], long_range[2], length.out = b), a)
  new.data <- as.data.frame(new.data)
  colnames(new.data) <- c("coordinates.y", "coordinates.x")
  N_T <- nrow(new.data)
  
  ### 4. GAM MODEL FOR Q (Y) [No time effect here; it is added later]
  y.names <- "Q"
  p <- 1
  max_Tpoints_q <- max(sapply(y.n, nrow))
  weight.y <- Y.star <- array(NA, dim = c(p, N_q, max_Tpoints_q),
                              dimnames = list(y.names, paste0("N_", 1:N_q), paste0("T_", 1:max_Tpoints_q)))
  models.plot <- vector("list", p)
  names(models.plot) <- y.names
  
  for(j in seq_len(p)) {
    models.plot[[y.names[j]]] <- list()
    for(i in seq_len(N_q)) {
      dat <- y.n[[i]]
      y_val <- as.data.frame(dat[, "Q"])
      x1 <- as.data.frame(dat$coordinates.x)
      x2 <- as.data.frame(dat$coordinates.y)
      data_temp <- data.frame(
        y = y_val, 
        x1 = as.numeric(unlist(x1)), 
        x2 = as.numeric(unlist(x2))
      )
      colnames(data_temp) <- c("y", "x1", "x2")
      
      unique_x1 <- length(unique(data_temp$x1))
      unique_x2 <- length(unique(data_temp$x2))
      if(unique_x1 > 1 && unique_x2 > 1) {
        k_val <- 25
        data_temp$y <- log1p(data_temp$y)
        model_i <- gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian())
        models.plot[[y.names[j]]][[i]] <- list(model = model_i)
        new_data_pred <- data.frame(x1 = data_temp$x1, x2 = data_temp$x2)
        preds <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        w_ji <- 1 / ((preds$se.fit)^2)
        w_ji[is.na(w_ji) | is.infinite(w_ji)] <- 0
        num_points <- length(preds$fit)
        weight.y[y.names[j], i, 1:num_points] <- w_ji
        Y.star[y.names[j], i, 1:num_points] <- preds$fit
      } else {
        weight.y[y.names[j], i, ] <- NA
        Y.star[y.names[j], i, ] <- NA
      }
    }
  }
  # Compute the mean Q prediction (Y.mean) across training dates:
  Y.mean <- matrix(0, p, max_Tpoints_q)
  rownames(Y.mean) <- y.names
  colnames(Y.mean) <- paste0("T_", 1:max_Tpoints_q)
  for(j in seq_len(p)) {
    for(i in seq_len(N_q)) {
      num_points <- sum(!is.na(weight.y[y.names[j], i, ]))
      if(num_points > 0) {
        Y.mean[y.names[j], 1:num_points] <- Y.mean[y.names[j], 1:num_points] +
          weight.y[y.names[j], i, 1:num_points] * Y.star[y.names[j], i, 1:num_points]
      }
    }
    weight_sum <- apply(weight.y[y.names[j], , ], 2, sum, na.rm = TRUE)
    weight_sum[weight_sum == 0] <- NA
    Y.mean[y.names[j], ] <- Y.mean[y.names[j], ] / weight_sum
  }
  
  ### 5. GAM MODEL FOR PRECIPITATION (X)
  x.names <- "prec"
  q_val <- 1
  max_Tpoints_x <- N_T  # Using grid size for precipitation predictions.
  weight.x <- X.star <- array(NA, dim = c(q_val, N_prec, max_Tpoints_x),
                              dimnames = list(x.names, paste0("N_", 1:N_prec), paste0("T_", 1:max_Tpoints_x)))
  models.plot.x <- vector("list", q_val)
  names(models.plot.x) <- x.names
  
  for(j in seq_len(q_val)) {
    models.plot.x[[x.names[j]]] <- list()
    for(i in seq_len(N_prec)) {
      dat <- x.n[[i]]
      y_val <- as.data.frame(dat[, "prec"])
      x1 <- as.data.frame(dat$coordinates.x)
      x2 <- as.data.frame(dat$coordinates.y)
      data_temp <- data.frame(
        y = y_val, 
        x1 = as.numeric(unlist(x1)), 
        x2 = as.numeric(unlist(x2))
      )
      colnames(data_temp) <- c("y", "x1", "x2")
      
      unique_x1 <- length(unique(data_temp$x1))
      unique_x2 <- length(unique(data_temp$x2))
      if(unique_x1 > 1 && unique_x2 > 1) {
        k_val <- 25
        data_temp$y <- log1p(data_temp$y)
        model_i <- gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian())
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        # For precipitation predictions on the spatial grid, use a representative time value.
        median_time <- median(as.numeric(prec_train_data$data))
        new_data_pred <- data.frame(x1 = as.numeric(new.data$coordinates.x),
                                    x2 = as.numeric(new.data$coordinates.y),
                                    time = median_time)
        preds <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        w_ji <- 1 / ((preds$se.fit)^2)
        w_ji[is.na(w_ji) | is.infinite(w_ji)] <- 0
        weight.x[x.names[j], i, ] <- w_ji
        X.star[x.names[j], i, ] <- preds$fit
      } else {
        weight.x[x.names[j], i, ] <- NA
        X.star[x.names[j], i, ] <- NA
      }
    }
  }
  
  ### 5b. GAM MODEL FOR SLOPE ----
  models.plot.slope <- list()
  if (!("slope" %in% colnames(new.data))) {
    new.data$slope <- mean(prec_train_data$slope, na.rm = TRUE)
  }
  weight.slope <- slope.star <- array(
    NA, dim = c(N_prec, N_T), 
    dimnames = list(paste0("N_", 1:N_prec), paste0("T_", 1:N_T))
  )
  
  for (i in seq_len(N_prec)) {
    dat <- x.n[[i]]
    if (!("slope" %in% names(dat))) next
    
    y  <- as.numeric(dat$slope)
    x1 <- as.numeric(dat$coordinates.x)
    x2 <- as.numeric(dat$coordinates.y)
    data_temp <- data.frame(y = y, x1 = x1, x2 = x2)
    
    if (length(unique(x1)) > 1 && length(unique(x2)) > 1 && length(unique(y)) > 1) {
      model_i <- gam(y ~ s(x1, x2, k = 25), data = data_temp)
      models.plot.slope[[i]] <- list(model = model_i)
      preds <- predict(model_i, newdata = data.frame(
        x1 = new.data$coordinates.x,
        x2 = new.data$coordinates.y
      ), se.fit = TRUE)
      
      w <- 1 / preds$se.fit^2
      w[is.na(w) | is.infinite(w)] <- 0
      
      weight.slope[i, ] <- w
      slope.star[i, ] <- preds$fit
    }
  }
  
  # Compute the mean precipitation prediction (X.mean) across dates:
  X.mean <- matrix(0, q_val, max_Tpoints_x)
  rownames(X.mean) <- x.names
  colnames(X.mean) <- paste0("T_", 1:max_Tpoints_x)
  for(j in seq_len(q_val)) {
    for(i in seq_len(N_prec)) {
      num_points <- sum(!is.na(weight.x[x.names[j], i, ]))
      if(num_points > 0) {
        X.mean[x.names[j], 1:num_points] <- X.mean[x.names[j], 1:num_points] +
          weight.x[x.names[j], i, 1:num_points] * X.star[x.names[j], i, 1:num_points]
      }
    }
    weight_sum <- apply(weight.x[x.names[j], , ], 2, sum)
    X.mean[x.names[j], ] <- X.mean[x.names[j], ] / weight_sum
  }
  
  ### 6. TWO BASIS SYSTEMS: Covariance matrices for Q and precipitation
  sigma.Y_j <- vector("list", p)
  names(sigma.Y_j) <- y.names
  for(j in seq_len(p)) {
    sigma.Y_ji_num <- sigma.Y_ji_den <- matrix(0, max_Tpoints_q, max_Tpoints_q)
    for(i in seq_len(N_q)) {
      curr_w <- weight.y[y.names[j], i, ]
      curr_Y <- Y.star[y.names[j], i, ]
      valid_indices <- which(!is.na(curr_w))
      if(length(valid_indices) > 0) {
        valid_w <- curr_w[valid_indices]
        valid_Y <- curr_Y[valid_indices]
        valid_Y_mean <- Y.mean[y.names[j], valid_indices]
        A <- matrix(valid_w * (valid_Y - valid_Y_mean), nrow = 1)
        temp_num <- t(A) %*% A  
        temp_den <- outer(valid_w, valid_w)
        embed_num <- matrix(0, max_Tpoints_q, max_Tpoints_q)
        embed_den <- matrix(0, max_Tpoints_q, max_Tpoints_q)
        embed_num[1:length(valid_indices), 1:length(valid_indices)] <- temp_num
        embed_den[1:length(valid_indices), 1:length(valid_indices)] <- temp_den
        sigma.Y_ji_num <- sigma.Y_ji_num + embed_num
        sigma.Y_ji_den <- sigma.Y_ji_den + embed_den
      }
    }
    sigma.Y_j[[y.names[j]]] <- sigma.Y_ji_num / sigma.Y_ji_den
  }
  
  sigma.X_j <- vector("list", q_val)
  names(sigma.X_j) <- x.names
  for(j in seq_len(q_val)) {
    sigma.X_ji_num <- sigma.X_ji_den <- matrix(0, max_Tpoints_x, max_Tpoints_x)
    for(i in seq_len(N_prec)) {
      curr_w <- weight.x[x.names[j], i, ]
      curr_X <- X.star[x.names[j], i, ]
      valid_indices <- which(!is.na(curr_w))
      if(length(valid_indices) > 0) {
        valid_w <- curr_w[valid_indices]
        valid_X <- curr_X[valid_indices]
        valid_X_mean <- X.mean[x.names[j], valid_indices]
        A <- matrix(valid_w * (valid_X - valid_X_mean), nrow = 1)
        temp_num <- t(A) %*% A
        temp_den <- outer(valid_w, valid_w)
        embed_num <- matrix(0, max_Tpoints_x, max_Tpoints_x)
        embed_den <- matrix(0, max_Tpoints_x, max_Tpoints_x)
        embed_num[1:length(valid_indices), 1:length(valid_indices)] <- temp_num
        embed_den[1:length(valid_indices), 1:length(valid_indices)] <- temp_den
        sigma.X_ji_num <- sigma.X_ji_num + embed_num
        sigma.X_ji_den <- sigma.X_ji_den + embed_den
      }
    }
    sigma.X_j[[x.names[j]]] <- sigma.X_ji_num / sigma.X_ji_den
  }
  
  # Mean covariance matrices:
  H.Y_mat <- matrix(0, max_Tpoints_q, max_Tpoints_q)
  for(j in seq_len(p)) {
    H.Y_mat <- H.Y_mat + sigma.Y_j[[y.names[j]]]
  }
  H.Y_mat <- H.Y_mat / p
  
  H.X_mat <- matrix(0, max_Tpoints_x, max_Tpoints_x)
  for(j in seq_len(q_val)) {
    H.X_mat <- H.X_mat + (1 / q_val) * sigma.X_j[[x.names[j]]]
  }
  
  ### 7. Orthonormal Bases (using splines)
  L <- 10       # Number of basis functions for Q
  H_basis <- 10 # Number of basis functions for precipitation
  phi.Y <- bs(seq(0, 1, length.out = nrow(H.Y_mat)), df = L, intercept = TRUE)
  phi.X <- bs(seq(0, 1, length.out = nrow(H.X_mat)), df = H_basis, intercept = TRUE)
  
  # Project slope onto phi.X basis
  row_idx <- which(rowSums(!is.na(slope.star)) > 0)[1]
  slope_surface <- slope.star[row_idx, ]
  slope_centered <- slope_surface - mean(slope_surface, na.rm = TRUE)
  slope_centered[is.na(slope_centered)] <- 0
  slope_scores <- as.vector(t(phi.X) %*% slope_centered)
  
  # Create static predictor matrix replicated for each Q observation
  x_reg_slope <- matrix(rep(slope_scores, each = N_q), nrow = N_q)
  colnames(x_reg_slope) <- paste0("slope_basis_", 1:H_basis)
  
  
  ### 8. SCORES: Compute gamma (for Q) and chi (for precipitation) ###
  gamma_arr <- array(NA, dim = c(p, N_q, L),
                     dimnames = list(y.names, paste0("N_", 1:N_q), paste0("L_", 1:L)))
  for(j in seq_len(p)) {
    for(i in seq_len(N_q)) {
      curr_w <- weight.y[y.names[j], i, ]
      curr_Y <- Y.star[y.names[j], i, ]
      valid_indices <- which(!is.na(curr_w))
      if(length(valid_indices) > 0) {
        gamma_arr[j, i, ] <- tryCatch({
          solve(t(phi.Y) %*% diag(curr_w[valid_indices]) %*% phi.Y) %*%
            t(phi.Y) %*% diag(curr_w[valid_indices]) %*% (curr_Y[valid_indices] - Y.mean[y.names[j], valid_indices])
        }, error = function(e) { rep(0, L) })
      } else {
        gamma_arr[j, i, ] <- rep(0, L)
      }
    }
  }
  
  chi_arr <- array(NA, dim = c(q_val, N_prec, H_basis),
                   dimnames = list(x.names, paste0("N_", 1:N_prec), paste0("H_", 1:H_basis)))
  for(j in seq_len(q_val)) {
    for(i in seq_len(N_prec)) {
      curr_w <- weight.x[x.names[j], i, ]
      curr_X <- X.star[x.names[j], i, ]
      valid_indices <- which(!is.na(curr_w))
      if(length(valid_indices) > 0 && sum(curr_w[valid_indices], na.rm = TRUE) != 0) {
        chi_arr[j, i, ] <- tryCatch({
          solve(t(phi.X) %*% diag(curr_w[valid_indices]) %*% phi.X) %*%
            t(phi.X) %*% diag(curr_w[valid_indices]) %*% (curr_X[valid_indices] - X.mean[x.names[j], valid_indices])
        }, error = function(e) { rep(0, H_basis) })
      } else {
        chi_arr[j, i, ] <- rep(0, H_basis)
      }
    }
  }
  
  ### 9. REGRESSION: Construct regression matrices and fit cglasso 
  # --- Build predictor matrix for regression ---
  # (i) Precipitation block from chi scores; using the first observation for each Q date.
  x_reg_chi <- matrix(NA, N_q, q_val * H_basis)
  colnames(x_reg_chi) <- paste0("h.", sort(rep(1:H_basis, q_val)), "_", rep(x.names, H_basis))
  for(h in seq_len(H_basis)) {
    x_reg_chi[, h] <- chi_arr[1, , h]
  }
  
  
  
  # Combine all predictor blocks:
  x_reg <- cbind(x_reg_chi, x_reg_slope)
  
  
  # --- Build response matrix from gamma scores (using first observation per Q date) ---
  y_reg <- matrix(NA, N_q, p * L)
  colnames(y_reg) <- paste0("l.", sort(rep(1:L, p)), "_", rep(y.names, L))
  for(l in seq_len(L)) {
    y_reg[, l] <- gamma_arr[1, , l]
  }
  
  # Package data for cglasso and fit the penalized regression.
  data_cg <- datacggm(Y = y_reg, X = x_reg)
  model_cg <- cglasso(. ~ ., data = data_cg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)
  
  # Extract the coefficient vector.
  B_vector <- model_cg$B[-1, , 2, 1]
  chi_coeff_count <- q_val * H_basis
  slope_coeff_count <- H_basis
  B_matrix <- matrix(B_vector[1:chi_coeff_count], nrow = H_basis, ncol = L)
  B_slope <- B_vector[(chi_coeff_count + 1):(chi_coeff_count + slope_coeff_count)]
  
  
  cat("Precipitation Coefficient Matrix (B_matrix):\n")
  print(B_matrix)
  cat("Estimated Slope Coefficients (B_slope):\n")
  print(B_slope)
  
  
  #############################
  # Return all relevant objects
  #############################
  return(list(
    B_matrix = B_matrix,
    phi.Y = phi.Y,
    phi.X = phi.X,
    Y.mean = Y.mean,
    models.plot = models.plot,      # GAM models for Q (by date)
    models.plot.x = models.plot.x,    # GAM models for prec (by date)
    new.data = new.data,
    max_Tpoints_q = max_Tpoints_q,
    max_Tpoints_x = max_Tpoints_x,
    L = L,
    H_basis = H_basis,
    q_train_dates = sorted_dates_q,   # sorted training dates for Q
    B_slope = B_slope,
    slope_scores = slope_scores
    
  ))
}

#############################
### PREDICTION FUNCTION (with slope)
#############################
predict_pipeline <- function(test_data, model_trained) {
  B_slope <- model_trained$B_slope
  slope_scores <- model_trained$slope_scores
  L_slope <- length(slope_scores)
  
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  
  for (i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data
    idx <- match(test_date, model_trained$q_train_dates)
    
    # --- GAM prediction for Q (spatial only)
    if (is.na(idx)) {
      gam_pred <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(
        x1 = test_row$coordinates.x, 
        x2 = test_row$coordinates.y
      )
      gam_pred <- predict(gam_model_Y, newdata = test_coords, type = "response")
    }
    
    # --- Static slope adjustment (projected slope_scores ⋅ B_slope)
    slope_adjustment <- sum(slope_scores * B_slope)
    
    # Combine GAM and slope-adjusted component
    preds[i] <- gam_pred + slope_adjustment
  }
  
  preds <- exp(preds) - 1  # Back-transform from log1p
  return(preds)
}


### CROSS-VALIDATION SIMULATION
rmse_results <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_with_slope$location_name)
set.seed(1234)  # For reproducibility

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_with_slope$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_data <- all_q_data_with_slope[-test_idx, ]
  q_test_data  <- all_q_data_with_slope[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_data, all_prec_data_with_slope)
  
  # Predict on the test set using the prediction function.
  predictions <- predict_pipeline(q_test_data, model_trained)
  
  actual <- q_test_data$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results <- rbind(rmse_results, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
}
cat("\n--- RMSE Results for each location ---\n")
print(rmse_results)
