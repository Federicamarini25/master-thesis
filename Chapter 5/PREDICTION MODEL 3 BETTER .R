rm(list = ls(all = TRUE))     
cat("\014")                    

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

# OPTION 1: Integrate Elevation into the 3D Smooth  ----

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
# TRY SOAP BELL / BALL bs = ... 
# if I put the 3 of them all together


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
# STEP 1: Construct current-day chi scores (chi_t)
x_now <- matrix(NA, N, q * H)
colnames(x_now) <- paste0("h.", sort(rep(seq_len(H), q)), "_", rep(x.names, H))
for (h in seq_len(H)) {
  x_now[, h] <- chi[1, , h]
}

# STEP 2: Construct lagged gamma scores (gamma_{t-1})
gamma_lag <- matrix(NA, N, p * L)
colnames(gamma_lag) <- paste0("lag.l.", sort(rep(seq_len(L), p)), "_", rep(y.names, L))
for (l in seq_len(L)) {
  gamma_lag[, l] <- gamma[1, , l]
}
gamma_lag <- rbind(rep(NA, p * L), gamma_lag[-N, ])  # shift down by one row

# STEP 3: Build the full design matrix x (NO lagged chi anymore)
x <- cbind(x_now, gamma_lag)
x <- x[-1, ]  # drop first row due to lagged NA
x <- cbind(1, x)  # add intercept

# STEP 4: Build response matrix y
y <- matrix(NA, N, p * L)
colnames(y) <- paste0("l.", sort(rep(seq_len(L), p)), "_", rep(y.names, L))
for (l in seq_len(L)) {
  y[, l] <- gamma[1, , l]
}
y <- y[-1, ]  # align with x

# STEP 5: Fit the penalized regression model
library(cglasso)
data <- datacggm(Y = y, X = x)
model <- cglasso(. ~ ., data = data, nlambda = 2, lambda.min.ratio = 0, nrho = 1)

# STEP 6: Extract and reshape the coefficient matrix
B_vector <- model$B[-1, , 2, 1]
B_matrix <- matrix(B_vector, nrow = ncol(x), ncol = ncol(y))

# STEP 7: Predict gamma and compute residuals
gamma_pred <- x %*% B_matrix
residuals <- y - gamma_pred #

# STEP 8: Compute MSE
MSE <- mean(residuals^2, na.rm = TRUE)
cat("Mean Squared Error (MSE):", MSE, "\n") # 0.03075721 

# STEP 9: Plotting
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
cat("RMSE from AR(1) innovations:", round(RMSE_AR1, 5), "\n") # 0.02956 

# Optional: summary of AR(1) coefficients
summary(phi)



# NEW ----
# Fetch new data
domain <- "surfacewater"
parameter <- "Q"
from_date_new <- as.Date("2024-01-01")
to_date_new <- as.Date("2024-08-31")

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

# DO EVERYTHING AGAIN FOR MY NEW DATA ----
lat_range <- range(all_q_data_new_imputed$coordinates.y)
long_range <- range(all_q_data_new_imputed$coordinates.x)

dim.grid <- 25
a <- dim.grid
b <- dim.grid
N_T <- a * b

new.data <- matrix(NA, a * b, 2,
                   dimnames = list("obs" = paste0("obs", seq(1, a*b)),
                                   "coord" = c("coordinates.y", "coordinates.x")))
new.data[, 1] <- sort(rep(seq(lat_range[1], lat_range[2], length.out = a), b))
new.data[, 2] <- rep(seq(long_range[1], long_range[2], length.out = b), a)
new.data <- as.data.frame(new.data)

colnames(new.data) <- c("coordinates.y", "coordinates.x")
N_T <- nrow(new.data)


###############################################################################
# 3) Restructure all_q_data_with_elev by Unique Dates
###############################################################################
y.n.2 <- list()
all_dates <- unique(all_q_data_new_with_elev$data)

for (i in seq_along(all_dates)) {
  current_date <- all_dates[i]
  subset_data <- all_q_data_new_with_elev[all_q_data_new_with_elev$data == current_date, ]
  
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
  
  y.n.2[[i]] <- sub_df
}

str(y.n.2[1:5])  # Check structure

###############################################################################
# 4) Restructure all_prec_data_with_elev by Unique Dates
###############################################################################
x.n.2 <- list()
all_dates_prec <- unique(all_prec_data_new_with_elev$data)

for (i in seq_along(all_dates_prec)) {
  current_date <- all_dates_prec[i]
  subset_data <- all_prec_data_new_with_elev[all_prec_data_new_with_elev$data == current_date, ]
  
  # Include elevation in sub_df
  sub_df <- data.frame(
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    elevation     = subset_data$elevation,      # NEW
    prec          = subset_data$prec,
    data          = subset_data$data
  )
  
  x.n.2[[i]] <- sub_df
}

str(x.n.2[1:5])  # Check structure

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
N <- length(x.n.2)  # Total number of dates in the x.n list

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
    subset_data <- x.n.2[[i]]
    
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
# TRY SOAP BELL / BALL bs = ... 
# if I put the 3 of them all together


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
range(H.Y) # 0.1296407 3.6871226
range(H.X) # -0.0000419009  0.0044038114

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
      # Extract relevant portions of the weight and data
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
    
    # 3) Otherwise proceed with your usual weight/X-star extraction
    current_weight_x <- weight.x[j, i, ]
    current_X_star   <- X.star[j, i, ]
    
    # 4) Now do the solve() … code you already have
    tryCatch({
      chi_new[j, i, ] <- solve(
        t(phi.X) %*% diag(current_weight_x) %*% phi.X
      ) %*% t(phi.X) %*% diag(current_weight_x) %*% (current_X_star - X.mean[j, ])
    }, error = function(e) {
      chi_new[j, i, ] <- rep(0, H)  # optional: keep your original error‐fallback
    })
  }
}


# --- Check Scores ---
which(apply(chi_new[1,,], 1, function(h) sum(h)==0)) # 0
which(apply(gamma_new[1,,], 1, function(h) sum(h)==0)) # 0

# REGRESSION ----
# Retrieve dimensions
H <- dim(chi_new)[3]
L <- dim(gamma_new)[3]
N <- dim(gamma_new)[2]  # Number of days in 2024

# Initialize prediction matrix
gamma_pred_new <- matrix(NA, nrow = N, ncol = L)

# Start with the last gamma from training (2023)
gamma_lag <- gamma[1, dim(gamma)[2], ]

# Loop over each day in 2024
for (t in 1:N) {
  # chi_t: current-day precipitation scores
  chi_t <- chi_new[1, t, ]
  
  # Build design vector for day t (NO chi_{t-1}, only chi_t + gamma_{t-1})
  X_reg_row <- c(
    1,         # intercept
    chi_t,     # current chi
    gamma_lag  # lagged gamma
  )
  
  # Predict gamma_t
  gamma_t_pred <- X_reg_row %*% B_matrix
  gamma_pred_new[t, ] <- gamma_t_pred
  
  # Update lag for next iteration
  gamma_lag <- gamma_t_pred
}

# Step 1: AR(1) residual correction (optional but recommended)
residuals_2024 <- matrix(0, nrow = N, ncol = L)
for (l in seq_len(L)) {
  for (t in 3:N) {
    residuals_2024[t, l] <- -phi[l] * (gamma_pred_new[t - 1, l] - gamma_pred_new[t - 2, l])
  }
}

# Add AR(1)-based correction to predictions
gamma_pred_corrected <- gamma_pred_new + residuals_2024

# Step 2: Reconstruct predicted Q (on log1p scale)
Q_pred_log <- gamma_pred_corrected %*% t(phi.Y) +
  matrix(rep(Y.mean, each = N), nrow = N, byrow = FALSE)

# Step 3: Add station-specific means and back-transform
station_means <- all_q_data_imputed %>%
  mutate(logQ = log1p(Q)) %>%
  group_by(location_name) %>%
  summarize(mean_logQ = mean(logQ, na.rm = TRUE)) %>%
  arrange(location_name)

mean_logQ_vec <- station_means$mean_logQ

Q_pred_real_ar1 <- exp(Q_pred_log + matrix(mean_logQ_vec, nrow = nrow(Q_pred_log), ncol = ncol(Q_pred_log), byrow = TRUE)) - 1


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

# 4. Optional: Check structure
head(results_df_ar1)






# PLOTS ----
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggpubr)


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

results_df_ar1 %>% 
  group_by(location_name) %>% 
  summarize(RMSE = sqrt(mean((Q_true - Q_pred)^2, na.rm = TRUE)))

# PLOT WITH TRAINING DATA (3 YEARS) AND PREDICTIONS ---
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggpubr)

# 1. Filter training data (2020–2023)
training_df <- all_q_data_with_elev %>%
  filter(data > as.Date("2023-01-01")) %>%
  rename(date = data, Q_true = Q) %>%
  mutate(type = "Training", Q_pred = NA)  # No predictions in training period

# 2. Prepare prediction data (2024)
prediction_df <- results_df_ar1 %>%
  mutate(type = "Prediction")

# 3. Combine the two datasets
combined_df <- bind_rows(training_df, prediction_df)

# 4. Plot loop for each location
location_names <- unique(combined_df$location_name)
plot_list <- list()

for (loc in location_names) {
  df_loc <- combined_df %>% filter(location_name == loc)
  
  p <- ggplot(df_loc, aes(x = date)) +
    geom_line(aes(y = Q_true, color = "True Q"), linewidth = 0.8) +
    geom_line(aes(y = Q_pred, color = "Predicted Q"), linewidth = 0.8, na.rm = TRUE) +
    labs(title = paste("Station:", loc), x = "Date", y = "River Flow (Q)") +
    scale_color_manual(values = c("True Q" = "black", "Predicted Q" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  plot_list[[loc]] <- p
}

# 5. Display plots
print(plot_list)


# CORRECTED VS NON CORRECTED PREDICTIONS
# Step 1: Compute non-corrected Q predictions again (just to be sure)
Q_pred_log_nocorr <- gamma_pred_new %*% t(phi.Y) +
  matrix(rep(Y.mean, each = N), nrow = N, byrow = FALSE)

Q_pred_real_nocorr <- exp(Q_pred_log_nocorr) - 1  # [N x Tpoints]

# Step 2: Flatten into a single vector (same order as rows in results_df_ar1)
Q_pred_nocorr_vector <- numeric()

for (i in seq_len(length(y.n.2))) {
  sub_data <- y.n.2[[i]]
  n_obs <- nrow(sub_data)
  
  Q_pred_nocorr_vector <- c(Q_pred_nocorr_vector, Q_pred_real_nocorr[i, 1:n_obs])
}

# Step 3: Add column to existing dataframe
results_df_ar1$Q_pred_nocorr <- Q_pred_nocorr_vector

# Compare both predictions
library(ggplot2)

for (loc in unique(results_df_ar1$location_name)) {
  df_loc <- results_df_ar1 %>% filter(location_name == loc)
  
  p <- ggplot(df_loc, aes(x = date)) +
    geom_line(aes(y = Q_true, color = "True Q")) +
    geom_line(aes(y = Q_pred, color = "Predicted Q (AR1)")) +
    geom_line(aes(y = Q_pred_nocorr, color = "Predicted Q (No AR1)")) +
    labs(title = loc, x = "Date", y = "River Flow (Q)") +
    scale_color_manual(values = c(
      "True Q" = "black",
      "Predicted Q (AR1)" = "blue",
      "Predicted Q (No AR1)" = "red"
    )) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
}


# Compute RMSE per location and sort in descending order
rmse_by_location <- results_df_ar1 %>%
  group_by(location_name) %>%
  summarize(RMSE = sqrt(mean((Q_true - Q_pred)^2, na.rm = TRUE))) %>%
  arrange(desc(RMSE))

print(rmse_by_location)
