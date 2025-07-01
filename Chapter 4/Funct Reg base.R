rm(list = ls(all = TRUE))     
cat("\014")   

# Load required libraries ---
library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(imputeTS)

library(sf)
library(leaflet)
library(osmdata)

library(ggplot2)
library(mgcv)
library(splines)
library(cglasso)

# Setting working directory
working_dir = "/Users/federicamarini/Desktop/Master Thesis"            
setwd(working_dir)    

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
      
      # ✅ NEW: Track missing values at fetch time
      data_df$was_missing <- is.na(data_df[[parameter]])
      
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


# Extract ranges from the datasets
lat_range <- range(all_q_data_imputed$coordinates.y)
long_range <- range(all_q_data_imputed$coordinates.x)

# Map of rivers ----
# 1. ---- Prepare Q locations ----
# Unique river flow station locations
unique_q_locations <- all_q_data_imputed %>%
  select(location_name, coordinates.x, coordinates.y) %>%
  distinct() %>%
  filter(!is.na(coordinates.x), !is.na(coordinates.y), coordinates.x != 0, coordinates.y != 0)

# Convert to sf and transform to WGS84
q_sf <- st_as_sf(unique_q_locations, coords = c("coordinates.x", "coordinates.y"), crs = 2056)
q_wgs84 <- st_transform(q_sf, crs = 4326)

# Extract lon/lat for leaflet
coords <- st_coordinates(q_wgs84)
q_wgs84$longitude <- coords[, 1]
q_wgs84$latitude <- coords[, 2]



# 2. ---- Get rivers from OpenStreetMap ----

# Bounding box for Ticino, Switzerland
bbox_ticino <- getbb("Ticino, Switzerland")

# Query rivers and streams
rivers_osm <- opq(bbox = bbox_ticino) %>%
  add_osm_feature(key = "waterway", value = c("river", "stream")) %>%
  osmdata_sf()

# Get river lines and convert to WGS84
river_lines <- st_transform(rivers_osm$osm_lines, crs = 4326)

# 3. ---- Create leaflet map ----

leaflet() %>%
  addTiles() %>%
  # Add river lines in blue
  addPolylines(data = river_lines, color = "blue", weight = 1, opacity = 0.6) %>%
  # Add Q stations as black dots
  addCircleMarkers(
    data = q_wgs84,
    lng = ~longitude, lat = ~latitude,
    popup = ~location_name,
    label = ~location_name,
    color = "black",
    radius = 5,
    stroke = FALSE,
    fillOpacity = 0.8
  )


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

# Plotting original points (river flow) and new grid points
ggplot() +
  geom_point(data = all_q_data_imputed, aes(x = coordinates.x, y = coordinates.y), color = "black", size = 0.7) +
  geom_point(data = new.data, aes(x = coordinates.x, y = coordinates.y), color = "red", size = 1.5) +  
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()

# RESTRUCTURE ALL_Q_DATA_IMPUTED ----
# Create a list to store sub-datasets
y.n <- list()

# Get all unique dates from the data
all_dates <- unique(all_q_data_imputed$data)

# Loop over each date to create sub-datasets
for (i in seq_along(all_dates)) {
  current_date <- all_dates[i]
  
  # Extract data for the current date
  subset_data <- all_q_data_imputed[all_q_data_imputed$data == current_date, ]
  
  # Select relevant columns and rename them for consistency
  sub_df <- data.frame(
    location_name = subset_data$location_name,
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    Q = subset_data$Q,
    data = subset_data$data
  )
  
  # Store the sub-dataset in the list
  y.n[[i]] <- sub_df
}

# Check the structure of the first few sub-datasets
str(y.n[1:5])

# RESTRUCTURE all_prec_data_imputed ----
# Create a list to store sub-datasets
x.n <- list()

# Get all unique dates from the precipitation data
all_dates_prec <- unique(all_prec_data_imputed$data)

# Loop over each date to create sub-datasets
for (i in seq_along(all_dates_prec)) {
  current_date <- all_dates_prec[i]
  
  # Extract data for the current date
  subset_data <- all_prec_data_imputed[all_prec_data_imputed$data == current_date, ]
  
  # Select relevant columns and rename them for consistency
  sub_df <- data.frame(
    coordinates.x = subset_data$coordinates.x,
    coordinates.y = subset_data$coordinates.y,
    prec = subset_data$prec,
    data = subset_data$data
  )
  
  # Store the sub-dataset in the list
  x.n[[i]] <- sub_df
}

# Check the structure of the first few sub-datasets
str(x.n[1:5])


# GAM MODEL FOR Y ----
# Define the dependent variable name
y.names <- "Q"
p <- length(y.names)
N <- length(y.n)  # Total number of dates in the y.n list

# Initialize arrays to store predictions, weights, and errors
models.plot <- vector(mode = "list", length = p)
names(models.plot) <- y.names

no.unit.y <- c()
h <- 0

# Dynamically setting the dimension of predictions per unit
max_Tpoints <- max(sapply(y.n, function(subset_data) nrow(subset_data)))

# Adjusting the size of the arrays
check.weight.y <- weight.y <- Y.star <- array(
  NA, 
  dim = c("var" = p, "unit" = N, "Tpoint" = max_Tpoints), 
  dimnames = list(y.names, paste0("N_", seq(1, N)), paste0("T_", seq(1, max_Tpoints)))
)

# Loop over variables (just one variable here: Q)
for (j in seq_len(p)) {
  models.plot[[y.names[j]]] <- list()  # Initialize the model storage
  
  # Loop over each date (i.e., sub-dataset in y.n)
  for (i in seq_len(N)) {
    
    # Extract the current sub-dataset for the date
    subset_data <- y.n[[i]]
    
    # Prepare the data for model fitting
    y <- as.data.frame(subset_data[, y.names[j]])
    x1 <- as.data.frame(subset_data$coordinates.x)
    x2 <- as.data.frame(subset_data$coordinates.y)
    
    # Create a data frame to fit the model
    data <- data.frame("y" = y, "x1" = as.numeric(unlist(x1)), "x2" = as.numeric(unlist(x2)))
    colnames(data) <- c("y", "x1", "x2")
    
    # Count the number of unique points for smoothing
    unique_x1 <- length(unique(data$x1))
    unique_x2 <- length(unique(data$x2))
    
    # Proceed only if we have enough unique points
    if (unique_x1 > 1 && unique_x2 > 1) { 
      
      # Set the smoothing parameters
      k_val <- 25
      data$y <- log1p(data$y)
      
      # Try fitting the model
      prova <- try(
        gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
      )
      
      
      if (!inherits(prova, "try-error")) {
        # Fit the model and store it
        model_i <- gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
        models.plot[[y.names[j]]][[i]] <- list(model = model_i)
        
        # Instead of predicting over the entire grid, predict only at river flow station coordinates
        new_data_pred <- data.frame(
          x1 = data$x1,
          x2 = data$x2
        )
        
        # Make predictions over the available river flow stations
        predictions <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        
        # Calculate weights based on prediction errors
        w_ji <- 1 / ((predictions$se.fit)^2)
        w_ji[is.na(w_ji) | is.infinite(w_ji)] <- 0  # Handling NaNs or Infs
        
        # Adjust the storage of predictions
        num_points <- length(predictions$fit)
        weight.y[y.names[j], i, 1:num_points] <- w_ji
        Y.star[y.names[j], i, 1:num_points] <- predictions$fit
        
        
      } else {  # If model fitting fails
        weight.y[y.names[j], i, ] <- NA
        Y.star[y.names[j], i, ] <- NA
        h <- h + 1
        no.unit.y[h] <- i
      }
      
      # Clean up temporary variables
      rm(unique_x1, unique_x2, k_val)
      
    } else {  # If not enough unique points
      weight.y[y.names[j], i, ] <- NA
      Y.star[y.names[j], i, ] <- NA
      h <- h + 1
      no.unit.y[h] <- i
    }
  }
}

# Print failed dates (if any)
print(no.unit.y) # Must be NULL


# GAM MODEL for X ----
# Define the dependent variable name (precipitation)
x.names <- "prec"
# x.names <- c("prec", "data")
q <- length(x.names)
N <- length(x.n)  # Total number of dates in the x.n list

# Initialize models and arrays to store predictions, weights, and errors
models.plot.x <- vector(mode = "list", length = q)
names(models.plot.x) <- x.names

no.unit.x <- c()
h <- 0

# Correct Initialization
N_T <- nrow(new.data)  

check.weight.x <- weight.x <- X.star <- array(
  NA, 
  dim = c("var" = q, "unit" = N, "Tpoint" = N_T), 
  dimnames = list(x.names, paste0("N_", seq(1, N)), paste0("T_", seq(1, N_T)))
)


# Loop over variables (just one variable here: prec)
for (j in seq_len(q)) {
  models.plot.x[[x.names[j]]] <- list()  # Initialize the model storage
  
  # Loop over each date (i.e., sub-dataset in x.n)
  for (i in seq_len(N)) {
    
    # Extract the current sub-dataset for the date
    subset_data <- x.n[[i]]
    
    # Prepare the data for model fitting
    y <- as.data.frame(subset_data[, x.names[j]])
    x1 <- as.data.frame(subset_data$coordinates.x)
    x2 <- as.data.frame(subset_data$coordinates.y)
    
    # Create a data frame to fit the model
    data <- data.frame("y" = y, "x1" = as.numeric(unlist(x1)), "x2" = as.numeric(unlist(x2)))
    colnames(data) <- c("y", "x1", "x2")
    
    # Count the number of unique points for smoothing
    unique_x1 <- length(unique(data$x1))
    unique_x2 <- length(unique(data$x2))
    
    # Proceed only if we have enough unique points
    if (unique_x1 > 1 && unique_x2 > 1) { 
      
      # Set the smoothing parameters
      k_val <- 25
      data$y <- log1p(data$y)
      
      # Try fitting the model
      prova <- try(
        gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
      )
      
      if (!inherits(prova, "try-error")) {
        # Fit the model and store it
        model_i <- gam(y ~ s(x1, x2, k = k_val), data = data, family = gaussian())
        models.plot[[y.names[j]]][[i]] <- list(model = model_i)
        
        # Prepare the prediction data frame (matching the structure of model data)
        new_data_pred <- data.frame(
          x1 = as.numeric(new.data$coordinates.x),
          x2 = as.numeric(new.data$coordinates.y)
        )
        
        # Make predictions over the grid
        predictions <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        
        # Calculate weights based on prediction errors
        w_ji <- 1 / ((predictions$se.fit)^2)
        w_ji[is.na(w_ji) | is.infinite(w_ji)] <- 0  # Handling NaNs or Infs
        
        # Store the results
        weight.x[x.names[j], i, ] <- w_ji
        X.star[x.names[j], i, ] <- predictions$fit
        
      } else {  # If model fitting fails
        weight.x[x.names[j], i, ] <- NA
        X.star[x.names[j], i, ] <- NA
        h <- h + 1
        no.unit.x[h] <- i
      }
      
      # Clean up temporary variables
      rm(unique_x1, unique_x2, k_val)
      
    } else {  # If not enough unique points
      weight.x[x.names[j], i, ] <- NA
      X.star[x.names[j], i, ] <- NA
      h <- h + 1
      no.unit.x[h] <- i
    }
  }
}

# Print failed dates (if any)
print(no.unit.x) # Must be NULL

sum( is.na( weight.x) ) # 0
sum( is.na( weight.y) ) # 0


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
H.Y <- H.Y / p  


for (j in seq_len(q)) {
  H.X <- H.X + (1 / q) * sigma.X_j[[j]]
}

# Display ranges to check for correctness
range(H.Y) # 0.006070222 0.215345285
range(H.X) # -0.0005721219  0.0053087805


# ORTHONORMAL BASIS ----
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
# --- Calculate Scores for X (Prec) with uniform‐score fallback ---
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



### EXPLAINED VARIABILITY ----
# --- Explained Variability Calculation for Y and X ---

# Initialize lists to store matrices for each basis function
n.Y.matrix.list <- vector(mode = "list", length = L)
n.X.matrix.list <- vector(mode = "list", length = H)

# Matrices to store standard deviations
sd.Y.matrix <- matrix(NA, p, L)
sd.X.matrix <- matrix(NA, q, H)

# Assign column names for clarity
colnames(sd.Y.matrix) <- paste("l_", seq(1, L, by = 1))
colnames(sd.X.matrix) <- paste("h_", seq(1, H, by = 1))

rownames(sd.Y.matrix) <- y.names
rownames(sd.X.matrix) <- x.names

# --- Calculate Explained Variability for Y (gamma) ---
for (l in seq_len(L)) {
  n.Y.matrix.list[[l]] <- matrix(NA, N, p)
  
  for (n in seq_len(N)) {
    if (!is.na(gamma[1,n,l])) {  # Check if the score is not missing
      n.Y.matrix.list[[l]][n,] <- gamma[,n,l]  # Extract scores for each unit
    }
  }
  
  # Calculate standard deviation for each basis function across all units
  sd.Y.matrix[,l] <- apply(n.Y.matrix.list[[l]], 2, sd, na.rm = TRUE)
}

# --- Calculate Explained Variability for X (chi) ---
for (h in seq_len(H)) {
  n.X.matrix.list[[h]] <- matrix(NA, N, q)
  
  for (n in seq_len(N)) {
    if (!is.na(chi[1,n,h])) {  # Check if the score is not missing
      n.X.matrix.list[[h]][n,] <- chi[,n,h]  # Extract scores for each unit
    }
  }
  
  # Calculate standard deviation for each basis function across all units
  sd.X.matrix[,h] <- apply(n.X.matrix.list[[h]], 2, sd, na.rm = TRUE)
}

# --- Plotting Explained Variability ---
par(mfrow = c(1,1))
par(mar = c(5, 5, 5, 5), cex.main = 2)

# Plot for Y (gamma) with title
matplot(t(sd.Y.matrix), type = "l", lwd = 2.3,
        xlab = "Terms of the expansion Y", 
        ylab = "Standard Deviation", 
        main = "Explained Variability of River Flow Scores",
        cex.axis = 1, cex.lab = 1, cex.main = 1, col = 1:p)

# Plot for X (chi) with title
matplot(t(sd.X.matrix), type = "l", lwd = 2.3,
        xlab = "Terms of the expansion X", 
        ylab = "Standard Deviation", 
        main = "Explained Variability of Precipitation Scores",
        cex.axis = 1, cex.lab = 1, cex.main = 1, col = 1:q)


# COMULATIVE VARIABILITY
### CUMULATIVE EXPLAINED VARIABILITY ----

# --- Calculate Cumulative Explained Variability for Y (gamma) ---
explained_var_Y <- matrix(0, p, L)
colnames(explained_var_Y) <- paste0("L_", seq(1, L))
rownames(explained_var_Y) <- y.names

for (j in seq_len(p)) {
  total_variance_Y <- sum(sd.Y.matrix[j, ]^2)
  cumulative_variance_Y <- cumsum(sd.Y.matrix[j, ]^2)
  explained_var_Y[j, ] <- cumulative_variance_Y / total_variance_Y
}

# --- Plotting Cumulative Explained Variability for Y (Q) ---
plot(1:L, explained_var_Y[1, ], type = "o", lwd = 2, col = "blue",
     xlab = "Number of Basis Functions (L)", ylab = "Cumulative Explained Variability",
     main = "Cumulative Explained Variability - Y (Q)", ylim = c(0, 1),
     cex.axis = 1, cex.lab = 1, cex.main= 1)
grid()


# --- Calculate Cumulative Explained Variability for X (chi) ---
explained_var_X <- matrix(0, q, H)
colnames(explained_var_X) <- paste0("H_", seq(1, H))
rownames(explained_var_X) <- x.names

for (j in seq_len(q)) {
  total_variance_X <- sum(sd.X.matrix[j, ]^2)
  cumulative_variance_X <- cumsum(sd.X.matrix[j, ]^2)
  explained_var_X[j, ] <- cumulative_variance_X / total_variance_X
}

# --- Plotting Cumulative Explained Variability for X (Prec) ---
plot(1:H, explained_var_X[1, ], type = "o", lwd = 2, col = "orange",
     xlab = "Number of Basis Functions (H)", ylab = "Cumulative Explained Variability",
     main = "Cumulative Explained Variability - X (Prec)", ylim = c(0, 1), 
     cex.axis = 1, cex.lab = 1, cex.main= 1)
grid()


# --- Determine the cutoff points for 90% explained variability ---
L_cutoff <- which(explained_var_Y[1, ] >= 0.95)[1]
H_cutoff <- which(explained_var_X[1, ] >= 0.95)[1]

cat("Number of Basis Functions for Y (L) capturing 90% variability:", L_cutoff, "\n")
cat("Number of Basis Functions for X (H) capturing 90% variability:", H_cutoff, "\n")

# REGRESSION ----
# Corrected x Matrix Construction
x <- matrix(NA, N, q * H)
colnames(x) <- paste0("h.", sort(rep(seq_len(H), q)), "_", rep(x.names, H))

for (h in seq_len(H)) {
  x[, h] <- chi[1, , h]
}

# Check Dimensions
print(dim(x))  

# Prepare the response matrix 'y'
y <- matrix(NA, N, p * L)
colnames(y) <- paste0("l.", sort(rep(seq_len(L), p)), "_", rep(y.names, L))

for (l in seq_len(L)) {
  y[, l] <- gamma[1, , l]
}

# Check Dimensions
print(dim(y))  

# Fitting Penalized Regression using cglasso
# Data Preparation for cglasso
data <- datacggm(Y = y, X = x)

# Fit the model using cglasso
model <- cglasso(. ~ ., data = data, nlambda = 2, lambda.min.ratio = 0, nrho = 1)

# Extract Coefficient Matrix B
B_vector <- model$B[, , 2, 1]
B_matrix <- matrix(B_vector, nrow = H+1, ncol = L)

# Display the B matrix
print(B_matrix)

x <- cbind("Intercept" = 1, x)
dim(x)
# Predict gamma using the fitted model
gamma_pred <- x %*% B_matrix

# Calculate Residuals
residuals <- y - gamma_pred

# Mean Squared Error (MSE) -  evaluating accuracy in score space
MSE <- mean(residuals^2, na.rm = TRUE)
cat("Mean Squared Error (MSE):", MSE, "\n") # 0.127146

# Plotting Predicted vs Actual gamma values
plot(y, gamma_pred,
     xlab = "Actual Gamma Values", ylab = "Predicted Gamma Values",
     main = "Predicted vs Actual Gamma Values",
     col = "blue", pch = 19, 
     cex.main = 1)
abline(0, 1, col = "red")

# Visualizing Model ----
n_locations <- length(unique(all_q_data_imputed$location_name)) # Should be 38
locations_names <- unique(all_q_data_imputed$location_name)

# Compute full spatial effects: [N_T x n_locations]
# (Precipitation grid locations) x (River flow locations)
spatial_effects_full <- phi.X %*% B_matrix[-1,] %*% t(phi.Y)

# Prepare DataFrame for plotting
spatial_effects_df <- data.frame(new.data, spatial_effects_full)
colnames(spatial_effects_df)[-(1:2)] <- locations_names

# HEATMAP FOR EVERY LOCATION ----
# Extract unique station coordinates
stations_coords <- all_q_data_imputed %>%
  select(location_name, coordinates.x, coordinates.y) %>%
  distinct()

# Global color scale for consistency across plots
global_effect_range <- range(spatial_effects_df[,-(1:2)], na.rm = TRUE)

# Loop over each river flow location
for (loc in colnames(spatial_effects_df)[-(1:2)]) {
  
  # Add the current location's effect to the dataframe
  spatial_effects_df$current_effect <- spatial_effects_df[[loc]]
  
  # Extract station coordinates for the current location
  station_point <- stations_coords %>% filter(location_name == loc)
  
  # Create heatmap
  heatm <- ggplot() +
    geom_tile(data = spatial_effects_df, 
              aes(x = coordinates.x, y = coordinates.y, fill = current_effect)) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = global_effect_range,
      oob = scales::squish,  # ensures values outside limits are shown as limit colors
      breaks = c(global_effect_range[1], seq(global_effect_range[1], global_effect_range[2], length.out = 5)[-c(1,5)], global_effect_range[2]),
      labels = c(
        paste0(round(global_effect_range[1], 2)),
        round(seq(global_effect_range[1], global_effect_range[2], length.out = 5)[-c(1,5)], 2),
        paste0(round(global_effect_range[2], 2))
      )
    ) + geom_point(data = station_point, 
               aes(x = coordinates.x, y = coordinates.y),
               color = "black", size = 3, inherit.aes = FALSE) +
    labs(
      title = paste("Spatial Effect of Precipitation on River Flow at", loc),
      x = "Longitude",
      y = "Latitude",
      fill = "Effect"
    ) +
    theme_minimal() +
    theme(
      plot.title.position = "plot",           # Makes hjust relative to the plot
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11)
    ) + coord_fixed()
  
  print(heatm)
}




# CROSS VALIDATION ----
# TRAINING PIPELINE FUNCTION
train_pipeline <- function(q_train_data, prec_train_data) {
  ### 1. Restructure Q training data by date (to create list y.n) 
  # Store the sorted unique dates for Q data
  sorted_dates_q <- sort(unique(q_train_data$data))
  y.n <- lapply(sorted_dates_q, function(dt) {
    q_train_data %>% 
      filter(data == dt) %>% 
      select(location_name, Q, data, coordinates.x, coordinates.y)
  })
  N_q <- length(y.n)
  
  ### 2. Restructure Precipitation training data by date (to create list x.n)
  sorted_dates_prec <- sort(unique(prec_train_data$data))
  x.n <- lapply(sorted_dates_prec, function(dt) {
    prec_train_data %>% 
      filter(data == dt) %>% 
      select(coordinates.x, coordinates.y, prec, data)
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
  
  ### 4. GAM MODEL FOR Q (Y) 
  y.names <- "Q"
  p <- 1
  max_Tpoints_q <- max(sapply(y.n, nrow))
  # Initialize arrays for weights and predictions, with dates indexed by the sorted dates
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
      data_temp <- data.frame(y = y_val, x1 = as.numeric(unlist(x1)), x2 = as.numeric(unlist(x2)))
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
      # Here we still compute a simple sum over non-NA values as before.
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
  max_Tpoints_x <- N_T  # use grid size for precipitation predictions
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
      data_temp <- data.frame(y = y_val, x1 = as.numeric(unlist(x1)), x2 = as.numeric(unlist(x2)))
      colnames(data_temp) <- c("y", "x1", "x2")
      
      unique_x1 <- length(unique(data_temp$x1))
      unique_x2 <- length(unique(data_temp$x2))
      if(unique_x1 > 1 && unique_x2 > 1) {
        k_val <- 25
        data_temp$y <- log1p(data_temp$y)
        model_i <- gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian())
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        new_data_pred <- data.frame(x1 = as.numeric(new.data$coordinates.x),
                                    x2 = as.numeric(new.data$coordinates.y))
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
  
  ### 6. TWO BASIS SYSTEMS: Covariance matrices and mean covariance 
  # For Q:
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
        A <- matrix(valid_w * (valid_Y - valid_Y_mean), nrow = 1)  # 1 x (length(valid_indices))
        temp_num <- t(A) %*% A  
        temp_den <- outer(valid_w, valid_w)
        # Embed these into full-size matrices:
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
  
  # For precipitation:
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
  L <- 10  # number of basis functions for Q
  H_basis <- 10  # for precipitation
  phi.Y <- bs(seq(0, 1, length.out = nrow(H.Y_mat)), df = L, intercept = TRUE)
  phi.X <- bs(seq(0, 1, length.out = nrow(H.X_mat)), df = H_basis, intercept = TRUE)
  
  ### 8. SCORES: Compute gamma (for Q) and chi (for prec) ###
  gamma_arr <- array(NA, dim = c(p, N_q, L), dimnames = list(y.names, paste0("N_", 1:N_q), paste0("L_", 1:L)))
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
  
  chi_arr <- array(NA, dim = c(q_val, N_prec, H_basis), dimnames = list(x.names, paste0("N_", 1:N_prec), paste0("H_", 1:H_basis)))
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
  # For X (from chi): use the first observation for each Q date
  x_reg <- matrix(NA, N_q, q_val * H_basis)
  colnames(x_reg) <- paste0("h.", sort(rep(1:H_basis, q_val)), "_", rep(x.names, H_basis))
  for(h in seq_len(H_basis)) {
    x_reg[, h] <- chi_arr[1, , h]
  }
  # For Y (from gamma): use the first observation for each Q date
  y_reg <- matrix(NA, N_q, p * L)
  colnames(y_reg) <- paste0("l.", sort(rep(1:L, p)), "_", rep(y.names, L))
  for(l in seq_len(L)) {
    y_reg[, l] <- gamma_arr[1, , l]
  }
  
  # Fit penalized regression using cglasso
  data_cg <- datacggm(Y = y_reg, X = x_reg)
  model_cg <- cglasso(. ~ ., data = data_cg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)
  B_vector <- model_cg$B[-1, , 2, 1]
  B_matrix <- matrix(B_vector, nrow = H_basis, ncol = L)
  
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
    q_train_dates = sorted_dates_q       # Return the sorted training dates
  ))
}

#############################
### PREDICTION FUNCTION
#############################
predict_pipeline <- function(test_data, model_trained) {
  # For each test observation, we use the corresponding GAM model (from training) for Q
  # to predict Q at the test location. If the test date is not present in the training dates,
  # we use the overall training mean.
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  for(i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data
    # Use match() to ensure the correct index is retrieved from the sorted training dates
    idx <- match(test_date, model_trained$q_train_dates)
    if(is.na(idx)) {
      preds[i] <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(x1 = test_row$coordinates.x, x2 = test_row$coordinates.y)
      pred_Y <- predict(gam_model_Y, newdata = test_coords, type = "response")
      preds[i] <- pred_Y
    }
  }
  
  # Inverse transformation to get back to original Q units
  preds <- exp(preds) - 1
  return(preds)
}

### CROSS-VALIDATION SIMULATION
rmse_results <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_imputed$location_name)
set.seed(1234)  # For reproducibility

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_imputed$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_data <- all_q_data_imputed[-test_idx, ]
  q_test_data  <- all_q_data_imputed[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_data, all_prec_data_imputed)
  
  # Predict on the test set using the prediction function
  predictions <- predict_pipeline(q_test_data, model_trained)
  
  actual <- q_test_data$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results <- rbind(rmse_results, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
}
# 1.039857, 1.037189, 6.112005, 3.446144



# ADDING THE TIME ----
# Construct predictor matrix from chi scores
x <- matrix(NA, N, q * H)
colnames(x) <- paste0("h.", sort(rep(seq_len(H), q)), "_", rep(x.names, H))
for (h in seq_len(H)) {
  x[, h] <- chi[1, , h]   # using the first row if q = 1
}

# Create time basis from the sorted training dates
# Convert dates to numeric, scale to [0,1], and build a B-spline basis
time_numeric <- as.numeric(all_dates)
min_time <- min(time_numeric)
range_time <- max(time_numeric) - min_time
time_scaled <- (time_numeric - min_time) / range_time
L_time <- 10
phi_time <- bs(time_scaled, df = L_time, intercept = TRUE)
colnames(phi_time) <- paste0("timeBasis_", seq_len(ncol(phi_time)))

# Augment the predictor matrix with the time basis functions
x_aug <- cbind(x, phi_time)

# Construct the response matrix from gamma scores
y <- matrix(NA, N, p * L)
colnames(y) <- paste0("l.", sort(rep(seq_len(L), p)), "_", rep(y.names, L))
for (l in seq_len(L)) {
  y[, l] <- gamma[1, , l]
}

# Package data for cglasso (assuming datacggm is available)
data_cg <- datacggm(Y = y, X = x_aug)

# Fit the penalized regression model using cglasso
model <- cglasso(. ~ ., data = data_cg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)

# Extract the coefficient vector for the chosen lambda.
# The first (q*H) coefficients come from the precipitation scores,
# and the remaining L_time coefficients represent the time effect.
B_vector <- model$B[, , 2, 1]
chi_coeff_count <- q * H + 1
B_matrix <- matrix(B_vector[1:chi_coeff_count], nrow = H + 1, ncol = L)
B_time <- B_vector[(chi_coeff_count + 1):(chi_coeff_count + L_time)]

cat("Estimated time coefficients (B_time):\n")
print(B_time)

x_aug <- cbind(1, x_aug)
# Compute predictions from the final regression
gamma_pred <- x_aug %*% matrix(model$B[, , 2, 1], nrow = ncol(x_aug), ncol = ncol(y))
residuals <- y - gamma_pred
MSE <- mean(residuals^2, na.rm = TRUE)
cat("Mean Squared Error (MSE):", MSE, "\n") # 0.0959994 

plot(y, gamma_pred,
     xlab = "Actual Gamma Values", ylab = "Predicted Gamma Values",
     main = "Predicted vs Actual Gamma Values",
     col = "blue", pch = 19, 
     cex.main = 1)
abline(0, 1, col = "red")


# CROSS VALIDATION ----
# TRAINING PIPELINE FUNCTION
train_pipeline <- function(q_train_data, prec_train_data) {
  
  ### 1. Restructure Q training data by date (to create list y.n)
  # Store the sorted unique dates for Q data for consistent ordering.
  sorted_dates_q <- sort(unique(q_train_data$data))
  y.n <- lapply(sorted_dates_q, function(dt) {
    q_train_data %>% 
      filter(data == dt) %>% 
      select(location_name, Q, data, coordinates.x, coordinates.y)
  })
  N_q <- length(y.n)
  
  ### 2. Restructure Precipitation training data by date (to create list x.n)
  sorted_dates_prec <- sort(unique(prec_train_data$data))
  x.n <- lapply(sorted_dates_prec, function(dt) {
    prec_train_data %>% 
      filter(data == dt) %>% 
      select(coordinates.x, coordinates.y, prec, data)
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
  
  ### 4. GAM MODEL FOR Q (Y)
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
        model_i <- model_i <- try(gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian()))
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
  max_Tpoints_x <- N_T  # Use grid size for precipitation predictions.
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
        model_i <- model_i <- try(gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian()))
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        # For precipitation, when predicting on the spatial grid, use a representative time value.
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
  
  ### 6. TWO BASIS SYSTEMS: Covariance matrices and mean covariance 
  # For Q:
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
  
  # For precipitation:
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
  L <- 10  # number of basis functions for Q
  H_basis <- 10  # for precipitation
  phi.Y <- bs(seq(0, 1, length.out = nrow(H.Y_mat)), df = L, intercept = TRUE)
  phi.X <- bs(seq(0, 1, length.out = nrow(H.X_mat)), df = H_basis, intercept = TRUE)
  
  ### 8. SCORES: Compute gamma (for Q) and chi (for prec) ###
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
  # For X (from chi): use the first observation for each Q date
  x_reg <- matrix(NA, N_q, q_val * H_basis)
  colnames(x_reg) <- paste0("h.", sort(rep(1:H_basis, q_val)), "_", rep(x.names, H_basis))
  for(h in seq_len(H_basis)) {
    x_reg[, h] <- chi_arr[1, , h]
  }
  
  # ----- Add Time as a Nonlinear Effect -----
  # Set the number of time basis functions
  L_time <- 10  
  # Convert sorted training dates to numeric and scale them
  time_numeric <- as.numeric(sorted_dates_q)
  min_time <- min(time_numeric)
  range_time <- max(time_numeric) - min_time
  time_scaled <- (time_numeric - min_time) / range_time
  # Create time basis for the training dates
  phi_time <- bs(time_scaled, df = L_time, intercept = TRUE)
  colnames(phi_time) <- paste0("timeBasis_", 1:L_time)
  # Append the time basis columns to x_reg
  x_reg <- cbind(x_reg, phi_time)
  
  # For Y (from gamma): use the first observation for each Q date
  y_reg <- matrix(NA, N_q, p * L)
  colnames(y_reg) <- paste0("l.", sort(rep(1:L, p)), "_", rep(y.names, L))
  for(l in seq_len(L)) {
    y_reg[, l] <- gamma_arr[1, , l]
  }
  
  # Fit penalized regression using cglasso
  data_cg <- datacggm(Y = y_reg, X = x_reg)
  model_cg <- cglasso(. ~ ., data = data_cg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)
  B_vector <- model_cg$B[, , 2, 1]
  # Determine the number of coefficients associated with chi
  chi_coeff_count <- q_val * H_basis +1
  # Extract coefficients for chi part and reshape into matrix
  B_matrix <- matrix(B_vector[1:chi_coeff_count], nrow = H_basis + 1, ncol = L)
  # Extract coefficients for time basis (nonlinear time effect)
  B_time <- B_vector[(chi_coeff_count + 1):(chi_coeff_count + L_time)]
  
  x_reg <- cbind(1, x_reg)
  
  cat("Time Coefficients (B_time):\n")
  print(B_time)
  
  return(list(
    B_matrix = B_matrix,
    B_time = B_time,
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
    L_time = L_time,
    q_train_dates = sorted_dates_q,   # sorted training dates
    time_scaling = list(min_time = min_time, range_time = range_time)
  ))
}

#############################
### PREDICTION FUNCTION
#############################
predict_pipeline <- function(test_data, model_trained) {
  # Retrieve time scaling parameters from model_trained
  min_time <- model_trained$time_scaling$min_time
  range_time <- model_trained$time_scaling$range_time
  L_time <- model_trained$L_time
  
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  for(i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data  # Assumes date can be converted to numeric
    # Use match() to retrieve the index for the GAM models
    idx <- match(test_date, model_trained$q_train_dates)
    if(is.na(idx)) {
      gam_pred <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(x1 = test_row$coordinates.x, 
                                x2 = test_row$coordinates.y)
      gam_pred <- predict(gam_model_Y, newdata = test_coords, type = "response")
    }
    
    # Compute the nonlinear time adjustment for the test observation:
    test_time <- as.numeric(test_date)
    test_time_scaled <- (test_time - min_time) / range_time
    phi_time_test <- bs(test_time_scaled, df = L_time, intercept = TRUE)
    # Convert the 1 x L_time basis to a simple numeric vector:
    phi_time_test_vec <- as.numeric(phi_time_test)
    # Compute time adjustment as the inner product of the time basis and B_time
    time_adjustment <- sum(phi_time_test_vec * model_trained$B_time)
    
    # Final prediction is the sum of the GAM prediction and the time adjustment:
    preds[i] <- gam_pred + time_adjustment
  }
  
  preds <- exp(preds) - 1
  return(preds)
}

### CROSS-VALIDATION SIMULATION
rmse_results <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_imputed$location_name)
set.seed(1234)  

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_imputed$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_data <- all_q_data_imputed[-test_idx, ]
  q_test_data  <- all_q_data_imputed[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_data, all_prec_data_imputed)
  
  # Predict on the test set using the prediction function
  predictions <- predict_pipeline(q_test_data, model_trained)
  
  actual <- q_test_data$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results <- rbind(rmse_results, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
} # 0.4766501, 1.330325, 9.019179, 7.301412, 7.301412, 0.4211794
# THE RMSE CHANGES A LOT IF I CHANGE THE NUMBER OF L_time. (MAYBE BETTER 8)

print(rmse_results)

# CROSS VALIDATION with l_time = 8  ----
# TRAINING PIPELINE FUNCTION
train_pipeline <- function(q_train_data, prec_train_data) {
  
  ### 1. Restructure Q training data by date (to create list y.n)
  # Store the sorted unique dates for Q data for consistent ordering.
  sorted_dates_q <- sort(unique(q_train_data$data))
  y.n <- lapply(sorted_dates_q, function(dt) {
    q_train_data %>% 
      filter(data == dt) %>% 
      select(location_name, Q, data, coordinates.x, coordinates.y)
  })
  N_q <- length(y.n)
  
  ### 2. Restructure Precipitation training data by date (to create list x.n)
  sorted_dates_prec <- sort(unique(prec_train_data$data))
  x.n <- lapply(sorted_dates_prec, function(dt) {
    prec_train_data %>% 
      filter(data == dt) %>% 
      select(coordinates.x, coordinates.y, prec, data)
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
  
  ### 4. GAM MODEL FOR Q (Y)
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
        model_i <- model_i <- try(gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian()))
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
  max_Tpoints_x <- N_T  # Use grid size for precipitation predictions.
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
        model_i <- model_i <- try(gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian()))
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        # For precipitation, when predicting on the spatial grid, use a representative time value.
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
  
  ### 6. TWO BASIS SYSTEMS: Covariance matrices and mean covariance 
  # For Q:
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
  
  # For precipitation:
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
  L <- 10  # number of basis functions for Q
  H_basis <- 10  # for precipitation
  phi.Y <- bs(seq(0, 1, length.out = nrow(H.Y_mat)), df = L, intercept = TRUE)
  phi.X <- bs(seq(0, 1, length.out = nrow(H.X_mat)), df = H_basis, intercept = TRUE)
  
  ### 8. SCORES: Compute gamma (for Q) and chi (for prec) ###
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
  # For X (from chi): use the first observation for each Q date
  x_reg <- matrix(NA, N_q, q_val * H_basis)
  colnames(x_reg) <- paste0("h.", sort(rep(1:H_basis, q_val)), "_", rep(x.names, H_basis))
  for(h in seq_len(H_basis)) {
    x_reg[, h] <- chi_arr[1, , h]
  }
  
  # ----- Add Time as a Nonlinear Effect -----
  # Set the number of time basis functions
  L_time <- 8 
  # Convert sorted training dates to numeric and scale them
  time_numeric <- as.numeric(sorted_dates_q)
  min_time <- min(time_numeric)
  range_time <- max(time_numeric) - min_time
  time_scaled <- (time_numeric - min_time) / range_time
  # Create time basis for the training dates
  phi_time <- bs(time_scaled, df = L_time, intercept = TRUE)
  colnames(phi_time) <- paste0("timeBasis_", 1:L_time)
  # Append the time basis columns to x_reg
  x_reg <- cbind(x_reg, phi_time)
  
  # For Y (from gamma): use the first observation for each Q date
  y_reg <- matrix(NA, N_q, p * L)
  colnames(y_reg) <- paste0("l.", sort(rep(1:L, p)), "_", rep(y.names, L))
  for(l in seq_len(L)) {
    y_reg[, l] <- gamma_arr[1, , l]
  }
  
  # Fit penalized regression using cglasso
  data_cg <- datacggm(Y = y_reg, X = x_reg)
  model_cg <- cglasso(. ~ ., data = data_cg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)
  B_vector <- model_cg$B[, , 2, 1]
  # Determine the number of coefficients associated with chi
  chi_coeff_count <- q_val * H_basis +1
  # Extract coefficients for chi part and reshape into matrix
  B_matrix <- matrix(B_vector[1:chi_coeff_count], nrow = H_basis + 1, ncol = L)
  # Extract coefficients for time basis (nonlinear time effect)
  B_time <- B_vector[(chi_coeff_count + 1):(chi_coeff_count + L_time)]
  
  x_reg <- cbind(1, x_reg)
  
  cat("Time Coefficients (B_time):\n")
  print(B_time)
  
  return(list(
    B_matrix = B_matrix,
    B_time = B_time,
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
    L_time = L_time,
    q_train_dates = sorted_dates_q,   # sorted training dates
    time_scaling = list(min_time = min_time, range_time = range_time)
  ))
}

#############################
### PREDICTION FUNCTION
#############################
predict_pipeline <- function(test_data, model_trained) {
  # Retrieve time scaling parameters from model_trained
  min_time <- model_trained$time_scaling$min_time
  range_time <- model_trained$time_scaling$range_time
  L_time <- model_trained$L_time
  
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  for(i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data  # Assumes date can be converted to numeric
    # Use match() to retrieve the index for the GAM models
    idx <- match(test_date, model_trained$q_train_dates)
    if(is.na(idx)) {
      gam_pred <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(x1 = test_row$coordinates.x, 
                                x2 = test_row$coordinates.y)
      gam_pred <- predict(gam_model_Y, newdata = test_coords, type = "response")
    }
    
    # Compute the nonlinear time adjustment for the test observation:
    test_time <- as.numeric(test_date)
    test_time_scaled <- (test_time - min_time) / range_time
    phi_time_test <- bs(test_time_scaled, df = L_time, intercept = TRUE)
    # Convert the 1 x L_time basis to a simple numeric vector:
    phi_time_test_vec <- as.numeric(phi_time_test)
    # Compute time adjustment as the inner product of the time basis and B_time
    time_adjustment <- sum(phi_time_test_vec * model_trained$B_time)
    
    # Final prediction is the sum of the GAM prediction and the time adjustment:
    preds[i] <- gam_pred + time_adjustment
  }
  
  preds <- exp(preds) - 1
  return(preds)
}

### CROSS-VALIDATION SIMULATION
rmse_results2 <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_imputed$location_name)
set.seed(1234)  # For reproducibility

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_imputed$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_data <- all_q_data_imputed[-test_idx, ]
  q_test_data  <- all_q_data_imputed[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_data, all_prec_data_imputed)
  
  # Predict on the test set using the prediction function
  predictions <- predict_pipeline(q_test_data, model_trained)
  
  actual <- q_test_data$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results2 <- rbind(rmse_results2, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
} 

print(rmse_results2)

# CROSS VALIDATION with l_time = 6 ----
# TRAINING PIPELINE FUNCTION
train_pipeline <- function(q_train_data, prec_train_data) {
  
  ### 1. Restructure Q training data by date (to create list y.n)
  # Store the sorted unique dates for Q data for consistent ordering.
  sorted_dates_q <- sort(unique(q_train_data$data))
  y.n <- lapply(sorted_dates_q, function(dt) {
    q_train_data %>% 
      filter(data == dt) %>% 
      select(location_name, Q, data, coordinates.x, coordinates.y)
  })
  N_q <- length(y.n)
  
  ### 2. Restructure Precipitation training data by date (to create list x.n)
  sorted_dates_prec <- sort(unique(prec_train_data$data))
  x.n <- lapply(sorted_dates_prec, function(dt) {
    prec_train_data %>% 
      filter(data == dt) %>% 
      select(coordinates.x, coordinates.y, prec, data)
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
  
  ### 4. GAM MODEL FOR Q (Y)
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
        model_i <- model_i <- try(gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian()))
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
  max_Tpoints_x <- N_T  # Use grid size for precipitation predictions.
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
        model_i <- model_i <- try(gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian()))
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        # For precipitation, when predicting on the spatial grid, use a representative time value.
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
  
  ### 6. TWO BASIS SYSTEMS: Covariance matrices and mean covariance 
  # For Q:
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
  
  # For precipitation:
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
  L <- 10  # number of basis functions for Q
  H_basis <- 10  # for precipitation
  phi.Y <- bs(seq(0, 1, length.out = nrow(H.Y_mat)), df = L, intercept = TRUE)
  phi.X <- bs(seq(0, 1, length.out = nrow(H.X_mat)), df = H_basis, intercept = TRUE)
  
  ### 8. SCORES: Compute gamma (for Q) and chi (for prec) ###
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
  # For X (from chi): use the first observation for each Q date
  x_reg <- matrix(NA, N_q, q_val * H_basis)
  colnames(x_reg) <- paste0("h.", sort(rep(1:H_basis, q_val)), "_", rep(x.names, H_basis))
  for(h in seq_len(H_basis)) {
    x_reg[, h] <- chi_arr[1, , h]
  }
  
  # ----- Add Time as a Nonlinear Effect -----
  # Set the number of time basis functions
  L_time <- 6  
  # Convert sorted training dates to numeric and scale them
  time_numeric <- as.numeric(sorted_dates_q)
  min_time <- min(time_numeric)
  range_time <- max(time_numeric) - min_time
  time_scaled <- (time_numeric - min_time) / range_time
  # Create time basis for the training dates
  phi_time <- bs(time_scaled, df = L_time, intercept = TRUE)
  colnames(phi_time) <- paste0("timeBasis_", 1:L_time)
  # Append the time basis columns to x_reg
  x_reg <- cbind(x_reg, phi_time)
  
  # For Y (from gamma): use the first observation for each Q date
  y_reg <- matrix(NA, N_q, p * L)
  colnames(y_reg) <- paste0("l.", sort(rep(1:L, p)), "_", rep(y.names, L))
  for(l in seq_len(L)) {
    y_reg[, l] <- gamma_arr[1, , l]
  }
  
  # Fit penalized regression using cglasso
  data_cg <- datacggm(Y = y_reg, X = x_reg)
  model_cg <- cglasso(. ~ ., data = data_cg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)
  B_vector <- model_cg$B[, , 2, 1]
  # Determine the number of coefficients associated with chi
  chi_coeff_count <- q_val * H_basis +1
  # Extract coefficients for chi part and reshape into matrix
  B_matrix <- matrix(B_vector[1:chi_coeff_count], nrow = H_basis + 1, ncol = L)
  # Extract coefficients for time basis (nonlinear time effect)
  B_time <- B_vector[(chi_coeff_count + 1):(chi_coeff_count + L_time)]
  
  x_reg <- cbind(1, x_reg)
  
  cat("Time Coefficients (B_time):\n")
  print(B_time)
  
  return(list(
    B_matrix = B_matrix,
    B_time = B_time,
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
    L_time = L_time,
    q_train_dates = sorted_dates_q,   # sorted training dates
    time_scaling = list(min_time = min_time, range_time = range_time)
  ))
}

#############################
### PREDICTION FUNCTION
#############################
predict_pipeline <- function(test_data, model_trained) {
  # Retrieve time scaling parameters from model_trained
  min_time <- model_trained$time_scaling$min_time
  range_time <- model_trained$time_scaling$range_time
  L_time <- model_trained$L_time
  
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  for(i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data  # Assumes date can be converted to numeric
    # Use match() to retrieve the index for the GAM models
    idx <- match(test_date, model_trained$q_train_dates)
    if(is.na(idx)) {
      gam_pred <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(x1 = test_row$coordinates.x, 
                                x2 = test_row$coordinates.y)
      gam_pred <- predict(gam_model_Y, newdata = test_coords, type = "response")
    }
    
    # Compute the nonlinear time adjustment for the test observation:
    test_time <- as.numeric(test_date)
    test_time_scaled <- (test_time - min_time) / range_time
    phi_time_test <- bs(test_time_scaled, df = L_time, intercept = TRUE)
    # Convert the 1 x L_time basis to a simple numeric vector:
    phi_time_test_vec <- as.numeric(phi_time_test)
    # Compute time adjustment as the inner product of the time basis and B_time
    time_adjustment <- sum(phi_time_test_vec * model_trained$B_time)
    
    # Final prediction is the sum of the GAM prediction and the time adjustment:
    preds[i] <- gam_pred + time_adjustment
  }
  
  preds <- exp(preds) - 1
  return(preds)
}

### CROSS-VALIDATION SIMULATION
rmse_results3 <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_imputed$location_name)
set.seed(1234)  

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_imputed$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_data <- all_q_data_imputed[-test_idx, ]
  q_test_data  <- all_q_data_imputed[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_data, all_prec_data_imputed)
  
  # Predict on the test set using the prediction function
  predictions <- predict_pipeline(q_test_data, model_trained)
  
  actual <- q_test_data$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results3 <- rbind(rmse_results3, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
}

print(rmse_results3)


# CROSS VALIDATION with l_time = 12 ----
# TRAINING PIPELINE FUNCTION
train_pipeline <- function(q_train_data, prec_train_data) {
  
  ### 1. Restructure Q training data by date (to create list y.n)
  # Store the sorted unique dates for Q data for consistent ordering.
  sorted_dates_q <- sort(unique(q_train_data$data))
  y.n <- lapply(sorted_dates_q, function(dt) {
    q_train_data %>% 
      filter(data == dt) %>% 
      select(location_name, Q, data, coordinates.x, coordinates.y)
  })
  N_q <- length(y.n)
  
  ### 2. Restructure Precipitation training data by date (to create list x.n)
  sorted_dates_prec <- sort(unique(prec_train_data$data))
  x.n <- lapply(sorted_dates_prec, function(dt) {
    prec_train_data %>% 
      filter(data == dt) %>% 
      select(coordinates.x, coordinates.y, prec, data)
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
  
  ### 4. GAM MODEL FOR Q (Y)
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
        model_i <- model_i <- try(gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian()))
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
  max_Tpoints_x <- N_T  # Use grid size for precipitation predictions.
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
        model_i <- model_i <- try(gam(y ~ s(x1, x2, k = k_val), data = data_temp, family = gaussian()))
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        # For precipitation, when predicting on the spatial grid, use a representative time value.
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
  
  ### 6. TWO BASIS SYSTEMS: Covariance matrices and mean covariance 
  # For Q:
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
  
  # For precipitation:
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
  L <- 10  # number of basis functions for Q
  H_basis <- 10  # for precipitation
  phi.Y <- bs(seq(0, 1, length.out = nrow(H.Y_mat)), df = L, intercept = TRUE)
  phi.X <- bs(seq(0, 1, length.out = nrow(H.X_mat)), df = H_basis, intercept = TRUE)
  
  ### 8. SCORES: Compute gamma (for Q) and chi (for prec) ###
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
  # For X (from chi): use the first observation for each Q date
  x_reg <- matrix(NA, N_q, q_val * H_basis)
  colnames(x_reg) <- paste0("h.", sort(rep(1:H_basis, q_val)), "_", rep(x.names, H_basis))
  for(h in seq_len(H_basis)) {
    x_reg[, h] <- chi_arr[1, , h]
  }
  
  # ----- Add Time as a Nonlinear Effect -----
  # Set the number of time basis functions
  L_time <- 12
  # Convert sorted training dates to numeric and scale them
  time_numeric <- as.numeric(sorted_dates_q)
  min_time <- min(time_numeric)
  range_time <- max(time_numeric) - min_time
  time_scaled <- (time_numeric - min_time) / range_time
  # Create time basis for the training dates
  phi_time <- bs(time_scaled, df = L_time, intercept = TRUE)
  colnames(phi_time) <- paste0("timeBasis_", 1:L_time)
  # Append the time basis columns to x_reg
  x_reg <- cbind(x_reg, phi_time)
  
  # For Y (from gamma): use the first observation for each Q date
  y_reg <- matrix(NA, N_q, p * L)
  colnames(y_reg) <- paste0("l.", sort(rep(1:L, p)), "_", rep(y.names, L))
  for(l in seq_len(L)) {
    y_reg[, l] <- gamma_arr[1, , l]
  }
  
  # Fit penalized regression using cglasso
  data_cg <- datacggm(Y = y_reg, X = x_reg)
  model_cg <- cglasso(. ~ ., data = data_cg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)
  B_vector <- model_cg$B[, , 2, 1]
  # Determine the number of coefficients associated with chi
  chi_coeff_count <- q_val * H_basis +1
  # Extract coefficients for chi part and reshape into matrix
  B_matrix <- matrix(B_vector[1:chi_coeff_count], nrow = H_basis + 1, ncol = L)
  # Extract coefficients for time basis (nonlinear time effect)
  B_time <- B_vector[(chi_coeff_count + 1):(chi_coeff_count + L_time)]
  
  x_reg <- cbind(1, x_reg)
  
  cat("Time Coefficients (B_time):\n")
  print(B_time)
  
  return(list(
    B_matrix = B_matrix,
    B_time = B_time,
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
    L_time = L_time,
    q_train_dates = sorted_dates_q,   # sorted training dates
    time_scaling = list(min_time = min_time, range_time = range_time)
  ))
}

#############################
### PREDICTION FUNCTION
#############################
predict_pipeline <- function(test_data, model_trained) {
  # Retrieve time scaling parameters from model_trained
  min_time <- model_trained$time_scaling$min_time
  range_time <- model_trained$time_scaling$range_time
  L_time <- model_trained$L_time
  
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  for(i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data  # Assumes date can be converted to numeric
    # Use match() to retrieve the index for the GAM models
    idx <- match(test_date, model_trained$q_train_dates)
    if(is.na(idx)) {
      gam_pred <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(x1 = test_row$coordinates.x, 
                                x2 = test_row$coordinates.y)
      gam_pred <- predict(gam_model_Y, newdata = test_coords, type = "response")
    }
    
    # Compute the nonlinear time adjustment for the test observation:
    test_time <- as.numeric(test_date)
    test_time_scaled <- (test_time - min_time) / range_time
    phi_time_test <- bs(test_time_scaled, df = L_time, intercept = TRUE)
    # Convert the 1 x L_time basis to a simple numeric vector:
    phi_time_test_vec <- as.numeric(phi_time_test)
    # Compute time adjustment as the inner product of the time basis and B_time
    time_adjustment <- sum(phi_time_test_vec * model_trained$B_time)
    
    # Final prediction is the sum of the GAM prediction and the time adjustment:
    preds[i] <- gam_pred + time_adjustment
  }
  
  preds <- exp(preds) - 1
  return(preds)
}

### CROSS-VALIDATION SIMULATION
rmse_results4 <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_imputed$location_name)
set.seed(1234)  # For reproducibility

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_imputed$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_data <- all_q_data_imputed[-test_idx, ]
  q_test_data  <- all_q_data_imputed[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_data, all_prec_data_imputed)
  
  # Predict on the test set using the prediction function
  predictions <- predict_pipeline(q_test_data, model_trained)
  
  actual <- q_test_data$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results4 <- rbind(rmse_results4, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
}

print(rmse_results4)


# PLOT OF INTERPOLATED PRECIPITATIONS ----
# Step 1: Combine coordinates and identify bounding box

# Extract coordinates from both datasets
coords_prec <- all_prec_data_imputed %>%
  distinct(location_name, coordinates.x, coordinates.y)

# Identify bounding box for the common domain
x_range <- range(coords_prec$coordinates.x)
y_range <- range(coords_prec$coordinates.y)

# Define interpolation grid resolution
x_seq <- seq(x_range[1], x_range[2], length.out = 100)
y_seq <- seq(y_range[1], y_range[2], length.out = 100)

# Create a spatial interpolation grid
interp_grid <- expand.grid(x = x_seq, y = y_seq)

library(gstat)
library(sp)
library(dplyr)

# Define the interpolation grid as spatial points
coordinates(interp_grid) <- ~ x + y
gridded(interp_grid) <- TRUE

# Step 2: Spatial interpolation for Precipitation (prec)
# Prepare precipitation data: aggregate mean precipitation by station
prec_stations <- all_prec_data_imputed %>%
  group_by(location_name, coordinates.x, coordinates.y) %>%
  summarize(mean_prec = mean(prec, na.rm = TRUE))

# Convert to spatial object
coordinates(prec_stations) <- ~ coordinates.x + coordinates.y

# Fit variogram for precipitation
prec_vgm <- variogram(mean_prec ~ 1, data = prec_stations)
prec_fit <- fit.variogram(prec_vgm, model = vgm(model = "Mat", kappa = 3))


# Plot the empirical variogram and overlay the fitted variogram model
plot(prec_vgm, prec_fit, 
     main = "Empirical Variogram and Fitted Model for Precipitation",
     xlab = "Distance", 
     ylab = "Semivariance")


# Perform kriging interpolation
prec_kriged <- krige(mean_prec ~ 1, prec_stations, interp_grid, model = prec_fit)

# Add interpolated precipitation to grid
interp_grid$prec <- prec_kriged$var1.pred


# Step 4: Visual check of interpolated results
library(ggplot2)

# Convert to data frame for visualization
interp_df <- as.data.frame(interp_grid)

# Plot interpolated precipitation
ggplot(interp_df, aes(x, y, fill = prec)) +
  geom_tile() +
  scale_fill_viridis_c() +
  ggtitle("Interpolated Precipitation")

library(scales)  # for squish()

# Calculate 5th and 95th percentiles for precipitation to clip colors
lower_limit <- quantile(interp_df$prec, 0.05, na.rm = TRUE)
upper_limit <- quantile(interp_df$prec, 0.95, na.rm = TRUE)

ggplot(interp_df, aes(x, y, fill = prec)) +
  geom_tile() +
  scale_fill_viridis_c(
    limits = c(lower_limit, upper_limit),
    oob = scales::squish  # squish values outside limits to min/max colors
  ) +
  ggtitle("Interpolated Precipitation (Color Scale Clipped)") +
  theme_minimal()





