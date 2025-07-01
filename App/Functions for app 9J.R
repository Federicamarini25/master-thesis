rm(list = ls(all = TRUE))     
cat("\014")   

# Load required libraries ---
library(dplyr)
library(tidyr)
library(imputeTS)
library(geodata)
library(terra)
library(mgcv)
library(splines)
library(cglasso)
library(stats)
library(ggplot2)
library(lubridate)
library(ggrepel)  
library(ggpubr)
library(ggforce)
library(gridExtra)
library(grid)
library(httr)
library(jsonlite)
library(caret)  
library(tseries)
library(forecast)
library(zoo)  
library(future.apply)
library(animation)
library(furrr)

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
from_date <- as.Date("2024-01-01")
to_date <- as.Date("2025-01-01")

# Parameters for precipitations 
domain_meteo <- "meteo"
parameter_prec <- "Prec"


# Fetch location data to obtain coordinates
locations_df_q <- fetch_locations_data(domain)
if (is.null(locations_df_q) || !"name" %in% names(locations_df_q) || !"coordinates.x" %in% names(locations_df_q) || !"coordinates.y" %in% names(locations_df_q)) {
  stop("Failed to fetch location data for river flow.")
}

location_names <- as.vector(locations_df_q$name)

# Fetch river data ----
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



### FUNCTIONS FOR FUNCTIONAL REGRESSION ----

# Fetch precipitation data
fetch_all_precipitation_data <- function(domain_meteo = "meteo", 
                                         parameter_prec = "Prec", 
                                         from_date, to_date,
                                         exclude_codes = c("auto_57"),
                                         verbose = TRUE) {
  # Step 1: Fetch station metadata
  locations_df_prec <- fetch_locations_data(domain_meteo)
  
  if (is.null(locations_df_prec) || 
      !"name" %in% names(locations_df_prec) || 
      !"coordinates.x" %in% names(locations_df_prec) || 
      !"coordinates.y" %in% names(locations_df_prec)) {
    stop("❌ Failed to fetch location data for precipitation.")
  }
  
  # Step 2: Filter unwanted stations (e.g. malfunctioning sensors)
  locations_df_prec <- locations_df_prec %>%
    filter(!code %in% exclude_codes)
  
  prec_locations <- as.vector(locations_df_prec$name)
  
  # Step 3: Initialize empty storage
  all_prec_data <- data.frame()
  
  # Step 4: Loop over each location
  for (prec_location in prec_locations) {
    location_row <- locations_df_prec[locations_df_prec$name == prec_location, ]
    if (nrow(location_row) == 0) next
    
    location_code <- location_row$code
    coord_x <- location_row$coordinates.x
    coord_y <- location_row$coordinates.y
    
    if (verbose) message("📡 Fetching data for: ", prec_location)
    
    # Fetch time series
    prec_data_response <- fetch_time_series_data(domain_meteo, location_code, parameter_prec, "d", from_date, to_date)
    if (is.null(prec_data_response)) {
      warning("⚠️ Failed to fetch data for ", prec_location)
      next
    }
    
    # Process
    data <- process_and_append_data(prec_data_response, parameter_prec, prec_location, data.frame())
    if (nrow(data) == 0) next
    
    data <- preprocess_data(data)
    
    # Format and append
    data <- data %>%
      select(prec = !!parameter_prec, data) %>%
      mutate(
        location_name  = prec_location,
        coordinates.x  = coord_x,
        coordinates.y  = coord_y
      )
    
    all_prec_data <- rbind(all_prec_data, data)
  }
  
  return(all_prec_data)
}

all_prec_data <- fetch_all_precipitation_data(domain_meteo = "meteo", 
                                              parameter_prec = "Prec", 
                                              from_date, to_date,
                                              exclude_codes = c("auto_57"),
                                              verbose = TRUE)
  
  
# Clean and impute data ----
clean_and_impute_inputs <- function(raw_prec_data, raw_q_data) {
  # Step 1: Remove outliers
  raw_prec_data$prec[raw_prec_data$prec > 500] <- NA
  raw_q_data$Q[raw_q_data$Q > 1000] <- NA
  
  # Step 2: Impute precipitation data (per location)
  prec_imputed <- raw_prec_data %>%
    group_by(location_name) %>%
    arrange(as.Date(data)) %>%
    mutate(
      prec = zoo::na.approx(prec, x = as.Date(data), na.rm = FALSE),
      prec = zoo::na.locf(prec, na.rm = FALSE),
      prec = zoo::na.locf(prec, fromLast = TRUE, na.rm = FALSE)
    ) %>%
    ungroup()
  
  # Step 3: Impute Q data (per location)
  q_imputed <- raw_q_data %>%
    group_by(location_name) %>%
    arrange(as.Date(data)) %>%
    mutate(
      Q = zoo::na.approx(Q, x = as.Date(data), na.rm = FALSE),
      Q = zoo::na.locf(Q, na.rm = FALSE),
      Q = zoo::na.locf(Q, fromLast = TRUE, na.rm = FALSE)
    ) %>%
    ungroup()
  
  # Step 4: Remove locations with all NA
  prec_imputed <- prec_imputed %>%
    group_by(location_name) %>%
    filter(!all(is.na(prec))) %>%
    ungroup()
  
  q_imputed <- q_imputed %>%
    group_by(location_name) %>%
    filter(!all(is.na(Q))) %>%
    ungroup()
  
  # Step 5: Impute remaining NAs in prec with station-specific mean
  prec_imputed <- prec_imputed %>%
    group_by(location_name) %>%
    mutate(prec = ifelse(is.na(prec), mean(prec, na.rm = TRUE), prec)) %>%
    ungroup()
  
  # Step 6: Center river flow (log1p transform + station-level centering)
  q_imputed <- q_imputed %>%
    mutate(logQ = log1p(Q)) %>%
    group_by(location_name) %>%
    mutate(Q_centered = logQ - mean(logQ, na.rm = TRUE)) %>%
    ungroup()
  
  # Output
  return(list(
    prec_imputed = prec_imputed,
    q_imputed = q_imputed
  ))
}

cleaned <- clean_and_impute_inputs(all_prec_data, all_q_data)
prec_ready <- cleaned$prec_imputed
q_ready <- cleaned$q_imputed


# Add elevation to precipitation data ----
add_elevation_to_prec_data <- function(prec_data, buffer = 20000, verbose = TRUE) {
  # STEP 1: Download SRTM tiles for Switzerland and Italy
  if (verbose) message("📥 Downloading elevation data...")
  elev_che <- elevation_30s(country = "CHE", path = tempdir())
  elev_ita <- elevation_30s(country = "ITA", path = tempdir())
  
  # STEP 2: Reproject both to EPSG:2056 (Swiss LV95)
  if (verbose) message("🗺️ Reprojecting rasters to EPSG:2056...")
  elev_che_lv95 <- project(elev_che, "EPSG:2056")
  elev_ita_lv95 <- project(elev_ita, "EPSG:2056")
  
  # STEP 3: Resample Italy to match Switzerland (bilinear for smoothness)
  elev_ita_lv95_resampled <- resample(elev_ita_lv95, elev_che_lv95, method = "bilinear")
  
  # STEP 4: Mosaic
  elev_combined_lv95 <- mosaic(elev_che_lv95, elev_ita_lv95_resampled)
  
  # STEP 5: Expand bounding box around precipitation points
  coords_mat <- as.matrix(prec_data[, c("coordinates.x", "coordinates.y")])
  pts_vect <- vect(coords_mat, type = "points", crs = "EPSG:2056")
  
  bbox <- ext(pts_vect)
  bbox_expanded <- ext(
    bbox$xmin - buffer, bbox$xmax + buffer,
    bbox$ymin - buffer, bbox$ymax + buffer
  )
  
  # STEP 6: Crop to expanded bounding box
  if (verbose) message("✂️ Cropping elevation raster to bounding box...")
  elev_expanded <- crop(elev_combined_lv95, bbox_expanded)
  
  # STEP 7: Extract elevation values at precipitation locations
  prec_points <- vect(prec_data, geom = c("coordinates.x", "coordinates.y"), crs = "EPSG:2056")
  prec_elev <- extract(elev_expanded, prec_points)
  
  # STEP 8: Bind elevation values into original data
  if (verbose) message("🔗 Merging elevation data into precipitation dataset...")
  prec_ready_with_elev <- prec_data %>%
    mutate(elevation = prec_elev[, 2, drop = TRUE])
  
  return(prec_ready_with_elev)
}

# Add elevation to river flow data ----
add_elevation_to_q_data <- function(q_data, elev_raster = NULL, buffer = 20000, verbose = TRUE) {
  # STEP 0: Load or build elevation raster if not supplied
  if (is.null(elev_raster)) {
    if (verbose) message("📥 Downloading and building elevation raster...")
    
    elev_che <- elevation_30s(country = "CHE", path = tempdir())
    elev_ita <- elevation_30s(country = "ITA", path = tempdir())
    elev_che_lv95 <- project(elev_che, "EPSG:2056")
    elev_ita_lv95 <- project(elev_ita, "EPSG:2056")
    elev_ita_lv95_resampled <- resample(elev_ita_lv95, elev_che_lv95, method = "bilinear")
    elev_raster <- mosaic(elev_che_lv95, elev_ita_lv95_resampled)
  }
  
  # STEP 1: Extract bounding box from river flow station coordinates
  coords_mat <- as.matrix(q_data[, c("coordinates.x", "coordinates.y")])
  q_pts_vect <- vect(coords_mat, type = "points", crs = "EPSG:2056")
  
  bbox <- ext(q_pts_vect)
  bbox_expanded <- ext(
    bbox$xmin - buffer, bbox$xmax + buffer,
    bbox$ymin - buffer, bbox$ymax + buffer
  )
  
  # STEP 2: Crop raster to expanded bounding box
  if (verbose) message("✂️ Cropping elevation raster to bounding box...")
  elev_expanded <- crop(elev_raster, bbox_expanded)
  
  # STEP 3: Convert river flow data to SpatVector and extract elevation
  q_points <- vect(q_data, geom = c("coordinates.x", "coordinates.y"), crs = "EPSG:2056")
  q_elev <- extract(elev_expanded, q_points)
  
  # STEP 4: Merge extracted elevation into original data
  if (verbose) message("🔗 Merging elevation data into river flow dataset...")
  q_data_with_elev <- q_data %>%
    mutate(elevation = q_elev[, 2, drop = TRUE])
  
  return(q_data_with_elev)
}

prec_ready_with_elev <- add_elevation_to_prec_data(prec_ready)

q_ready_with_elev <-add_elevation_to_q_data(q_ready)


# Restructure data ----
restructure_daily_lists <- function(prec_data, q_data) {
  # --- Build x.n ---
  all_dates_prec <- unique(prec_data$data)
  x.n <- vector("list", length(all_dates_prec))
  
  for (i in seq_along(all_dates_prec)) {
    current_date <- all_dates_prec[i]
    subset_data <- prec_data[prec_data$data == current_date, ]
    
    x.n[[i]] <- data.frame(
      coordinates.x = subset_data$coordinates.x,
      coordinates.y = subset_data$coordinates.y,
      elevation     = subset_data$elevation,
      prec          = subset_data$prec,
      data          = subset_data$data
    )
  }
  
  # --- Build y.n ---
  all_dates_q <- unique(q_data$data)
  y.n <- vector("list", length(all_dates_q))
  
  for (i in seq_along(all_dates_q)) {
    current_date <- all_dates_q[i]
    subset_data <- q_data[q_data$data == current_date, ]
    
    # Keep rows even if Q is NA!
    y.n[[i]] <- data.frame(
      location_name = subset_data$location_name,
      coordinates.x = subset_data$coordinates.x,
      coordinates.y = subset_data$coordinates.y,
      elevation     = subset_data$elevation,
      Q             = subset_data$Q,  # can be NA
      data          = subset_data$data,
      Q_centered    = subset_data$Q_centered
    )
  }
  
  return(list(x.n = x.n, y.n = y.n))
}

structured <- restructure_daily_lists(prec_ready_with_elev, q_ready_with_elev)
x.n <- structured$x.n
y.n <- structured$y.n

# GAM MODEL FOR Y ----
fit_gam_weights_y <- function(y.n, y.name = "Q_centered", k_val = 25) {
  p <- length(y.name)
  N <- length(y.n)
  max_Tpoints <- max(sapply(y.n, nrow))
  
  models.plot <- vector("list", length = p)
  names(models.plot) <- y.name
  
  weight.y <- Y.star <- array(NA,
                              dim = c(p, N, max_Tpoints),
                              dimnames = list(y.name, paste0("N_", seq_len(N)), paste0("T_", seq_len(max_Tpoints)))
  )
  no.unit.y <- c()
  h <- 0
  
  for (j in seq_len(p)) {
    models.plot[[y.name[j]]] <- list()
    for (i in seq_len(N)) {
      subset_data <- y.n[[i]]
      y <- as.numeric(subset_data[[y.name[j]]])
      x1 <- as.numeric(subset_data$coordinates.x)
      x2 <- as.numeric(subset_data$coordinates.y)
      elev <- as.numeric(subset_data$elevation)
      
      data <- data.frame(y = y, x1 = x1, x2 = x2, elev = elev)
      
      if (length(unique(x1)) > 1 && length(unique(x2)) > 1 && length(unique(elev)) > 1) {
        data$y <- log1p(data$y)
        try_fit <- try(gam(y ~ s(x1, x2, elev, k = k_val), data = data), silent = TRUE)
        if (!inherits(try_fit, "try-error")) {
          model_i <- try_fit
          models.plot[[y.name[j]]][[i]] <- list(model = model_i)
          
          preds <- predict(model_i, newdata = data[, c("x1", "x2", "elev")], type = "response", se.fit = TRUE)
          w <- 1 / (preds$se.fit^2); w[!is.finite(w)] <- 0
          
          weight.y[y.name[j], i, 1:length(w)] <- w
          Y.star[y.name[j], i, 1:length(preds$fit)] <- preds$fit
        } else {
          no.unit.y[h <- h + 1] <- i
        }
      } else {
        no.unit.y[h <- h + 1] <- i
      }
    }
  }
  
  return(list(weight.y = weight.y, Y.star = Y.star, models.plot = models.plot, failed = no.unit.y))
}

# GAM MODEL FOR X ----
fit_gam_weights_x <- function(x.n, x.name = "prec", new.data, k_val = 25) {
  q <- length(x.name)
  N <- length(x.n)
  N_T <- nrow(new.data)
  
  models.plot.x <- vector("list", length = q)
  names(models.plot.x) <- x.name
  
  weight.x <- X.star <- array(NA,
                              dim = c("var" = q, "unit" = N, "Tpoint" = N_T),
                              dimnames = list(x.name, paste0("N_", seq_len(N)), paste0("T_", seq_len(N_T)))
  )
  no.unit.x <- c()
  h <- 0
  
  for (j in seq_len(q)) {
    models.plot.x[[x.name[j]]] <- list()
    
    for (i in seq_len(N)) {
      subset_data <- x.n[[i]]
      y <- as.numeric(subset_data[[x.name[j]]])
      x1 <- as.numeric(subset_data$coordinates.x)
      x2 <- as.numeric(subset_data$coordinates.y)
      elev <- as.numeric(subset_data$elevation)
      
      data <- data.frame(y = y, x1 = x1, x2 = x2, elev = elev)
      
      if (!("elevation" %in% colnames(new.data))) {
        stop("❌ new.data must contain an 'elevation' column.")
      }
      
      
      if (length(unique(x1)) > 1 && length(unique(x2)) > 1 && length(unique(elev)) > 1) {
        data$y <- log1p(data$y)
        try_fit <- try(gam(y ~ s(x1, x2, elev, k = k_val), data = data), silent = TRUE)
        
        if (!inherits(try_fit, "try-error")) {
          model_i <- try_fit
          models.plot.x[[x.name[j]]][[i]] <- list(model = model_i)
          
          new_data_pred <- data.frame(
            x1   = as.numeric(new.data$coordinates.x),
            x2   = as.numeric(new.data$coordinates.y),
            elev = as.numeric(new.data$elevation)
          )
          
          preds <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
          w <- 1 / (preds$se.fit^2); w[!is.finite(w)] <- 0
          
          weight.x[x.name[j], i, ] <- w
          X.star[x.name[j], i, ] <- preds$fit
        } else {
          no.unit.x[h <- h + 1] <- i
        }
      } else {
        no.unit.x[h <- h + 1] <- i
      }
    }
  }
  
  return(list(weight.x = weight.x, X.star = X.star, models.plot.x = models.plot.x, failed = no.unit.x))
}

res_y <- fit_gam_weights_y(y.n)
weight.y <- res_y$weight.y
Y.star   <- res_y$Y.star


res_x <- fit_gam_weights_x(x.n, new.data = new.data)
weight.x <- res_x$weight.x
X.star   <- res_x$X.star

# Gamma scores ----
compute_gamma_scores <- function(Y.star, weight.y, Y.mean, phi.Y) {
  p <- dim(Y.star)[1]
  N <- dim(Y.star)[2]
  L <- ncol(phi.Y)
  
  gamma <- array(NA, dim = c(p, N, L),
                 dimnames = list(dimnames(Y.star)[[1]], paste0("N_", seq_len(N)), paste0("L_", seq_len(L))))
  
  for (j in seq_len(p)) {
    for (i in seq_len(N)) {
      w <- weight.y[j, i, ]
      y_star <- Y.star[j, i, ]
      y_mean <- Y.mean[j, ]
      
      num_points <- sum(!is.na(w))
      if (num_points > 0) {
        w_valid <- w[1:num_points]
        y_star_valid <- y_star[1:num_points]
        y_mean_valid <- y_mean[1:num_points]
        
        tryCatch({
          gamma[j, i, ] <- solve(
            t(phi.Y) %*% diag(w_valid) %*% phi.Y
          ) %*% t(phi.Y) %*% diag(w_valid) %*% (y_star_valid - y_mean_valid)
        }, error = function(e) {
          gamma[j, i, ] <- rep(0, L)
        })
      } else {
        gamma[j, i, ] <- rep(0, L)
      }
    }
  }
  
  return(gamma)
}

# Chi scores ----
compute_chi_scores <- function(X.star, weight.x, X.mean, phi.X, x.n) {
  q <- dim(X.star)[1]
  N <- dim(X.star)[2]
  H <- ncol(phi.X)
  
  chi <- array(NA, dim = c(q, N, H),
               dimnames = list(dimnames(X.star)[[1]], paste0("N_", seq_len(N)), paste0("H_", seq_len(H))))
  
  for (j in seq_len(q)) {
    for (i in seq_len(N)) {
      subset_data <- x.n[[i]]
      
      # Fallback: no spatial variation
      if (var(subset_data$prec, na.rm = TRUE) == 0) {
        chi[j, i, ] <- c(1, rep(0, H - 1))
        next
      }
      
      w <- weight.x[j, i, ]
      x_star <- X.star[j, i, ]
      x_mean <- X.mean[j, ]
      
      tryCatch({
        chi[j, i, ] <- solve(
          t(phi.X) %*% diag(w) %*% phi.X
        ) %*% t(phi.X) %*% diag(w) %*% (x_star - x_mean)
      }, error = function(e) {
        chi[j, i, ] <- rep(0, H)
      })
    }
  }
  
  return(chi)
}

# Chi (from prec)
chi <- compute_chi_scores(X.star, weight.x, X.mean, phi.X, x.n)

# Gamma (from Q)
gamma <- compute_gamma_scores(Y.star, weight.y, Y.mean, phi.Y)

# Predict river flow ----
predict_river_flow_df <- function(chi, gamma, B_matrix, phi.Y, Y.mean,
                                  y.n.new, all_q_data_imputed, phi,
                                  gamma_sd, return_log = FALSE) {
  # Dimensions
  N <- dim(chi)[2]             # number of prediction days
  L <- ncol(B_matrix)          # number of Q basis functions
  Tpoints <- length(Y.mean)    # number of stations
  
  # Step 1: Predict gamma scores recursively
  gamma_pred <- matrix(NA, nrow = N, ncol = L)
  
  for (t in 2:N) {
    chi_t   <- chi[, t, ]
    chi_tm1 <- chi[, t - 1, ]
    gamma_tm1 <- gamma[1, t - 1, ]  # TRUE gamma scores from training
    
    X_row <- c(1, chi_tm1, chi_t, gamma_tm1)
    gamma_pred[t, ] <- X_row %*% B_matrix
  }
  
  
  # Step 2: AR(1) correction
  residuals_ar1 <- matrix(0, nrow = N, ncol = L)
  for (l in seq_len(L)) {
    if (N >= 3) {
      for (t in 3:N) {
        residuals_ar1[t, l] <- -phi[l] * (gamma_pred[t - 1, l] - gamma_pred[t - 2, l])
      }
    } else if (N == 2) {
      residuals_ar1[2, l] <- -phi[l] * (gamma_pred[2, l] - gamma_pred[1, l])
    } else if (N == 1) {
      gamma_last <- gamma[1, dim(gamma)[2], l]
      gamma_prev <- gamma[1, dim(gamma)[2] - 1, l]
      residuals_ar1[1, l] <- -phi[l] * (gamma_last - gamma_prev)
    }
  }
  
  gamma_corrected <- gamma_pred + residuals_ar1
  
  # Step 3: Reconstruct Q on log1p scale
  Q_pred_log <- gamma_corrected %*% t(phi.Y) +
    matrix(rep(Y.mean, each = N), nrow = N, byrow = FALSE)
  
  # Step 3.1: Compute standard confidence intervals using gamma_sd
  Q_var_log <- rowSums((phi.Y^2) * matrix(gamma_sd^2, nrow = Tpoints, ncol = L, byrow = TRUE))
  Q_se_log <- matrix(rep(sqrt(Q_var_log), each = N), nrow = N)  # [N x Tpoints]
  
  z_score <- qnorm(0.975)
  Q_pred_log_lower <- Q_pred_log - z_score * Q_se_log
  Q_pred_log_upper <- Q_pred_log + z_score * Q_se_log
  
  # Step 4: Add back station-level log(Q) means
  station_means <- all_q_data_imputed %>%
    mutate(logQ = log1p(Q)) %>%
    group_by(location_name) %>%
    summarize(mean_logQ = mean(logQ, na.rm = TRUE)) %>%
    arrange(location_name)
  
  mean_logQ_vec <- station_means$mean_logQ
  
  Q_pred_real <- exp(Q_pred_log + matrix(mean_logQ_vec, nrow = N, ncol = Tpoints, byrow = TRUE)) - 1
  Q_lower     <- exp(Q_pred_log_lower + matrix(mean_logQ_vec, nrow = N, ncol = Tpoints, byrow = TRUE)) - 1
  Q_upper     <- exp(Q_pred_log_upper + matrix(mean_logQ_vec, nrow = N, ncol = Tpoints, byrow = TRUE)) - 1
  
  Q_pred_real[Q_pred_real < 0] <- 0
  Q_lower[Q_lower < 0] <- 0
  Q_upper[Q_upper < 0] <- 0
  
  # Optional: keep log1p version if requested
  Q_output <- if (return_log) Q_pred_log else Q_pred_real
  
  # Step 5: Format results
  results_list <- vector("list", length = N)
  
  for (i in seq_len(N)) {
    sub_data <- y.n.new[[i]]
    n_obs <- nrow(sub_data)
    
    df_day <- data.frame(
      date           = as.Date(sub_data$data),
      location_name  = sub_data$location_name,
      coordinates.x  = sub_data$coordinates.x,
      coordinates.y  = sub_data$coordinates.y,
      Q_true         = sub_data$Q,
      Q_pred         = Q_output[i, 1:n_obs],
      Q_lower        = Q_lower[i, 1:n_obs],
      Q_upper        = Q_upper[i, 1:n_obs]
    )
    
    if (return_log) {
      df_day$Q_pred_log <- Q_pred_log[i, 1:n_obs]
    }
    
    results_list[[i]] <- df_day
  }
  
  results_df <- do.call(rbind, results_list)
  return(results_df)
}

# IMPUTATION FUNCTION ----
# Impute missing data based on selected method
impute_missing_data <- function(data_frame, method) {
  data_frame <- data_frame[order(data_frame$data), ]
  data_frame$time_num <- as.numeric(data_frame$data)
  
  if (method == "linear") {
    data_frame$Q_imputed <- zoo::na.approx(data_frame$Q, x = data_frame$time_num, na.rm = FALSE)
    
  } else if (method == "gam") {
    if (sum(!is.na(data_frame$Q)) > 2) {
      fit <- gam(Q ~ s(time_num), data = data_frame, na.action = na.exclude)
      data_frame$Q_imputed <- predict(fit, newdata = data_frame)
    } else {
      warning("Not enough data points for GAM imputation. No imputation applied.")
      data_frame$Q_imputed <- data_frame$Q
    }
    
  } else if (method == "functional") {
    # --- Load trained model components ---
    load("functional.RData")  # Should contain: phi.X, phi.Y, X.mean, Y.mean, B_matrix, gamma_sd, phi
    
    # Ensure required columns exist
    required_cols <- c("data", "location_name", "coordinates.x", "coordinates.y", "Q")
    # if (!all(required_cols %in% colnames(data_frame))) {
    #   stop("Functional imputation requires data_frame to include: ", paste(required_cols, collapse = ", "))}
    
    # --- Step 1: Fetch precipitation data for same date range ---
    from_date <- min(data_frame$data, na.rm = TRUE)
    to_date   <- max(data_frame$data, na.rm = TRUE)
    
    all_prec_data <- fetch_all_precipitation_data(
      from_date = from_date,
      to_date   = to_date
    )
    
    all_q_data <- data_frame  # use fetched river flow
    
    # --- Step 2: Clean and impute inputs ---
    cleaned <- clean_and_impute_inputs(all_prec_data, all_q_data)
    prec_ready <- cleaned$prec_imputed
    q_ready <- cleaned$q_imputed
    
    # --- Step 3: Add elevation ---
    # Step 2: Add elevation
    # --- Step 3: Add elevation to precipitation ---
    if (exists("prec_ready") && nrow(prec_ready) > 0) {
      prec_ready_with_elev <- tryCatch({
        add_elevation_to_prec_data(prec_ready)
      }, error = function(e) {
        stop("❌ Failed to create prec_ready_with_elev: ", e$message)
      })
    } else {
      stop("❌ prec_ready is empty or missing.")
    }
    
    
    q_ready_with_elev <- add_elevation_to_q_data(q_ready)
    
    # --- Step 4: Restructure by day ---
    structured <- restructure_daily_lists(prec_ready_with_elev, q_ready_with_elev)
    x.n <- structured$x.n
    y.n <- structured$y.n
    
    # --- Step 5: Compute GAM weights ---
    res_y <- fit_gam_weights_y(y.n)
    weight.y <- res_y$weight.y
    Y.star   <- res_y$Y.star
    
    # --- Defensive check for prec_ready_with_elev ---
    if (!exists("prec_ready_with_elev")) {
      stop("❌ prec_ready_with_elev was never created. Check previous elevation merge step.")
    }
    
    res_x <- fit_gam_weights_x(x.n, new.data = new.data)
    weight.x <- res_x$weight.x
    X.star   <- res_x$X.star
    
    # --- Step 6: Compute scores ---
    chi   <- compute_chi_scores(X.star, weight.x, X.mean, phi.X, x.n)
    gamma <- compute_gamma_scores(Y.star, weight.y, Y.mean, phi.Y)
    
    # --- Step 7: Predict Q using functional regression ---
    results_df_ar1 <- predict_river_flow_df(
      chi = chi,
      gamma = gamma,
      B_matrix = B_matrix,
      phi.Y = phi.Y,
      Y.mean = Y.mean,
      y.n.new = y.n,
      all_q_data_imputed = q_ready_with_elev,
      phi = phi,
      gamma_sd = gamma_sd,
      return_log = FALSE
    )
    
    # --- Step 8: Replace only missing Q ---
    # merged <- right_join(data_frame, results_df_ar1, by = c("data" = "date", "location_name"))
    # merged$Q_imputed <- ifelse(is.na(merged$Q), merged$Q_pred, merged$Q)
    # data_frame <- merged
    
    data_frame1 <- data_frame %>%
      mutate(was_missing = is.na(Q))
    
    # Now perform the right join — retain only rows from results_df_ar1
    merged <- right_join(
      data_frame1 %>% mutate(data = as.Date(data)),  # truncate time to match
      results_df_ar1,
      by = c("data" = "date", "location_name", "coordinates.x", "coordinates.y")
    )
    
    data_frame <- merged %>%
      mutate(Q_imputed_f = ifelse(was_missing, Q_pred, Q))
    
  } else {
    stop("Unknown imputation method.")
  }
  
  return(data_frame)
}


imp_func <- impute_missing_data(all_q_data, method ="functional")

filtered_data <- imp_func %>%
  filter(location_name == "Maggia - Lodano")

# 📈 Plot only for that station
ggplot(filtered_data, aes(x = data, y = Q_imputed_f, color = was_missing)) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_color_manual(
    name = "Type",
    values = c("FALSE" = "blue", "TRUE" = "red"),
    labels = c("FALSE" = "Original", "TRUE" = "Imputed")
  ) +
  labs(
    title = "Imputation Plot: Calancasca - Buseno",
    x = "Date", y = "Q Values"
  ) +
  theme_minimal()


# PLOTS with CI ----
# Get unique locations
location_names <- unique(imp_func$location_name)

# Create an empty list to store all plots
plot_list <- list()

# Loop over each location and create ggplot
for (loc in location_names) {
  df_loc <- imp_func %>% filter(location_name == loc)
  
  p <- ggplot(df_loc, aes(x = data)) +
    # Confidence interval ribbon
    geom_ribbon(aes(ymin = Q_lower, ymax = Q_upper), fill = "grey80", alpha = 0.4) +
    
    # True and predicted lines
    geom_line(aes(y = Q_true, color = "True Q")) +
    geom_line(aes(y = Q_pred, color = "Predicted Q")) +
    
    # Labels and theme
    labs(title = loc, x = "Date", y = "River Flow (Q)") +
    scale_color_manual(values = c("True Q" = "black", "Predicted Q" = "red")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  plot_list[[loc]] <- p
}

print(plot_list)
