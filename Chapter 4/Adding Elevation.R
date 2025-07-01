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
library(geodata)
library(terra)
library(dplyr)
library(ggplot2)
library(mgcv)
library(splines)
library(cglasso)
library(imputeTS)

# Step 1: Download SRTM elevation raster for Switzerland (WGS84 by default)
elev <- elevation_30s(country = "CHE", path = tempdir())  # ~1km resolution

# Step 2: Get Swiss cantonal boundaries
swiss_admin <- gadm(country = "CHE", level = 1, path = tempdir())

# Step 3: Filter to Ticino
ticino_vect <- swiss_admin[swiss_admin$NAME_1 == "Ticino", ]

# Step 4: Reproject both elevation and boundary to LV95 (EPSG:2056)
elev_lv95 <- project(elev, "EPSG:2056")
ticino_lv95 <- project(ticino_vect, "EPSG:2056")

# Step 5: Crop and mask elevation to Ticino boundary
elev_ticino_lv95 <- crop(elev_lv95, ticino_lv95)
elev_ticino_lv95 <- mask(elev_ticino_lv95, ticino_lv95)

# Step 6: Plot
plot(elev_ticino_lv95, main = "Elevation in Ticino")

# Step 7: Convert to data frame with LV95 coordinates
elev_df_lv95 <- as.data.frame(elev_ticino_lv95, xy = TRUE, na.rm = TRUE)
colnames(elev_df_lv95) <- c("coordinates.x", "coordinates.y", "elevation")

head(elev_df_lv95)

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
coords_mat <- as.matrix(all_prec_data_imputed[, c("coordinates.x", "coordinates.y")])
pts_vect <- vect(coords_mat, type = "points", crs = "EPSG:2056")

bbox <- ext(pts_vect)
bbox_expanded <- ext(bbox$xmin - 20000, bbox$xmax + 20000,
                     bbox$ymin - 20000, bbox$ymax + 20000)

# STEP 5: Crop to expanded bounding box (includes parts of Italy)
elev_expanded <- crop(elev_combined_lv95, bbox_expanded)

# STEP 6: Plot (should show full coverage)
plot(elev_expanded, main = "Expanded Elevation")

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



###############################################################################
# 2) Extract Ranges from the Datasets and Create a Grid
###############################################################################
lat_range <- range(all_q_data_with_elev$coordinates.y)
long_range <- range(all_q_data_with_elev$coordinates.x)

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
    elevation     = subset_data$elevation,      
    Q             = subset_data$Q,
    data          = subset_data$data
  )
  
  y.n[[i]] <- sub_df
}

str(y.n[1:5])  

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
    elevation     = subset_data$elevation,    
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
y.names <- "Q"
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
    elev <- as.data.frame(subset_data$elevation) 
    
    data <- data.frame("y" = y,
                       "x1" = as.numeric(unlist(x1)),
                       "x2" = as.numeric(unlist(x2)),
                       "elev" = as.numeric(unlist(elev)))  
    colnames(data) <- c("y", "x1", "x2", "elev")
    
    unique_x1 <- length(unique(data$x1))
    unique_x2 <- length(unique(data$x2))
    unique_elev <- length(unique(data$elev))  
    
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
          elev = data$elev  
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
H.Y <- H.Y / p 


for (j in seq_len(q)) {
  H.X <- H.X + (1 / q) * sigma.X_j[[j]]
}

# Display ranges
range(H.Y) # -0.02588097  0.24368886
range(H.X) # -0.0007313096  0.0074061257

# ORTHONORMAL BASIS ----
# Define the number of basis functions
L <- 10  # For Q 
H <-10  # For Prec 

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
L_cutoff <- which(explained_var_Y[1, ] >= 0.90)[1]
H_cutoff <- which(explained_var_X[1, ] >= 0.90)[1]

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

x <- cbind(1, x)

# Predict gamma using the fitted model
gamma_pred <- x %*% B_matrix

# Calculate Residuals
residuals <- y - gamma_pred

# Mean Squared Error (MSE)
MSE <- mean(residuals^2, na.rm = TRUE)
cat("Mean Squared Error (MSE):", MSE, "\n") # 0.1259373 

# Plotting Predicted vs Actual gamma values
plot(y, gamma_pred,
     xlab = "Actual Gamma Values", ylab = "Predicted Gamma Values",
     main = "Predicted vs Actual Gamma Values",
     col = "blue", pch = 19, 
     cex.main = 1)
abline(0, 1, col = "red")



# CROSS VALIDATION ( no intercept) ----
# TRAINING PIPELINE FUNCTION (modified to include elevation)
train_pipeline <- function(q_train_data, prec_train_data) {
  ### 1. Restructure Q training data by date (to create list y.n) 
  sorted_dates_q <- sort(unique(q_train_data$data))
  y.n <- list()
  for(i in seq_along(sorted_dates_q)) {
    dt <- sorted_dates_q[i]
    sub <- q_train_data %>% 
      filter(data == dt) %>% 
      select(location_name, Q, data, coordinates.x, coordinates.y, elevation)  # Include elevation
    y.n[[i]] <- sub
  }
  N_q <- length(y.n)
  
  ### 2. Restructure Precipitation training data by date (to create list x.n)
  unique_dates_prec <- sort(unique(prec_train_data$data))
  x.n <- list()
  for(i in seq_along(unique_dates_prec)) {
    dt <- unique_dates_prec[i]
    sub <- prec_train_data %>% 
      filter(data == dt) %>% 
      select(coordinates.x, coordinates.y, prec, data, elevation)  # Include elevation
    x.n[[i]] <- sub
  }
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
  # Add elevation to the grid (using the mean elevation from Q training data)
  new.data$elevation <- mean(q_train_data$elevation, na.rm = TRUE)
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
      elev <- as.data.frame(dat$elevation)  
      data_temp <- data.frame(
        y = as.numeric(unlist(y_val)),
        x1 = as.numeric(unlist(x1)),
        x2 = as.numeric(unlist(x2)),
        elev = as.numeric(unlist(elev))
      )
      colnames(data_temp) <- c("y", "x1", "x2", "elev")
      
      unique_x1   <- length(unique(data_temp$x1))
      unique_x2   <- length(unique(data_temp$x2))
      unique_elev <- length(unique(data_temp$elev))
      if(unique_x1 > 1 && unique_x2 > 1 && unique_elev > 1) {
        k_val  <- 25
        data_temp$y <- log1p(data_temp$y)
        model_i <- gam(y ~ s(x1, x2, elev, k = k_val), data = data_temp, family = gaussian())
        models.plot[[y.names[j]]][[i]] <- list(model = model_i)
        # Predict at the same locations (now including elevation)
        new_data_pred <- data.frame(
          x1 = data_temp$x1,
          x2 = data_temp$x2,
          elev = data_temp$elev
        )
        preds <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        w_ji <- 1/((preds$se.fit)^2)
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
      if(num_points > 0){
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
      elev <- as.data.frame(dat$elevation)  # new: elevation
      data_temp <- data.frame(
        y = as.numeric(unlist(y_val)),
        x1 = as.numeric(unlist(x1)),
        x2 = as.numeric(unlist(x2)),
        elev = as.numeric(unlist(elev))
      )
      colnames(data_temp) <- c("y", "x1", "x2", "elev")
      
      unique_x1   <- length(unique(data_temp$x1))
      unique_x2   <- length(unique(data_temp$x2))
      unique_elev <- length(unique(data_temp$elev))
      if(unique_x1 > 1 && unique_x2 > 1 && unique_elev > 1) {
        k_val  <- 25
        data_temp$y <- log1p(data_temp$y)
        model_i <- gam(y ~ s(x1, x2, elev, k = k_val),data = data_temp, family = gaussian())
        models.plot.x[[x.names[j]]][[i]] <- list(model = model_i)
        new_data_pred <- data.frame(
          x1 = as.numeric(new.data$coordinates.x),
          x2 = as.numeric(new.data$coordinates.y),
          elev = as.numeric(new.data$elevation)
        )
        preds <- predict(model_i, newdata = new_data_pred, type = "response", se.fit = TRUE)
        w_ji <- 1/((preds$se.fit)^2)
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
      if(num_points > 0){
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
      num_points <- sum(!is.na(curr_w))
      if(num_points > 0) {
        valid_w <- curr_w[1:num_points]
        valid_Y <- curr_Y[1:num_points]
        valid_Y_mean <- Y.mean[y.names[j], 1:num_points]
        A <- matrix(valid_w * (valid_Y - valid_Y_mean), nrow = 1)
        temp_num <- t(A) %*% A
        temp_den <- outer(valid_w, valid_w)
        embed_num <- matrix(0, max_Tpoints_q, max_Tpoints_q)
        embed_den <- matrix(0, max_Tpoints_q, max_Tpoints_q)
        embed_num[1:num_points, 1:num_points] <- temp_num
        embed_den[1:num_points, 1:num_points] <- temp_den
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
      num_points <- sum(!is.na(curr_w))
      if(num_points > 0) {
        valid_w <- curr_w[1:num_points]
        valid_X <- curr_X[1:num_points]
        valid_X_mean <- X.mean[x.names[j], 1:num_points]
        A <- matrix(valid_w * (valid_X - valid_X_mean), nrow = 1)
        temp_num <- t(A) %*% A
        temp_den <- outer(valid_w, valid_w)
        embed_num <- matrix(0, max_Tpoints_x, max_Tpoints_x)
        embed_den <- matrix(0, max_Tpoints_x, max_Tpoints_x)
        embed_num[1:num_points, 1:num_points] <- temp_num
        embed_den[1:num_points, 1:num_points] <- temp_den
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
    H.X_mat <- H.X_mat + (1/q_val) * sigma.X_j[[x.names[j]]]
  }
  
  ### 7. Orthonormal Bases (using splines) 
  L <- 10  # number of basis functions for Q
  H_basis <- 10 # for precipitation
  phi.Y <- bs(seq(0, 1, length.out = max_Tpoints_q), df = L, intercept = TRUE)
  phi.X <- bs(seq(0, 1, length.out = max_Tpoints_x), df = H_basis, intercept = TRUE)
  
  ### 8. SCORES: Compute gamma (for Q) and chi (for prec) ###
  gamma_arr <- array(NA, dim = c(p, N_q, L), dimnames = list(y.names, paste0("N_", 1:N_q), paste0("L_", 1:L)))
  for(j in seq_len(p)) {
    for(i in seq_len(N_q)) {
      curr_w <- weight.y[y.names[j], i, ]
      curr_Y <- Y.star[y.names[j], i, ]
      num_points <- sum(!is.na(curr_w))
      if(num_points > 0) {
        valid_w <- curr_w[1:num_points]
        valid_Y <- curr_Y[1:num_points]
        valid_Y_mean <- Y.mean[y.names[j], 1:num_points]
        gamma_arr[j, i, ] <- tryCatch({
          solve(t(phi.Y) %*% diag(valid_w) %*% phi.Y) %*% 
            t(phi.Y) %*% diag(valid_w) %*% (valid_Y - valid_Y_mean)
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
      num_points <- sum(!is.na(curr_w))
      if(num_points > 0 && sum(curr_w, na.rm = TRUE) != 0) {
        chi_arr[j, i, ] <- tryCatch({
          solve(t(phi.X) %*% diag(curr_w[1:num_points]) %*% phi.X) %*% 
            t(phi.X) %*% diag(curr_w[1:num_points]) %*% (curr_X[1:num_points] - X.mean[x.names[j], 1:num_points])
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
  B_vector <- model_cg$B[, , 2, 1]
  B_matrix <- matrix(B_vector, nrow = H_basis + 1, ncol = L)

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
    q_train_dates = sorted_dates_q,
    gamma = gamma_arr,                # For later use if needed
    chi = chi_arr                     # For later use if needed
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
  for(i in 1:n_test) {
    test_row <- test_data[i, ]
    test_date <- test_row$data
    idx <- which(model_trained$q_train_dates == test_date)
    if(length(idx) == 0) {
      preds[i] <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      # Supply elevation along with x1 and x2 for prediction
      test_coords <- data.frame(
        x1 = test_row$coordinates.x, 
        x2 = test_row$coordinates.y,
        elev = test_row$elevation
      )
      pred_Y <- predict(gam_model_Y, newdata = test_coords, type = "response")
      preds[i] <- pred_Y
    }
  }
  preds <- exp(preds) - 1
  return(preds)
}

### CROSS-VALIDATION SIMULATION
# Use all_q_data_with_elev and all_prec_data_with_elev in place of q_train_data and prec_train_data.
rmse_results <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_with_elev$location_name)
set.seed(1234)  # For reproducibility

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_with_elev$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_subset <- all_q_data_with_elev[-test_idx, ]
  q_test_subset  <- all_q_data_with_elev[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_subset, all_prec_data_with_elev)
  
  # Predict on the test set using the prediction function
  predictions <- predict_pipeline(q_test_subset, model_trained)
  
  actual <- q_test_subset$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results <- rbind(rmse_results, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
}

cat("\n--- RMSE Results for each location ---\n")
print(rmse_results)

# SECOND APPROACH ----
# GAM MODEL FOR Y  ----
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

# Display ranges
range(H.Y) # 0.006079362 0.215365595
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


# 1. Create a normalized index from 0 to 1
x_seq <- seq(0, 1, length.out = nrow(H.Y))

# 2. Recompute basis matrix
L <- 10
phi.Y <- bs(x_seq, df = L, intercept = TRUE)

# 3. Plot the first two basis functions
matplot(x_seq,
        phi.Y[, 1:2],
        type = "l", lty = 1,                    # lines, solid
        col = c("steelblue","firebrick"),
        xlab = "Normalized position (0–1)",
        ylab = "Basis value",
        main = "First two B‐spline bases (River‐flow)")
legend("topright",
       legend = c("Basis 1","Basis 2"),
       col = c("steelblue","firebrick"),
       lty = 1)

###############################################################################
# GAM MODEL for Elev (Elevation) ----
###############################################################################
# This block fits a GAM for the static elevation map.
# Because elevation is static (does not depend on time), we model it as:
#     elev ~ s(coordinates.x, coordinates.y, k = k_val_elev)
# Set the smoothing parameter for the elevation model.
k_val_elev <- 25

# Fit the GAM using the column 'elev'
model.elev <- gam(elevation ~ s(coordinates.x, coordinates.y, k = k_val_elev),
                  data = elev_df_lv95)

# Define a prediction grid for the elevation model.
# Here we use the unique coordinates provided in elev_df_lv95.
new.data.elev <- data.frame(
  coordinates.x = elev_df_lv95$coordinates.x,
  coordinates.y = elev_df_lv95$coordinates.y
)

# Obtain predictions and their standard errors on the grid.
pred_elev <- predict(model.elev, newdata = new.data.elev, type = "response", se.fit = TRUE)

# Compute weights based on the precision (inverse squared standard error).
w_elev <- 1 / ((pred_elev$se.fit)^2)
w_elev[is.na(w_elev) | is.infinite(w_elev)] <- 0

# Extract predicted elevation values.
E.star <- pred_elev$fit

# Calculate the overall mean elevation (using the predicted values).
E.mean <- rep(mean(E.star, na.rm = TRUE), length(E.star))

# Generate B-spline basis functions over a sequence from 0 to 1.
# The length of the sequence is equal to the number of prediction points.
L_elev <- 10  # Number of basis functions for elevation; adjust as needed.
phi.elev <- bs(seq(0, 1, length.out = length(E.star)), df = L_elev, intercept = TRUE)

# Compute the elevation scores via weighted least squares projection.
alpha_elev <- tryCatch({
  solve(t(phi.elev) %*% diag(w_elev) %*% phi.elev) %*% 
    t(phi.elev) %*% diag(w_elev) %*% (E.star - E.mean)
}, error = function(e) {
  rep(0, L_elev)
})

# Print the elevation scores.
print(alpha_elev)

eta <- alpha_elev  # elevation functional scores (length L_elev)

# Retrieve dimensions from the previously calculated arrays:
H <- dim(chi)[3]         # number of basis functions for precipitation
L <- dim(gamma)[3]       # number of basis functions for river flow
N <- dim(gamma)[2]       # number of units/dates (assumed same across arrays)

# Build the design matrix for the regression.
# The predictor matrix X_reg will have (H + L_elev) columns:
#   - The first H columns come from the precipitation scores (chi[1, , h]).
#   - The next L_elev columns are the elevation scores (eta), repeated across units.
X_reg <- matrix(NA, nrow = N, ncol = H + L_elev)
colnames(X_reg) <- c(paste0("prec_basis_", 1:H),
                     paste0("elev_basis_", 1:L_elev))
# Fill in the precipitation scores.
for (h in seq_len(H)) {
  X_reg[, h] <- chi[1, , h]
}
# Fill in the elevation scores (replicated for each unit).
for (j in seq_len(L_elev)) {
  X_reg[, H + j] <- rep(eta[j], N)
}

# The response matrix Y_reg is taken as the river flow functional scores.
# Its dimensions are N x L.
Y_reg <- matrix(NA, nrow = N, ncol = L)
colnames(Y_reg) <- paste0("Q_basis_", 1:L)
for (l in seq_len(L)) {
  Y_reg[, l] <- gamma[1, , l]
}

# Now, we fit a penalized regression model.
data_reg <- datacggm(Y = Y_reg, X = X_reg)

# Fit the penalized regression.
model_reg <- cglasso(. ~ ., data = data_reg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)

# Extract the coefficient matrix.
B_vector<- model_reg$B[-1, , 2, 1]
B_matrix <- matrix(B_vector, nrow = (H + L_elev), ncol = L)

cat("Regression Coefficient Matrix (including elevation):\n")
print(B_matrix)

# Predict the river flow scores and compute residuals.
gamma_pred <- X_reg %*% B_matrix
residuals <- Y_reg - gamma_pred
MSE <- mean(residuals^2, na.rm = TRUE) 
cat("Mean Squared Error (MSE) for the regression:", MSE, "\n") # 0.1320576 

# Plotting Predicted vs Actual gamma values:
plot(Y_reg, gamma_pred,
     xlab = "Actual Gamma Values", 
     ylab = "Predicted Gamma Values",
     main = "Predicted vs Actual Gamma Values",
     col = "blue", pch = 19, cex.main = 1)
abline(0, 1, col = "red")


# Visualizing the effect of elevation----
# Extract unique station coordinates
stations_coords <- all_q_data_with_elev %>%
  select(location_name, coordinates.x, coordinates.y) %>%
  distinct()

n_locations <- length(unique(all_q_data_with_elev$location_name)) # Should be 38
locations_names <- unique(all_q_data_with_elev$location_name)


# 1) Compute elevation effects on the SAME 25×25 grid:

# 1a) Predict elev on that grid
pred_grid <- predict(model.elev,
                     newdata = new.data[, c("coordinates.x","coordinates.y")],
                     type    = "response",
                     se.fit  = TRUE)

E_star_grid  <- pred_grid$fit
se_grid      <- pred_grid$se.fit
w_elev_grid  <- 1/(se_grid^2); w_elev_grid[!is.finite(w_elev_grid)] <- 0

# 1b) Build elevation basis on that grid
phi.elev_grid <- bs(seq(0,1,length.out=nrow(new.data)), df=L_elev, intercept=TRUE)

# 1c) Extract only the elevation‐rows from B_matrix_reg
B_elev <- B_matrix[(H+1):(H+L_elev), , drop=FALSE] 

# 1d) Project to get a matrix [n_grid × n_stations]
spatial_effects_elev_full <- phi.elev_grid %*% B_elev %*% t(phi.Y)

# 2) Build a clean data.frame with no duplicate names
grid_coords <- new.data[, c("coordinates.x","coordinates.y")]
spatial_effects_elev_df <- cbind(
  grid_coords,
  spatial_effects_elev_full
)
colnames(spatial_effects_elev_df) <- c(
  "coordinates.x", 
  "coordinates.y", 
  locations_names
)

# 3) Color scale common to all stations
global_elev_range <- range(spatial_effects_elev_df[, locations_names], na.rm=TRUE)

# 4) Station points
stations_coords <- all_q_data_with_elev %>%
  select(location_name, coordinates.x, coordinates.y) %>%
  distinct()

# 5) Loop and plot
for(loc in locations_names) {
  spatial_effects_elev_df$current_effect <- spatial_effects_elev_df[[loc]]
  station_pt <- stations_coords %>% filter(location_name==loc)
  
  ggplot() +
    geom_tile(
      data = spatial_effects_elev_df,
      aes(x = coordinates.x, y = coordinates.y, fill = current_effect)
    ) +
    scale_fill_gradient2(
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0,
      limits   = global_elev_range
    ) +
    geom_point(
      data        = station_pt,
      aes(x=coordinates.x, y=coordinates.y),
      inherit.aes = FALSE,
      color       = "black",
      size        = 3
    ) +
    labs(
      title = paste("Spatial Effect of Elevation on River Flow at", loc),
      x     = "Longitude",
      y     = "Latitude",
      fill  = "Elev Effect"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10)) +
    coord_fixed() -> elev_heatmap
  
  print(elev_heatmap)
}


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
  
  # (ii) Elevation block:
  # Fit a GAM on the external elevation data (assumed available as elev_df_lv95 with columns: elev, coordinates.x, coordinates.y)
  elev_df_lv95$elev <- as.numeric(unlist(elev_df_lv95$elev))
  elev_df_lv95$coordinates.x <- as.numeric(elev_df_lv95$coordinates.x)
  elev_df_lv95$coordinates.y <- as.numeric(elev_df_lv95$coordinates.y)
  k_val_elev <- 25
  elev_model <- gam(elev ~ s(coordinates.x, coordinates.y, k = k_val_elev),
                    data = elev_df_lv95)
  new.data.elev <- data.frame(coordinates.x = elev_df_lv95$coordinates.x,
                              coordinates.y = elev_df_lv95$coordinates.y)
  pred_elev <- predict(elev_model, newdata = new.data.elev, type = "response", se.fit = TRUE)
  w_elev <- 1 / ((pred_elev$se.fit)^2)
  w_elev[is.na(w_elev) | is.infinite(w_elev)] <- 0
  E.star <- pred_elev$fit
  E.mean <- rep(mean(E.star, na.rm = TRUE), length(E.star))
  L_elev <- 10  # Number of basis functions for elevation
  phi.elev <- bs(seq(0, 1, length.out = length(E.star)), df = L_elev, intercept = TRUE)
  alpha_elev <- tryCatch({
    solve(t(phi.elev) %*% diag(w_elev) %*% phi.elev) %*% 
      t(phi.elev) %*% diag(w_elev) %*% (E.star - E.mean)
  }, error = function(e) {
    rep(0, L_elev)
  })
  eta <- alpha_elev
  # Store elevation scaling parameters for later predictions:
  min_elev <- min(E.star, na.rm = TRUE)
  range_elev <- max(E.star, na.rm = TRUE) - min_elev
  # Create an elevation design matrix by replicating eta for each Q date.
  x_reg_elev <- matrix(rep(eta, each = N_q), nrow = N_q, ncol = L_elev)
  colnames(x_reg_elev) <- paste0("elev_basis_", 1:L_elev)
  
  # Combine all predictor blocks:
  x_reg <- cbind(x_reg_chi, x_reg_elev)
  
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
  elev_coeff_count <- L_elev
  B_matrix <- matrix(B_vector[1:chi_coeff_count], nrow = H_basis, ncol = L)
  B_elev <- B_vector[(chi_coeff_count + 1):(chi_coeff_count + elev_coeff_count)]
  
  cat("Precipitation Coefficient Matrix (B_matrix):\n")
  print(B_matrix)
  cat("Estimated Elevation Coefficients (B_elev):\n")
  print(B_elev)
  
  #############################
  # Return all relevant objects
  #############################
  return(list(
    B_matrix = B_matrix,
    B_elev = B_elev,
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
    L_elev = L_elev,
    q_train_dates = sorted_dates_q,   # sorted training dates for Q
    elev_model = elev_model,          # the static elevation GAM model
    elev_scaling = list(min_elev = min_elev, range_elev = range_elev)
  ))
}

#############################
### PREDICTION FUNCTION
#############################
predict_pipeline <- function(test_data, model_trained) {
  min_elev <- model_trained$elev_scaling$min_elev
  range_elev <- model_trained$elev_scaling$range_elev
  L_elev <- model_trained$L_elev
  
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  for(i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data  # Assumes date can be converted to numeric
    idx <- match(test_date, model_trained$q_train_dates)
    if(is.na(idx)) {
      gam_pred <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(x1 = test_row$coordinates.x, 
                                x2 = test_row$coordinates.y)
      gam_pred <- predict(gam_model_Y, newdata = test_coords, type = "response")
    }
    
    # --- Elevation Adjustment ---
    # Predict elevation using the static elevation GAM model.
    test_elev_pred <- predict(model_trained$elev_model, 
                              newdata = data.frame(coordinates.x = test_row$coordinates.x, 
                                                   coordinates.y = test_row$coordinates.y), 
                              type = "response")
    # Scale the elevation prediction to [0,1] using training scaling.
    elev_scaled_test <- (test_elev_pred - min_elev) / range_elev
    if (is.na(elev_scaled_test) || is.nan(elev_scaled_test) || is.infinite(elev_scaled_test)) {
      phi_elev_test <- rep(0, L_elev)
    } else {
      # Clamp elevation to [0 + ε, 1 - ε] to avoid boundary collisions
      elev_scaled_test <- min(max(elev_scaled_test, 1e-6), 1 - 1e-6)
      phi_elev_test <- bs(elev_scaled_test, df = L_elev, intercept = TRUE)
    }
    
    elev_adjustment <- sum(as.numeric(phi_elev_test) * model_trained$B_elev)
    
    # Final prediction is the sum of the GAM prediction and elevation adjustment.
    preds[i] <- gam_pred + elev_adjustment
  }
  
  preds <- exp(preds) - 1
  return(preds)
}

### CROSS-VALIDATION SIMULATION
rmse_results <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_with_elev$location_name)
set.seed(1234)  # For reproducibility

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_with_elev$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_data <- all_q_data_with_elev[-test_idx, ]
  q_test_data  <- all_q_data_with_elev[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_data, all_prec_data_imputed)
  
  # Predict on the test set using the prediction function.
  predictions <- predict_pipeline(q_test_data, model_trained)
  
  actual <- q_test_data$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results <- rbind(rmse_results, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
}
cat("\n--- RMSE Results for each location ---\n")
print(rmse_results)
# 1.047747, 1.04021 (k = 15)


# ADDING THE TIME ----
###############################################################################
# GAM MODEL for Elev (Elevation) ----
###############################################################################
# Set the smoothing parameter for the elevation model.
k_val_elev <- 25

# Fit the GAM using the column 'elevation'
model.elev <- gam(elevation ~ s(coordinates.x, coordinates.y, k = k_val_elev),
                  data = elev_df_lv95)

# Define a prediction grid for the elevation model.
new.data.elev <- data.frame(
  coordinates.x = elev_df_lv95$coordinates.x,
  coordinates.y = elev_df_lv95$coordinates.y
)

# Obtain predictions and their standard errors on the grid.
pred_elev <- predict(model.elev, newdata = new.data.elev, type = "response", se.fit = TRUE)

# Compute weights based on the precision (inverse squared standard error).
w_elev <- 1 / ((pred_elev$se.fit)^2)
w_elev[is.na(w_elev) | is.infinite(w_elev)] <- 0

# Extract predicted elevation values.
E.star <- pred_elev$fit

# Calculate the overall mean elevation (using the predicted values).
E.mean <- rep(mean(E.star, na.rm = TRUE), length(E.star))

# Generate B-spline basis functions over a sequence from 0 to 1.
L_elev <- 10  # Number of basis functions for elevation
phi.elev <- bs(seq(0, 1, length.out = length(E.star)), df = L_elev, intercept = TRUE)

# Compute the elevation scores via weighted least squares projection.
alpha_elev <- tryCatch({
  solve(t(phi.elev) %*% diag(w_elev) %*% phi.elev) %*% 
    t(phi.elev) %*% diag(w_elev) %*% (E.star - E.mean)
}, error = function(e) {
  rep(0, L_elev)
})

# Print the elevation scores.
print(alpha_elev)
eta <- alpha_elev  # Elevation functional scores (length L_elev)

###############################################################################
# Construct Predictor Matrix with Precipitation, Elevation, and Time --------
###############################################################################
# Retrieve dimensions from the previously calculated arrays:
H <- dim(chi)[3]         # number of basis functions for precipitation
L <- dim(gamma)[3]       # number of basis functions for river flow
N <- dim(gamma)[2]       # number of observational units/dates

# Build the base design matrix with precipitation and elevation scores.
# The predictor matrix X_reg will have (H + L_elev) columns:
X_reg <- matrix(NA, nrow = N, ncol = H + L_elev)
colnames(X_reg) <- c(paste0("prec_basis_", 1:H),
                     paste0("elev_basis_", 1:L_elev))

# Fill in the precipitation scores (assuming q = 1 so we use the first row).
for (h in seq_len(H)) {
  X_reg[, h] <- chi[1, , h]
}

# Fill in the elevation scores (replicated for all units).
for (j in seq_len(L_elev)) {
  X_reg[, H + j] <- rep(eta[j], N)
}

# -----------------------------
# Incorporate the Time Component
# -----------------------------
# (Assuming 'all_dates' is a vector of training dates of class Date)
time_numeric <- as.numeric(all_dates)
min_time <- min(time_numeric)
range_time <- max(time_numeric) - min_time
time_scaled <- (time_numeric - min_time) / range_time

# Create a time basis using B-spline functions.
L_time <- 10  # Number of basis functions for time
phi_time <- bs(time_scaled, df = L_time, intercept = TRUE)
colnames(phi_time) <- paste0("timeBasis_", 1:ncol(phi_time))

# Augment the predictor matrix with the time basis functions.
X_reg_aug <- cbind(X_reg, phi_time)

###############################################################################
# Construct Response Matrix -----------------------------------------------
###############################################################################
# Build the response matrix Y_reg from the river flow functional scores.
Y_reg <- matrix(NA, nrow = N, ncol = L)
colnames(Y_reg) <- paste0("Q_basis_", 1:L)
for (l in seq_len(L)) {
  Y_reg[, l] <- gamma[1, , l]
}

###############################################################################
# Penalized Regression Model using cglasso -------------------------------
###############################################################################
data_reg <- datacggm(Y = Y_reg, X = X_reg_aug)
model_reg <- cglasso(. ~ ., data = data_reg, nlambda = 2, lambda.min.ratio = 0, nrho = 1)

# Extract the coefficient vector.
# The first (H + L_elev) coefficients correspond to precipitation and elevation.
# The next L_time coefficients represent the time effect.
B_vector_reg <- model_reg$B[-1, , 2, 1]
base_count <- H + L_elev
time_coeff <- B_vector_reg[(base_count + 1):(base_count + L_time)]
B_matrix_reg <- matrix(B_vector_reg[1:base_count], nrow = (H + L_elev), ncol = L)

cat("Regression Coefficient Matrix (Precipitation and Elevation):\n")
print(B_matrix_reg)
cat("Estimated Time Coefficients:\n")
print(time_coeff)

###############################################################################
# Predictions and Evaluation -----------------------------------------------
###############################################################################
Y_pred <- X_reg_aug %*% matrix(model_reg$B[-1, , 2, 1], nrow = ncol(X_reg_aug), ncol = ncol(Y_reg))
residuals <- Y_reg - Y_pred
MSE <- mean(residuals^2, na.rm = TRUE)
cat("Mean Squared Error (MSE) for the augmented regression:", MSE, "\n")  #0.1008567 

# Plotting Predicted vs Actual Gamma Values.
plot(as.vector(Y_reg), as.vector(Y_pred),
     xlab = "Actual Gamma Values", 
     ylab = "Predicted Gamma Values",
     main = "Predicted vs Actual Gamma Values (with Time)",
     col = "blue", pch = 19)
abline(0, 1, col = "red")



# CROSS VALIDATION ----
# --- Load Required Libraries ---
#############################
# TRAINING PIPELINE FUNCTION
#############################
train_pipeline <- function(q_train_data, prec_train_data) {
  
  ### 1. Restructure Q training data by date (to create list y.n)
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
  
  # (ii) Elevation block:
  # Fit a GAM on the external elevation data (assumed available as elev_df_lv95 with columns: elev, coordinates.x, coordinates.y)
  elev_df_lv95$elev <- as.numeric(unlist(elev_df_lv95$elev))
  elev_df_lv95$coordinates.x <- as.numeric(elev_df_lv95$coordinates.x)
  elev_df_lv95$coordinates.y <- as.numeric(elev_df_lv95$coordinates.y)
  k_val_elev <- 25
  elev_model <- gam(elev ~ s(coordinates.x, coordinates.y, k = k_val_elev),
                    data = elev_df_lv95)
  new.data.elev <- data.frame(coordinates.x = elev_df_lv95$coordinates.x,
                              coordinates.y = elev_df_lv95$coordinates.y)
  pred_elev <- predict(elev_model, newdata = new.data.elev, type = "response", se.fit = TRUE)
  w_elev <- 1 / ((pred_elev$se.fit)^2)
  w_elev[is.na(w_elev) | is.infinite(w_elev)] <- 0
  E.star <- pred_elev$fit
  E.mean <- rep(mean(E.star, na.rm = TRUE), length(E.star))
  L_elev <- 10  # Number of basis functions for elevation
  phi.elev <- bs(seq(0, 1, length.out = length(E.star)), df = L_elev, intercept = TRUE)
  alpha_elev <- tryCatch({
    solve(t(phi.elev) %*% diag(w_elev) %*% phi.elev) %*% 
      t(phi.elev) %*% diag(w_elev) %*% (E.star - E.mean)
  }, error = function(e) {
    rep(0, L_elev)
  })
  eta <- alpha_elev
  # Store elevation scaling parameters for later predictions:
  min_elev <- min(E.star, na.rm = TRUE)
  range_elev <- max(E.star, na.rm = TRUE) - min_elev
  # Create an elevation design matrix by replicating eta for each Q date.
  x_reg_elev <- matrix(rep(eta, each = N_q), nrow = N_q, ncol = L_elev)
  colnames(x_reg_elev) <- paste0("elev_basis_", 1:L_elev)
  
  # (iii) Time block:
  # Create a time basis from the sorted training dates.
  time_numeric <- as.numeric(sorted_dates_q)
  min_time <- min(time_numeric)
  range_time <- max(time_numeric) - min_time
  time_scaled <- (time_numeric - min_time) / range_time
  L_time <- 10  # Number of time basis functions
  phi_time <- bs(time_scaled, df = L_time, intercept = TRUE)
  colnames(phi_time) <- paste0("timeBasis_", 1:ncol(phi_time))
  
  # Combine all predictor blocks:
  x_reg <- cbind(x_reg_chi, x_reg_elev, phi_time)
  
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
  elev_coeff_count <- L_elev
  time_coeff_count <- L_time
  B_matrix <- matrix(B_vector[1:chi_coeff_count], nrow = H_basis, ncol = L)
  B_elev <- B_vector[(chi_coeff_count + 1):(chi_coeff_count + elev_coeff_count)]
  B_time <- B_vector[(chi_coeff_count + elev_coeff_count + 1):(chi_coeff_count + elev_coeff_count + time_coeff_count)]
  
  cat("Precipitation Coefficient Matrix (B_matrix):\n")
  print(B_matrix)
  cat("Estimated Elevation Coefficients (B_elev):\n")
  print(B_elev)
  cat("Estimated Time Coefficients (B_time):\n")
  print(B_time)
  
  #############################
  # Return all relevant objects
  #############################
  return(list(
    B_matrix = B_matrix,
    B_elev = B_elev,
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
    L_elev = L_elev,
    q_train_dates = sorted_dates_q,   # sorted training dates for Q
    time_scaling = list(min_time = min_time, range_time = range_time),
    elev_model = elev_model,          # the static elevation GAM model
    elev_scaling = list(min_elev = min_elev, range_elev = range_elev)
  ))
}

#############################
### PREDICTION FUNCTION
#############################
predict_pipeline <- function(test_data, model_trained) {
  # Retrieve time and elevation scaling parameters.
  min_time <- model_trained$time_scaling$min_time
  range_time <- model_trained$time_scaling$range_time
  L_time <- model_trained$L_time
  min_elev <- model_trained$elev_scaling$min_elev
  range_elev <- model_trained$elev_scaling$range_elev
  L_elev <- model_trained$L_elev
  
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  for(i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data  # Assumes date can be converted to numeric
    idx <- match(test_date, model_trained$q_train_dates)
    if(is.na(idx)) {
      gam_pred <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(x1 = test_row$coordinates.x, 
                                x2 = test_row$coordinates.y)
      gam_pred <- predict(gam_model_Y, newdata = test_coords, type = "response")
    }
    
    # --- Time Adjustment ---
    test_time <- as.numeric(test_date)
    test_time_scaled <- (test_time - min_time) / range_time
    phi_time_test <- bs(test_time_scaled, df = L_time, intercept = TRUE)
    time_adjustment <- sum(as.numeric(phi_time_test) * model_trained$B_time)
    
    # --- Elevation Adjustment ---
    # Predict elevation using the static elevation GAM model.
    test_elev_pred <- predict(model_trained$elev_model, 
                              newdata = data.frame(coordinates.x = test_row$coordinates.x, 
                                                   coordinates.y = test_row$coordinates.y), 
                              type = "response")
    # Scale the elevation prediction to [0,1] using training scaling.
    elev_scaled_test <- (test_elev_pred - min_elev) / range_elev
    phi_elev_test <- bs(elev_scaled_test, df = L_elev, intercept = TRUE)
    elev_adjustment <- sum(as.numeric(phi_elev_test) * model_trained$B_elev)
    
    # Final prediction is the sum of the GAM prediction, time, and elevation adjustments.
    preds[i] <- gam_pred + time_adjustment + elev_adjustment
  }
  preds <- exp(preds) - 1
  return(preds)
}

### CROSS-VALIDATION SIMULATION
rmse_results2 <- data.frame(location = character(), RMSE = numeric(), stringsAsFactors = FALSE)
unique_locations <- unique(all_q_data_with_elev$location_name)
set.seed(1234)  # For reproducibility

for(loc in unique_locations) {
  loc_indices <- which(all_q_data_with_elev$location_name == loc)
  if(length(loc_indices) < 5) next  # Skip if not enough data
  test_idx <- sample(loc_indices, size = floor(0.2 * length(loc_indices)))
  q_train_data <- all_q_data_with_elev[-test_idx, ]
  q_test_data  <- all_q_data_with_elev[test_idx, ]
  
  # Fit the full functional spatial regression model on training Q data and full precipitation data.
  model_trained <- train_pipeline(q_train_data, all_prec_data_imputed)
  
  # Predict on the test set using the prediction function.
  predictions <- predict_pipeline(q_test_data, model_trained)
  
  actual <- q_test_data$Q
  rmse_loc <- sqrt(mean((actual - predictions)^2, na.rm = TRUE))
  rmse_results2 <- rbind(rmse_results2, data.frame(location = loc, RMSE = rmse_loc, stringsAsFactors = FALSE))
  cat("Completed cross-validation for location:", loc, "with RMSE:", rmse_loc, "\n")
}




