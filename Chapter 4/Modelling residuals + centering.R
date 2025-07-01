# Load required libraries ---
library(stats)
library(dplyr)
library(mgcv)
library(splines)
library(imputeTS)
library(cglasso)
library(tidyr)
library(ggplot2)

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
    data = subset_data$data, 
    Q_centered =  subset_data$Q_centered
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


# GAM MODEL FOR Y (Q IS SPATIALLY DISCRETE) ----
# Define the dependent variable name
y.names <- "Q_centered"
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
H.Y <- H.Y / p  # This is now a proper average, based on the actual sizes of the covariance matrices


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
cat("Mean Squared Error (MSE):", MSE, "\n") # 0.1162345 

# MODEL RESIDUALS FROM MODEL 1 (BASIC WITHOUT ALTITUDE) ----
# Try to model residuals with 
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
cat("RMSE from AR(1) innovations:", round(RMSE_AR1, 5), "\n") #0.03595 

# summary of AR(1) coefficients
summary(phi) 



# CV ----
# CROSS VALIDATION ----
# --- Load Required Libraries ---
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
  y.names <- "Q_centered"
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
  
  # Predict gamma using the regression model
  x_pred <- x_reg
  gamma_pred <- x_pred %*% B_matrix
  
  # Compute residuals
  residuals <- y_reg - gamma_pred
  
  # AR(1) residual modeling
  phi_vec <- numeric(ncol(residuals))  # AR(1) coefficients
  tau_vec <- numeric(ncol(residuals))  # Innovation SDs
  eta_matrix <- matrix(NA, nrow = nrow(residuals) - 1, ncol = ncol(residuals))
  
  for (l in seq_len(ncol(residuals))) {
    res_l <- residuals[, l]
    
    fit <- try(arima(res_l, order = c(1, 0, 0), include.mean = FALSE), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      phi_vec[l] <- fit$coef["ar1"]
      tau_vec[l] <- sqrt(fit$sigma2)
      eta_matrix[, l] <- res_l[2:length(res_l)] - phi_vec[l] * res_l[1:(length(res_l) - 1)]
    } else {
      phi_vec[l] <- 0
      tau_vec[l] <- NA
      eta_matrix[, l] <- NA
    }
  }
  
  # Compute RMSE on the innovations
  RMSE_AR1 <- mean(eta_matrix^2, na.rm = TRUE)
  
  # Include AR(1) info in return
  return(list(
    B_matrix = B_matrix,
    phi.Y = phi.Y,
    phi.X = phi.X,
    Y.mean = Y.mean,
    models.plot = models.plot,
    models.plot.x = models.plot.x,
    new.data = new.data,
    max_Tpoints_q = max_Tpoints_q,
    max_Tpoints_x = max_Tpoints_x,
    L = L,
    H_basis = H_basis,
    q_train_dates = sorted_dates_q,
    residuals = residuals,
    gamma_pred = gamma_pred,
    phi_ar1 = phi_vec,
    tau_ar1 = tau_vec,
    eta_matrix = eta_matrix,
    RMSE_AR1 = RMSE_AR1
  ))
  
}

station_means <- all_q_data_imputed %>%
  mutate(logQ = log1p(Q)) %>%
  group_by(location_name) %>%
  summarize(mean_logQ = mean(logQ, na.rm = TRUE)) %>%
  arrange(location_name)  # Ensure same order as columns in Q_pred_log

mean_logQ_vec <- station_means$mean_logQ


predict_pipeline <- function(test_data, model_trained) {
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  
  # Store last known residuals for AR(1) correction
  last_resid <- model_trained$residuals[nrow(model_trained$residuals), , drop = FALSE]
  phi <- model_trained$phi_ar1
  
  for (i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data
    idx <- match(test_date, model_trained$q_train_dates)
    
    if (is.na(idx) || is.null(model_trained$models.plot[["Q"]][[idx]])) {
      preds[i] <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      gam_model_Y <- model_trained$models.plot[["Q"]][[idx]]$model
      test_coords <- data.frame(x1 = test_row$coordinates.x, x2 = test_row$coordinates.y)
      
      # 1. Predict using GAM
      pred_Y_log <- tryCatch({
        predict(gam_model_Y, newdata = test_coords, type = "response")
      }, error = function(e) {
        log1p(mean(model_trained$Y.mean, na.rm = TRUE))
      })
      
      # 2. AR(1) correction in score space
      if (!is.null(phi) && !all(is.na(phi))) {
        phi_1 <- phi[1]  # first basis
        resid_1 <- last_resid[1]  # previous residual for this basis
        gamma_correction <- phi_1 * resid_1
        
        # Reconstruct functional prediction using basis function
        phi_1_vec <- model_trained$phi.Y[, 1]
        q_correction <- gamma_correction * phi_1_vec[1]  # Assume test point ~ first basis element
        
        # Add correction to log-space prediction
        pred_Y_log <- pred_Y_log + q_correction
      }
      
      # 3. Inverse transform
      preds[i] <- exp(pred_Y_log + matrix(mean_logQ_vec, nrow = nrow(Q_pred_log), ncol = ncol(pred_Y_log), byrow = TRUE)) - 1
      
    }
  }
  
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
