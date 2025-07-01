# Load required libraries ---
library(stats)
library(dplyr)
library(mgcv)
library(splines)
library(imputeTS)
library(cglasso)

# MODEL RESIDUALS FROM MODEL 3 ----
# FIT AR(1) TO THE RESIDUALS 
# to model 3, ie model with elevation option 1 
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
RMSE_AR1 <- mean(eta_matrix^2, na.rm = TRUE)
cat("RMSE from AR(1) innovations:", round(RMSE_AR1, 5), "\n") # 0.04529

# Summary of AR(1) coefficients
summary(phi)




# CV ----
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
  
  # --- AR(1) Residual Modeling ---
  # Add intercept to x
  x_with_intercept <- cbind(1, x_reg)
  
  # Predict gamma
  gamma_pred <- x_with_intercept %*% B_matrix
  
  # Compute residuals
  residuals <- y_reg - gamma_pred
  
  # Initialize containers
  phi_ar1 <- numeric(ncol(residuals))
  tau_ar1 <- numeric(ncol(residuals))
  eta_matrix <- matrix(NA, nrow = nrow(residuals) - 1, ncol = ncol(residuals))
  
  # Fit AR(1) to each score time series
  for (l in seq_len(ncol(residuals))) {
    res_l <- residuals[, l]
    fit <- try(arima(res_l, order = c(1, 0, 0), include.mean = FALSE), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      phi_ar1[l] <- fit$coef["ar1"]
      tau_ar1[l] <- sqrt(fit$sigma2)
      eps_t <- res_l[2:length(res_l)]
      eps_tm1 <- res_l[1:(length(res_l) - 1)]
      eta_matrix[, l] <- eps_t - phi_ar1[l] * eps_tm1
    } else {
      phi_ar1[l] <- 0
      tau_ar1[l] <- NA
      eta_matrix[, l] <- NA
    }
  }
  
  # Compute innovation-based RMSE
  RMSE_AR1 <- mean(eta_matrix^2, na.rm = TRUE)
  cat("AR(1) innovation-based RMSE:", round(RMSE_AR1, 5), "\n")
  
  
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
    gamma = gamma_arr,
    chi = chi_arr,
    phi_ar1 = phi_ar1,
    tau_ar1 = tau_ar1,
    eta_matrix = eta_matrix,
    RMSE_AR1 = RMSE_AR1
  ))
  
}

predict_pipeline <- function(test_data, model_trained) {
  n_test <- nrow(test_data)
  preds <- numeric(n_test)
  
  for (i in seq_len(n_test)) {
    test_row <- test_data[i, ]
    test_date <- test_row$data
    
    idx <- which(model_trained$q_train_dates == test_date)
    if (length(idx) == 0) {
      # No matching GAM model: fall back to mean
      preds[i] <- mean(model_trained$Y.mean, na.rm = TRUE)
    } else {
      # --- Step 1: Predict precipitation (X) at test location
      gam_model_x <- model_trained$models.plot.x[["prec"]][[idx]]$model
      test_coords <- data.frame(
        x1 = test_row$coordinates.x,
        x2 = test_row$coordinates.y,
        elev = test_row$elevation
      )
      pred_prec <- predict(gam_model_x, newdata = test_coords, type = "response")
      
      # --- Step 2: Project prec prediction onto X basis (φ_X)
      phi_X <- model_trained$phi.X
      H <- model_trained$H_basis
      X_mean <- model_trained$X.mean["prec", ]
      weights <- rep(1, nrow(phi_X))  # fallback: uniform weights
      
      chi_new <- tryCatch({
        solve(t(phi_X) %*% diag(weights) %*% phi_X) %*%
          t(phi_X) %*% diag(weights) %*% (pred_prec - X_mean)
      }, error = function(e) rep(0, H))
      
      # --- Step 3: Predict gamma using cglasso model
      x_with_intercept <- c(1, chi_new)
      gamma_pred <- x_with_intercept %*% model_trained$B_matrix  # (1 x L)
      
      # --- Step 4: AR(1) correction
      phi_ar1 <- model_trained$phi_ar1
      gamma_last <- model_trained$gamma[1, idx, ]  # Last known gamma before this prediction
      gamma_pred_ar1 <- gamma_pred + phi_ar1 * (gamma_last - gamma_pred)
      
      # --- Step 5: Project γ̂ onto φ_Y basis and reconstruct predicted Q
      phi_Y <- model_trained$phi.Y
      Q_pred_components <- phi_Y %*% matrix(gamma_pred_ar1, ncol = 1)
      Q_pred_val <- mean(Q_pred_components, na.rm = TRUE)
      
      preds[i] <- Q_pred_val
    }
  }
  
  preds <- exp(preds) - 1  # Inverse of log1p transform
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

print(rmse_results)

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
cat("RMSE from AR(1) innovations:", round(RMSE_AR1, 5), "\n") #0.04529

# Summary of AR(1) coefficients
summary(phi) 



# CV ----
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
      
      # 2. AR(1) correction in score space (first component only, generalizable)
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
      preds[i] <- exp(pred_Y_log) - 1
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
