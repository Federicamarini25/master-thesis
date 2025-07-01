rm(list = ls(all = TRUE))     
cat("\014")                    

# Setting working directory
working_dir = "/Users/federicamarini/Desktop/Master Thesis"            
setwd(working_dir)     

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

# --- Define Parameters ---
domain <- "surfacewater"
parameter <- "Q"
from_date <- as.Date("2020-01-01")
to_date <- as.Date("2023-12-31")


# Fetch location codes
locations_df <- fetch_locations_data(domain)
if (is.null(locations_df) || !"name" %in% names(locations_df)) {
  stop("Failed to fetch location data.")
}

location_names <- as.vector(locations_df$name)

# ---- LINEAR MODEL ----

# Initialize vectors to store MSE values and Q ranges
rmse_values <- c()
q_ranges <- c()

# Loop through each location
for (location_name in location_names) {
  
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if location not found
  
  # Fetch time series data
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  # Skip if data fetch fails
  
  # Process data
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next  # Skip if no data available
  
  # Preprocess data (convert negative values to NA)
  data <- preprocess_data(data)
  
  # Select relevant columns and modify time_num
  data <- data %>%
    select(data, Q) %>%
    mutate(
      time_num = as.numeric(data),  
      X1 = as.numeric(format(data, "%j"))  
    )
  
  # Remove missing values
  data <- na.omit(data)
  
  if (nrow(data) < 10) next  # Skip locations with too few data points
  
  # Calculate Q range
  q_range <- max(data$Q, na.rm = TRUE) - min(data$Q, na.rm = TRUE)
  q_ranges <- c(q_ranges, q_range)
  
  # Define the linear model formula
  formula <- Q ~ time_num  
  
  # Set up cross-validation
  set.seed(1234)
  cv_control <- trainControl(method = "cv", number = 5)  
  
  # Train the model using cross-validation
  model <- train(formula, data = data, method = "lm", trControl = cv_control)
  
  # Make predictions using the trained model
  predictions <- predict(model, data)
  
  # Calculate MSE
  rmse_manual <- sqrt(mean((data$Q - predictions)^2, na.rm = TRUE))
  
  # Store MSE
  rmse_values <- c(rmse_values, rmse_manual)
  
  print(paste("Location:", location_name, "- RMSE:", round(rmse_manual, 3), "- Q Range:", round(q_range, 3)))
}


# ---- GAM MODEL CV ----
rmse_values_gam <- c()
q_ranges <- c()

set.seed(1234)  
for (location_name in location_names) {
  
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if location not found
  
  # Fetch time series data
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  # Skip if data fetch fails
  
  # Process data
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next  # Skip if no data available
  
  # Preprocess data
  data <- preprocess_data(data)
  
  # Select relevant columns and modify time_num
  data <- data %>%
    select(data, Q) %>%
    mutate(
      time_num = as.numeric(data)  
    )
  
  # Remove missing values
  data <- na.omit(data)
  
  if (nrow(data) < 10) next  # Skip locations with too few data points
  
  # Calculate Q range
  q_range <- max(data$Q, na.rm = TRUE) - min(data$Q, na.rm = TRUE)
  q_ranges <- c(q_ranges, q_range)
  
  # --- MANUAL 5-FOLD CROSS-VALIDATION ---
  folds <- createFolds(data$Q, k = 5, list = TRUE, returnTrain = TRUE)
  rmse_folds <- c()
  
  for (fold in folds) {
    train_data <- data[fold, ]
    test_data <- data[-fold, ]
    
    # Train GAM model
    gam_model <- gam(Q ~ s(time_num), data = train_data)
    
    # Predict on test set
    predictions <- predict(gam_model, newdata = test_data)
    
    # Calculate RMSE for this fold
    rmse_fold <- sqrt(mean((test_data$Q - predictions)^2, na.rm = TRUE))
    rmse_folds <- c(rmse_folds, rmse_fold)
  }
  
  # Average RMSE across folds
  rmse_gam <- mean(rmse_folds, na.rm = TRUE)
  rmse_values_gam <- c(rmse_values_gam, rmse_gam)
  
  print(paste("Location:", location_name, "- GAM RMSE:", round(rmse_gam, 3), "- Q Range:", round(q_range, 3)))
}


# DIFFERENT SMOOTH FUNCTIONS ----
# GAM MODEL with Different Basis Functions 
mse_values_gam_default <- c()
mse_values_gam_cyclic <- c()
mse_values_gam_tensor <- c()
mse_values_gam_cubic <- c()

q_ranges <- c()

for (location_name in location_names) {
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if location not found
  
  # Fetch time series data
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  # Skip if data fetch fails
  
  # Process data
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next  # Skip if no data available
  
  # Preprocess data 
  data <- preprocess_data(data)
  
  # Select relevant columns and modify time_num
  data <- data %>%
    select(data, Q) %>%
    mutate(
      time_num = as.numeric(data),  
      X1 = as.numeric(format(data, "%j"))  
    )
  
  # Remove missing values
  data <- na.omit(data)
  
  if (nrow(data) < 10) next  # Skip locations with too few data points
  
  # Calculate Q range
  q_range <- max(data$Q, na.rm = TRUE) - min(data$Q, na.rm = TRUE)
  q_ranges <- c(q_ranges, q_range)
  
  # Train GAM models with different basis functions
  gam_model_default <- gam(Q ~ s(time_num), data = data, na.action = na.exclude)
  gam_model_cyclic <- gam(Q ~ s(time_num, bs = "cc"), data = data, na.action = na.exclude)
  gam_model_tensor <- gam(Q ~ te(time_num), data = data, na.action = na.exclude)
  gam_model_cubic <- gam(Q ~ s(time_num, bs = "cr"), data = data, na.action = na.exclude)
  
  # Predictions & MSE for each GAM model
  predictions_gam_default <- predict(gam_model_default, data)
  predictions_gam_cyclic <- predict(gam_model_cyclic, data)
  predictions_gam_tensor <- predict(gam_model_tensor, data)
  predictions_gam_cubic <- predict(gam_model_cubic, data)
  
  mse_gam_default <- mean((data$Q - predictions_gam_default)^2, na.rm = TRUE)
  mse_gam_cyclic <- mean((data$Q - predictions_gam_cyclic)^2, na.rm = TRUE)
  mse_gam_tensor <- mean((data$Q - predictions_gam_tensor)^2, na.rm = TRUE)
  mse_gam_cubic <- mean((data$Q - predictions_gam_cubic)^2, na.rm = TRUE)
  
  mse_values_gam_default <- c(mse_values_gam_default, mse_gam_default)
  mse_values_gam_cyclic <- c(mse_values_gam_cyclic, mse_gam_cyclic)
  mse_values_gam_tensor <- c(mse_values_gam_tensor, mse_gam_tensor)
  mse_values_gam_cubic <- c(mse_values_gam_cubic, mse_gam_cubic)
  
  print(paste("Location:", location_name, "- Default GAM MSE:", round(mse_gam_default, 3), 
              "- Cyclic GAM MSE:", round(mse_gam_cyclic, 3), "- Tensor GAM MSE:", round(mse_gam_tensor, 3),
              "- Cubic GAM MSE:", round(mse_gam_cubic, 3),
              "- Q Range:", round(q_range, 3)))
}

# Compute the average MSE across all locations for each GAM model
if (length(mse_values_gam_default) > 0) {
  avg_mse_gam_default <- mean(mse_values_gam_default, na.rm = TRUE)
  avg_mse_gam_cyclic <- mean(mse_values_gam_cyclic, na.rm = TRUE)
  avg_mse_gam_tensor <- mean(mse_values_gam_tensor, na.rm = TRUE)
  avg_mse_gam_cubic <- mean(mse_values_gam_cubic, na.rm = TRUE)
  
  print(paste("Average Default GAM MSE across all locations:", round(avg_mse_gam_default, 3)))
  print(paste("Average Cyclic GAM MSE across all locations:", round(avg_mse_gam_cyclic, 3)))
  print(paste("Average Tensor GAM MSE across all locations:", round(avg_mse_gam_tensor, 3)))
  print(paste("Average Cubic GAM MSE across all locations:", round(avg_mse_gam_cubic, 3)))
  
} else {
  print("No valid GAM MSE values were computed.")
}
# no significant changes 


# ---- TIME SERIES ----
# Select the optimal ARIMA parameters for each time series (WITHOUT ASSUMING STATIONARITY)
# Initialize a vector to store AIC values
aic_values <- c()
aic_results <- data.frame(location = character(), AIC = numeric(), stringsAsFactors = FALSE)

# Initialize a data frame to store the best ARIMA models
arima_results <- data.frame(
  Location = character(),
  p = integer(),
  d = integer(),
  q = integer(),
  stringsAsFactors = FALSE
)

# Loop through each location to fit ARIMA model
for (location_name in location_names) {
  
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if location not found
  
  # Fetch time series data
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  # Skip if data fetch fails
  
  # Process data
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next  # Skip if no data available
  
  # Preprocess data (convert negative values to NA)
  data <- preprocess_data(data)
  
  # Ensure we are selecting the relevant column (Q) and a time column (e.g., 'date')
  data <- data %>%
    select(data, Q) %>%  # Make sure 'date' is the correct column for time
    mutate(time_num = as.numeric(as.Date(data)))  # Convert date to numeric
  
  # Remove missing values
  data <- na.omit(data)
  
  if (nrow(data) < 10) next  # Skip locations with too few data points
  
  # Fit ARIMA model using auto.arima() for the time series
  arima_model <- auto.arima(data$Q)
  
  aic_results <- bind_rows(aic_results, data.frame(location = location_name, AIC = arima_model$aic))
  
  # Store the AIC value of the model
  aic_values <- c(aic_values, arima_model$aic)
  
  # Print the results
  cat("Location:", location_name, "\n")
  cat("ARIMA Model Summary:\n")
  print(summary(arima_model))  # Print the ARIMA model summary
  cat("\n")
  
  # Extract ARIMA order (p, d, q)
  arima_order <- arimaorder(arima_model)
  
  # Append results to the data frame
  arima_results <- rbind(arima_results, data.frame(
    Location = location_name,
    p = arima_order[1],
    d = arima_order[2],
    q = arima_order[3]
  ))
}

# Print the final ARIMA model table
print(arima_results)



# ARIMA CV ----
# Cross-validation for ARIMA model
cv_arima_results <- list()
for (location_name in location_names) {
  
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if location not found
  
  # Fetch and process data
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  # Skip if data fetch fails
  
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  data <- preprocess_data(data) %>% select(Q) %>% na.omit()
  
  if (nrow(data) < 10) next  # Skip locations with too few data points
  
  # Create folds for cross-validation
  folds <- createFolds(data$Q, k = 5, list = TRUE)
  
  # Store errors for each fold
  fold_errors <- c()
  
  # Perform cross-validation for ARIMA
  for (fold_index in 1:5) {
    # Split the data into training and testing sets
    train_data <- data[folds[[fold_index]], , drop = FALSE]  # Preserve as data frame
    test_data <- data[-folds[[fold_index]], , drop = FALSE]  # Preserve as data frame
    
    # Fit ARIMA model to the training data
    arima_model <- auto.arima(train_data$Q)
    arima_predictions <- predict(arima_model, n.ahead = nrow(test_data))$pred
    
    # Calculate RMSE or AIC for the fold (RMSE shown here)
    fold_rmse <- sqrt(mean((test_data$Q - arima_predictions)^2))
    fold_errors <- c(fold_errors, fold_rmse)
  }
  
  # Store the average fold error for the location
  cv_arima_results[[location_name]] <- mean(fold_errors)
}

# Output cross-validation results
cv_arima_results_df <- data.frame(
  Location = names(cv_arima_results),
  RMSE = unlist(cv_arima_results),
  stringsAsFactors = FALSE
)

print(cv_arima_results_df)



# ---- Hybrid approach: GAM + TIME SERIES ----
# Initialize an empty list to store cross-validation results
cv_results <- list()

# Number of folds for cross-validation (e.g., K = 5)
k_folds <- 5

# Loop through each location to apply cross-validation
set.seed(1234)
for (location_name in location_names) {
  
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if location not found
  
  # Fetch time series data
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  # Skip if data fetch fails
  
  # Process data
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next  # Skip if no data available
  
  # Preprocess data (convert negative values to NA)
  data <- preprocess_data(data)
  
  # Select relevant columns and modify time_num
  data <- data %>%
    select(data, Q) %>%
    mutate(
      time_num = as.numeric(data),  
      X1 = as.numeric(format(data, "%j"))
    )
  
  # Remove missing values
  data <- na.omit(data)
  
  # Skip if not enough data for cross-validation
  if (nrow(data) < 10) next
  
  # Define the folds for cross-validation
  folds <- createFolds(data$Q, k = k_folds, list = TRUE)
  
  # Store the fold errors for this location
  location_cv_errors <- c()
  
  # Perform cross-validation
  for (fold_index in 1:k_folds) {
    # Split the data into training and testing sets
    train_data <- data[folds[[fold_index]], ]
    test_data <- data[-folds[[fold_index]], ]
    
    # Step 1: Fit the GAM model to capture long-term trends on training data
    gam_model <- gam(Q ~ s(time_num), data = train_data, na.action = na.exclude)
    gam_predictions <- predict(gam_model, newdata = test_data)
    
    # Step 2: Fit the ARIMA model to the non-missing values on training data
    non_na_train_values <- !is.na(train_data$Q)
    arima_model <- auto.arima(train_data$Q[non_na_train_values], stationary = TRUE)
    arima_predictions <- predict(arima_model, n.ahead = nrow(test_data))$pred
    
    # Step 3: Combine GAM and ARIMA predictions using weighted averaging
    combined_predictions <- (gam_predictions + arima_predictions)/2
    
    # Step 4: Calculate the error for this fold (e.g., RMSE or MAE)
    fold_error <- sqrt(mean((test_data$Q - combined_predictions)^2))  # RMSE
    location_cv_errors <- c(location_cv_errors, fold_error)
  }
  
  # Average the errors over all folds for this location
  avg_cv_error <- mean(location_cv_errors)
  
  # Store the result for this location
  cv_results[[location_name]] <- avg_cv_error
  
  # Print the cross-validation results for this location
  cat("Location:", location_name, "\n")
  cat("Average Cross-Validation Error (RMSE):", avg_cv_error, "\n")
  cat("\n")
}

# Create an empty data frame to store results
cv_results_df <- data.frame(
  Location = character(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

# Convert the list to a data frame
for (location_name in names(cv_results)) {
  cv_results_df <- rbind(cv_results_df, data.frame(
    Location = location_name,
    RMSE = cv_results[[location_name]]
  ))
}

# Print the final RMSE table
print(cv_results_df)

# IMPUTATION METHODS COMPARISON ----
set.seed(1234)  # Ensure reproducibility
missing_percentage <- 0.2  # 20% missing data simulation

# Initialize vectors to store RMSE values
rmse_mean <- c()
rmse_linear <- c()
rmse_kalman <- c()
rmse_arima <- c()
rmse_gam <- c()
rmse_hybrid <- c()

for (location_name in location_names) {
  
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  
  
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  
  
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next  
  
  data <- preprocess_data(data)
  
  data <- data %>%
    select(Q) %>%
    mutate(time_num = as.numeric(as.Date(data$data)))
  
  data <- na.omit(data)
  if (nrow(data) < 10) next  
  
  missing_indices <- sample(1:nrow(data), size = floor(missing_percentage * nrow(data)), replace = FALSE)
  data_missing <- data
  data_missing$Q[missing_indices] <- NA  
  
  # Mean Imputation
  data_mean <- data_missing
  data_mean$Q <- na_mean(data_missing$Q)
  
  # Linear Interpolation
  data_interpolation <- data_missing
  data_interpolation$Q <- na_interpolation(data_missing$Q)
  
  # Kalman Filter
  data_kalman <- data_missing
  data_kalman$Q <- na_kalman(data_missing$Q)
  
  # Recursive ARIMA Imputation
  data_arima <- data_missing
  for (i in which(is.na(data_arima$Q))) {
    train_data <- data_arima$Q[1:(i - 1)]
    if (length(train_data) < 5) next  
    model <- auto.arima(train_data)
    data_arima$Q[i] <- forecast(model, h = 1)$mean
  }
  
  # GAM Imputation
  data_gam <- data_missing
  fit_gam <- gam(Q ~ s(time_num), data = data_gam, na.action = na.exclude)
  data_gam$Q <- predict(fit_gam, newdata = data_gam)
  
  # Hybrid Imputation (Average of ARIMA & GAM where both exist)
  data_hybrid <- data_missing
  hybrid_values <- (data_arima$Q + data_gam$Q) / 2
  data_hybrid$Q[is.na(data_hybrid$Q)] <- hybrid_values[is.na(data_hybrid$Q)]
  
  # Compute RMSE for all methods
  original_values <- data$Q[missing_indices]
  rmse_mean <- c(rmse_mean, sqrt(mean((original_values - data_mean$Q[missing_indices])^2, na.rm = TRUE)))
  rmse_linear <- c(rmse_linear, sqrt(mean((original_values - data_interpolation$Q[missing_indices])^2, na.rm = TRUE)))
  rmse_kalman <- c(rmse_kalman, sqrt(mean((original_values - data_kalman$Q[missing_indices])^2, na.rm = TRUE)))
  rmse_arima <- c(rmse_arima, sqrt(mean((original_values - data_arima$Q[missing_indices])^2, na.rm = TRUE)))
  rmse_gam <- c(rmse_gam, sqrt(mean((original_values - data_gam$Q[missing_indices])^2, na.rm = TRUE)))
  rmse_hybrid <- c(rmse_hybrid, sqrt(mean((original_values - data_hybrid$Q[missing_indices])^2, na.rm = TRUE)))

  print(paste("Location:", location_name))
  print(paste("RMSE - Mean:", round(tail(rmse_mean, 1), 3)))
  print(paste("RMSE - Linear:", round(tail(rmse_linear, 1), 3)))
  print(paste("RMSE - Kalman:", round(tail(rmse_kalman, 1), 3)))
  print(paste("RMSE - ARIMA:", round(tail(rmse_arima, 1), 3)))
  print(paste("RMSE - GAM:", round(tail(rmse_gam, 1), 3)))
  print(paste("RMSE - Hybrid:", round(tail(rmse_hybrid, 1), 3)))
}

# Compute overall RMSE summary
total_rmse <- data.frame(
  Method = c("Mean", "Linear", "Kalman", "ARIMA", "GAM", "Hybrid"),
  Avg_RMSE = c(
    mean(rmse_mean, na.rm = TRUE),
    mean(rmse_linear, na.rm = TRUE),
    mean(rmse_kalman, na.rm = TRUE),
    mean(rmse_arima, na.rm = TRUE),
    mean(rmse_gam, na.rm = TRUE),
    mean(rmse_hybrid, na.rm = TRUE),
  )
)

best_method <- total_rmse$Method[which.min(total_rmse$Avg_RMSE)]
print(paste("Best Imputation Method Based on RMSE:", best_method))

