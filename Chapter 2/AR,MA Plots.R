# AR ----
# ---- Libraries ----
library(caret)
library(dplyr)
library(ggplot2)
library(tidyr)

location_names <- c("Bolletta - Porto Ceresio", 
               "Calcaccia - Airolo", "Canale Bonifica - Quartino", 
               "Faloppia - Chiasso", 
               "Laveggio - Mendrisio", "Laveggio - Riva S. Vitale") 

# Set maximum AR lag to test 
max_lag <- 10

# Create an empty results table to store AR model metrics for each location and lag
results <- data.frame(Location = character(),
                      AR_Lag = integer(),
                      AIC = numeric(),
                      BIC = numeric(),
                      RMSE = numeric(),
                      stringsAsFactors = FALSE)

# Loop over each location
for (location_name in location_names) {
  
  # Get the location code from the locations data
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if not found
  
  # Fetch time series data (for parameter Q) for the current location
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  # Skip if data fetch fails
  
  # Process the raw CSV data and preprocess it
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next
  
  data <- preprocess_data(data)
  
  # Keep only the date and Q columns and sort by date (ascending)
  data <- data %>%
    select(data, Q) %>%
    arrange(data)
  
  # Extract the Q values as our time series vector
  Q_ts <- data$Q
  
  # Extract the Q values as our time series vector and remove missing values
  Q_ts <- na.omit(data$Q)
  
  # Skip locations with too few observations (here we require at least max_lag+5 observations)
  if (length(Q_ts) < (max_lag + 5)) next
  
  # Loop over different AR lags (from 1 to max_lag)
  for (lag in 1:max_lag) {
    
    # ---- 5-Fold Cross-Validation for AR(lag) ----
    folds <- createFolds(Q_ts, k = 5, list = TRUE, returnTrain = TRUE)
    rmse_folds <- c()
    
    for (fold in folds) {
      train_idx <- fold
      test_idx <- setdiff(seq_along(Q_ts), train_idx)
      
      train_series <- Q_ts[train_idx]
      test_series <- Q_ts[test_idx]
      
      # Ensure the training series is long enough for an AR model of order 'lag'
      if (length(train_series) <= lag) next
      
      # Fit an AR model
      ar_model <- try(arima(train_series, order = c(lag, 0, 0)), silent = TRUE)
      if (inherits(ar_model, "try-error")) next
      
      # Predict the next n values (n = length(test_series))
      preds <- try(predict(ar_model, n.ahead = length(test_series))$pred, silent = TRUE)
      if (inherits(preds, "try-error")) next
      
      # Compute RMSE for the current fold
      fold_rmse <- sqrt(mean((test_series - preds)^2, na.rm = TRUE))
      rmse_folds <- c(rmse_folds, fold_rmse)
    }
    
    # Average the RMSE across the folds for this AR lag
    avg_rmse <- mean(rmse_folds, na.rm = TRUE)
    
    # ---- Full Model Fit on Entire Series ----
    # Fit the AR model on the full series to extract AIC and BIC values
    full_model <- try(arima(Q_ts, order = c(lag, 0, 0)), silent = TRUE)
    if (inherits(full_model, "try-error")) next
    model_aic <- AIC(full_model)
    model_bic <- BIC(full_model)
    
    # Append these results to our results table
    results <- rbind(results, data.frame(Location = location_name,
                                         AR_Lag = lag,
                                         AIC = model_aic,
                                         BIC = model_bic,
                                         RMSE = avg_rmse,
                                         stringsAsFactors = FALSE))
  }
}

# Print the final results table
print(results)

# ---- Plotting ----
# Each plot will have AR Lag on the x-axis and the corresponding metric on the y-axis.
# Different locations are shown as differently coloured lines.

# Plot AIC vs. AR Lag
plot_aic <- ggplot(results, aes(x = AR_Lag, y = AIC, color = Location)) +
  geom_line() +
  geom_point() +
  labs(title = "AIC vs. AR Lag", x = "AR Lag", y = "AIC") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal()

# Plot BIC vs. AR Lag
plot_bic <- ggplot(results, aes(x = AR_Lag, y = BIC, color = Location)) +
  geom_line() +
  geom_point() +
  labs(title = "BIC vs. AR Lag", x = "AR Lag", y = "BIC") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal()

# Plot RMSE vs. AR Lag (evaluated via 5-Fold CV)
plot_rmse <- ggplot(results, aes(x = AR_Lag, y = RMSE, color = Location)) +
  geom_line() +
  geom_point() +
  labs(title = "RMSE vs. AR Lag (5-Fold CV)", x = "AR Lag", y = "RMSE") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal()

# Display the plots
print(plot_aic)
print(plot_bic)
print(plot_rmse)

# MA ----
max_ma <- 10

# Create an empty results table to store MA model metrics for each location and order
ma_results <- data.frame(Location = character(),
                         MA_Order = integer(),
                         AIC = numeric(),
                         BIC = numeric(),
                         RMSE = numeric(),
                         stringsAsFactors = FALSE)

for (location_name in location_names) {
  
  # Get the location code from the locations data
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if location not found
  
  # Fetch time series data (for parameter Q) for the current location
  ma_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(ma_data_response)) next  # Skip if data fetch fails
  
  # Process and preprocess the data
  data <- process_and_append_data(ma_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next
  
  data <- preprocess_data(data)
  
  # Keep only the date and Q columns and sort by date (ascending)
  data <- data %>%
    select(data, Q) %>%
    arrange(data)
  
  # Extract the Q values as our time series vector and remove any missing values
  Q_ts <- na.omit(data$Q)
  
  # Skip locations with too few observations (here we require at least max_ma + 5 observations)
  if (length(Q_ts) < (max_ma + 5)) next
  
  # Loop over different MA orders (from 1 to max_ma)
  for (q in 1:max_ma) {
    
    # ---- 5-Fold Cross-Validation for MA(q) ----
    folds <- createFolds(Q_ts, k = 5, list = TRUE, returnTrain = TRUE)
    rmse_folds <- c()
    
    for (fold in folds) {
      train_idx <- fold
      test_idx <- setdiff(seq_along(Q_ts), train_idx)
      
      train_series <- Q_ts[train_idx]
      test_series <- Q_ts[test_idx]
      
      # Ensure the training series is long enough for an MA model of order 'q'
      if (length(train_series) <= q) next
      
      # Fit an MA model using arima (order = c(0, 0, q))
      ma_model <- try(arima(train_series, order = c(0, 0, q)), silent = TRUE)
      if (inherits(ma_model, "try-error")) next
      
      # Predict the next n values (n = length(test_series))
      preds <- try(predict(ma_model, n.ahead = length(test_series))$pred, silent = TRUE)
      if (inherits(preds, "try-error")) next
      
      # Compute RMSE for the current fold
      fold_rmse <- sqrt(mean((test_series - preds)^2, na.rm = TRUE))
      rmse_folds <- c(rmse_folds, fold_rmse)
    }
    
    # Average the RMSE across the folds for this MA order
    avg_rmse <- mean(rmse_folds, na.rm = TRUE)
    
    # ---- Full Model Fit on Entire Series ----
    # Fit the MA model on the full series to extract AIC and BIC values
    full_model <- try(arima(Q_ts, order = c(0, 0, q)), silent = TRUE)
    if (inherits(full_model, "try-error")) next
    model_aic <- AIC(full_model)
    model_bic <- BIC(full_model)
    
    # Append these results to our results table
    ma_results <- rbind(ma_results, data.frame(Location = location_name,
                                               MA_Order = q,
                                               AIC = model_aic,
                                               BIC = model_bic,
                                               RMSE = avg_rmse,
                                               stringsAsFactors = FALSE))
  }
}

# Print the final results table
print(ma_results)

# ---- Plotting ----
# Plot AIC vs. MA Order
plot_ma_aic <- ggplot(ma_results, aes(x = MA_Order, y = AIC, color = Location)) +
  geom_line() +
  geom_point() +
  labs(title = "AIC vs. MA Order", x = "MA Order", y = "AIC") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal()

# Plot BIC vs. MA Order
plot_ma_bic <- ggplot(ma_results, aes(x = MA_Order, y = BIC, color = Location)) +
  geom_line() +
  geom_point() +
  labs(title = "BIC vs. MA Order", x = "MA Order", y = "BIC") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal()

# Plot RMSE vs. MA Order (evaluated via 5-Fold CV)
plot_ma_rmse <- ggplot(ma_results, aes(x = MA_Order, y = RMSE, color = Location)) +
  geom_line() +
  geom_point() +
  labs(title = "RMSE vs. MA Order (5-Fold CV)", x = "MA Order", y = "RMSE") +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal()

# Display the plots
print(plot_ma_aic)
print(plot_ma_bic)
print(plot_ma_rmse)


# PLOT AR----
# Create an empty results table (only AIC and RMSE)
results <- data.frame(
  Location = character(),
  AR_Lag = integer(),
  AIC = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

for (location_name in location_names) {
  
  # Get the location code from the locations data
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if not found
  
  # Fetch time series data for parameter Q for the current location
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(q_data_response)) next  # Skip if data fetch fails
  
  # Process and preprocess the data
  data <- process_and_append_data(q_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next
  
  data <- preprocess_data(data)
  
  # Keep only the date and Q columns and sort by date (ascending)
  data <- data %>%
    select(data, Q) %>%
    arrange(data)
  
  # Extract the Q values as our time series vector and remove missing values
  Q_ts <- na.omit(data$Q)
  
  # Skip locations with too few observations (here we require at least max_lag + 5 observations)
  if (length(Q_ts) < (max_lag + 5)) next
  
  # Loop over different AR lags (from 1 to max_lag)
  for (lag in 1:max_lag) {
    
    # ---- 5-Fold Cross-Validation for AR(lag) ----
    folds <- createFolds(Q_ts, k = 5, list = TRUE, returnTrain = TRUE)
    rmse_folds <- c()
    
    for (fold in folds) {
      train_idx <- fold
      test_idx <- setdiff(seq_along(Q_ts), train_idx)
      
      train_series <- Q_ts[train_idx]
      test_series <- Q_ts[test_idx]
      
      # Ensure the training series is long enough for an AR model of order 'lag'
      if (length(train_series) <= lag) next
      
      # Fit an AR model (ARIMA with no differencing or MA part)
      ar_model <- try(arima(train_series, order = c(lag, 0, 0)), silent = TRUE)
      if (inherits(ar_model, "try-error")) next
      
      # Predict the next n values (n = length(test_series))
      preds <- try(predict(ar_model, n.ahead = length(test_series))$pred, silent = TRUE)
      if (inherits(preds, "try-error")) next
      
      # Compute RMSE for the current fold
      fold_rmse <- sqrt(mean((test_series - preds)^2, na.rm = TRUE))
      rmse_folds <- c(rmse_folds, fold_rmse)
    }
    
    # Average the RMSE across the folds for this AR lag
    avg_rmse <- mean(rmse_folds, na.rm = TRUE)
    
    # ---- Full Model Fit on Entire Series ----
    full_model <- try(arima(Q_ts, order = c(lag, 0, 0)), silent = TRUE)
    if (inherits(full_model, "try-error")) next
    model_aic <- AIC(full_model)
    
    # Append the results (only AIC and RMSE)
    results <- rbind(results, data.frame(
      Location = location_name,
      AR_Lag = lag,
      AIC = model_aic,
      RMSE = avg_rmse,
      stringsAsFactors = FALSE
    ))
  }
}

# Print the final results table
print(results)

# ---- Combine Metrics and Plot ----
# Reshape the data so that AIC and RMSE appear in one column for faceting
results_melt <- pivot_longer(results, cols = c(AIC, RMSE), 
                             names_to = "Metric", values_to = "Value")

# Create a combined plot with facets for AIC and RMSE
plot_combined <- ggplot(results_melt, aes(x = AR_Lag, y = Value, color = Location)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 1:10) +
  labs(
       x = "AR Lag",
       y = "Metric Value") +
  theme_minimal() +
  theme(legend.position = "bottom")  # single legend at the bottom

# Display the combined plot
print(plot_combined)


# PLOT MA ----
# Set maximum MA order to test
max_ma <- 10

# Create an empty results table to store MA model metrics for each location and order (AIC and RMSE only)
ma_results <- data.frame(
  Location = character(),
  MA_Order = integer(),
  AIC = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

for (location_name in location_names) {
  
  # Get the location code from the locations data
  location_code <- locations_df[locations_df$name == location_name, "code"]
  if (length(location_code) == 0) next  # Skip if location not found
  
  # Fetch time series data (for parameter Q) for the current location
  ma_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date, to_date)
  if (is.null(ma_data_response)) next  # Skip if data fetch fails
  
  # Process and preprocess the data
  data <- process_and_append_data(ma_data_response, parameter, location_name, data.frame())
  if (nrow(data) == 0) next
  
  data <- preprocess_data(data)
  
  # Keep only the date and Q columns and sort by date (ascending)
  data <- data %>%
    select(data, Q) %>%
    arrange(data)
  
  # Extract the Q values as our time series vector and remove any missing values
  Q_ts <- na.omit(data$Q)
  
  # Skip locations with too few observations (here we require at least max_ma + 5 observations)
  if (length(Q_ts) < (max_ma + 5)) next
  
  # Loop over different MA orders (from 1 to max_ma)
  for (q in 1:max_ma) {
    
    # ---- 5-Fold Cross-Validation for MA(q) ----
    folds <- createFolds(Q_ts, k = 5, list = TRUE, returnTrain = TRUE)
    rmse_folds <- c()
    
    for (fold in folds) {
      train_idx <- fold
      test_idx <- setdiff(seq_along(Q_ts), train_idx)
      
      train_series <- Q_ts[train_idx]
      test_series <- Q_ts[test_idx]
      
      # Ensure the training series is long enough for an MA model of order 'q'
      if (length(train_series) <= q) next
      
      # Fit an MA model using arima (order = c(0, 0, q))
      ma_model <- try(arima(train_series, order = c(0, 0, q)), silent = TRUE)
      if (inherits(ma_model, "try-error")) next
      
      # Predict the next n values (n = length(test_series))
      preds <- try(predict(ma_model, n.ahead = length(test_series))$pred, silent = TRUE)
      if (inherits(preds, "try-error")) next
      
      # Compute RMSE for the current fold
      fold_rmse <- sqrt(mean((test_series - preds)^2, na.rm = TRUE))
      rmse_folds <- c(rmse_folds, fold_rmse)
    }
    
    # Average the RMSE across the folds for this MA order
    avg_rmse <- mean(rmse_folds, na.rm = TRUE)
    
    # ---- Full Model Fit on Entire Series ----
    # Fit the MA model on the full series to extract AIC
    full_model <- try(arima(Q_ts, order = c(0, 0, q)), silent = TRUE)
    if (inherits(full_model, "try-error")) next
    model_aic <- AIC(full_model)
    
    # Append these results to our results table
    ma_results <- rbind(ma_results, data.frame(
      Location = location_name,
      MA_Order = q,
      AIC = model_aic,
      RMSE = avg_rmse,
      stringsAsFactors = FALSE
    ))
  }
}

# Print the final results table
print(ma_results)

# ---- Combine Metrics and Plot ----
# Reshape the data so that AIC and RMSE appear in one column for faceting
ma_results_melt <- pivot_longer(ma_results, cols = c(AIC, RMSE), 
                                names_to = "Metric", values_to = "Value")

# Create a combined plot with facets for AIC and RMSE
plot_combined_ma <- ggplot(ma_results_melt, aes(x = MA_Order, y = Value, color = Location)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 1:10) +
  labs(
       x = "MA Order",
       y = "Metric Value") +
  theme_minimal() +
  theme(legend.position = "bottom")  

# Display the combined plot
print(plot_combined_ma)

