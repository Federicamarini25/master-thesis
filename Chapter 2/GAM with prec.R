# Load necessary libraries
library(dplyr)
library(mgcv)
library(caret)
library(future)
library(future.apply)
library(tidyr)

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


# PREPARE THE DATASET----
# Reshape precipitation data: wide format (date-wise matching)
prec_wide <- all_prec_data %>%
  select(data, location_name, prec) %>%
  pivot_wider(names_from = location_name, values_from = prec)

prec_wide <- prec_wide %>%
  rename_with(~ gsub("[-.,]", "", gsub(" ", "_", .x)))  # Replace spaces with underscores


# Merge river flow data (`Q`) with precipitation (`Prec`) data by date
model_data <- left_join(all_q_data, prec_wide, by = "data")

# Remove all rows in which Q is NA
model_data <- model_data %>% filter(!is.na(Q))

# Remove columns made only by NA
model_data <- model_data %>% 
  select(where(~ !all(is.na(.))))

colSums(is.na(model_data))

predictor_cols <- names(model_data)[6:ncol(model_data)]  # Select from column 6 onward

# Select only relevant columns (Q as dependent variable, precipitation as predictors)
model_data <- model_data %>%
  mutate(time_num = as.numeric(data)) %>%  # Convert 'data' column to numeric
  select(Q, time_num, location_name, all_of(predictor_cols))  # Select relevant columns


# MODEL 1 ----
# Helper function to clean and validate location-specific data
clean_location_data <- function(data, loc) {
  message(paste("Cleaning data for location:", loc))
  
  # Remove columns with any NA values
  data_clean <- data %>% select(where(~ all(!is.na(.))))
  
  # Remove columns with constant values (only one unique value)
  constant_cols <- names(which(sapply(data_clean, function(col) length(unique(col)) == 1)))
  if (length(constant_cols) > 0) {
    message(paste("  Removing constant columns:", paste(constant_cols, collapse = ", ")))
    data_clean <- data_clean %>% select(-all_of(constant_cols))
  }
  
  # Check that 'Q' exists in the dataset; if not, return NULL
  if (!"Q" %in% names(data_clean)) {
    message(paste("  Skipping", loc, "- 'Q' variable missing"))
    return(NULL)
  }
  
  # Identify and remove non-numeric columns
  non_numeric <- names(data_clean)[!sapply(data_clean, is.numeric)]
  if (length(non_numeric) > 0) {
    message(paste("  Skipping", loc, "- non-numeric columns present:", paste(non_numeric, collapse = ", ")))
    return(NULL)
  }
  
  # Check if enough observations exist
  if (nrow(data_clean) < 10) {
    message(paste("  Skipping", loc, "- insufficient data (<10 observations)"))
    return(NULL)
  }
  
  return(data_clean)
}

# --- GAM MODEL WITH JUST PRECIPITATION PREDICTORS ---
# Initialize lists to store fitted GAM models and any errors encountered
gam_models1 <- list()
gam_errors1 <- list()

# Get unique location names
unique_locations <- unique(model_data$location_name)

for (loc in unique_locations) {
  message("------------------------------------------------------")
  message(paste("Processing location:", loc))
  
  # Subset data for the current location
  loc_data <- model_data %>% filter(location_name == loc)
  
  # Clean the data using the helper function
  loc_data <- clean_location_data(loc_data, loc)
  if (is.null(loc_data)) {
    message(paste("Data cleaning failed for", loc, "- model fitting skipped."))
    next
  }
  
  # Identify predictor columns: exclude 'Q', 'location_name', and 'time_num'
  predictor_cols <- setdiff(names(loc_data), c("Q", "location_name", "time_num"))
  if (length(predictor_cols) == 0) {
    message(paste("Skipping", loc, "- no valid predictors remain"))
    next
  }
  
  # Dynamically build the GAM formula using only precipitation predictors.
  # Using cubic regression splines ("cr") with k = 5 basis functions.
  gam_formula_str <- paste("Q ~", paste(sprintf("s(%s, bs = 'cr', k = 5)", predictor_cols), collapse = " + "))
  gam_formula <- as.formula(gam_formula_str)
  
  # Instead of printing the full formula, print a simpler confirmation
  message("GAM formula successfully built for precipitation predictors.")
  
  # Fit GAM model with error handling
  gam_fit <- tryCatch({
    gam(gam_formula, data = loc_data, method = "GCV.Cp", select = TRUE)
  }, error = function(e) {
    message(paste("Error in location", loc, ":", e$message))
    gam_errors1[[loc]] <<- e$message  # store error message for future reference
    return(NULL)
  })
  
  # Store the model if fitting is successful
  if (!is.null(gam_fit)) {
    gam_models1[[loc]] <- gam_fit
    message(paste("Model successfully fitted for:", loc))
  }
}

# Summary of the first successfully fitted model
if (length(gam_models1) > 0) {
  first_model <- gam_models1[[1]]
  cat("------------------------------------------------------\n")
  cat("Summary of the first fitted GAM model:\n")
  print(summary(first_model))
} else {
  cat("No GAM models were successfully fitted.\n")
} # 0.539

# MODEL 1 CV ----
# --- Set Up Parallel Processing ---
plan(multisession)

# --- Initialize Storage for Results ---
rmse_values1  <- c()    # To store averaged RMSE for each location
q_ranges1     <- c()    # To store the range difference of Q per location
log_messages1 <- c()    # To store log messages during the simulation

set.seed(1234)  # For reproducibility

# --- Loop Over Each Location ---
for (loc in unique(model_data$location_name)) {
  
  log_messages1 <- c(log_messages1, paste("Processing location:", loc))
  
  # Subset the data for the current location
  data_loc <- model_data %>% filter(location_name == loc)
  
  # --- Remove Columns with Any Missing Values (using tidyverse) ---
  data_loc <- data_loc %>% select(where(~ !any(is.na(.))))
  
  # --- Remove Constant Predictors (columns with 1 unique value) ---
  constant_cols <- names(which(sapply(data_loc, function(x) length(unique(x)) <= 1)))
  if (length(constant_cols) > 0) {
    data_loc <- data_loc %>% select(-all_of(constant_cols))
  }
  
  # --- Ensure 'Q' Exists ---
  if (!"Q" %in% names(data_loc)) {
    log_messages1 <- c(log_messages1, paste("Skipping", loc, "- Q variable missing"))
    next
  }
  
  # --- Simulate 20% Missing Values in Q ---
  set.seed(1234)
  na_indices <- sample(seq_len(nrow(data_loc)), size = floor(0.2 * nrow(data_loc)))
  data_loc$Q_missing <- data_loc$Q  # Retain original Q values
  data_loc$Q[na_indices] <- NA      # Introduce 20% missing values
  
  # --- Identify Numeric Predictor Columns ---
  # Exclude Q, Q_missing, location_name and time
  all_predictors <- setdiff(names(data_loc), c("Q", "Q_missing", "location_name", "time_num"))
  numeric_predictors <- all_predictors[sapply(data_loc[all_predictors], is.numeric)]
  
  # --- Build GAM Formula Without a Time Term ---
  if (length(numeric_predictors) == 0) {
    # Use an intercept-only model if no additional numeric predictors exist.
    gam_formula <- as.formula("Q ~ 1")
    log_messages1 <- c(log_messages1, paste("Location", loc, ": No additional numeric predictors; using intercept-only model"))
  } else {
    smooth_terms <- paste(sapply(numeric_predictors, function(pred) {
      paste0("s(", pred, ", bs = 'cr', k = 5)")
    }), collapse = " + ")
    gam_formula <- as.formula(paste("Q ~", smooth_terms))
    # Log a simple confirmation without printing the full formula.
    log_messages1 <- c(log_messages1, paste("Location", loc, ": GAM formula built using precipitation predictors"))
  }
  
  # --- Ensure Sufficient Data ---
  if (nrow(data_loc) < 10) {
    log_messages1 <- c(log_messages1, paste("Skipping", loc, "- insufficient data"))
    next
  }
  
  # --- Compute Q Range Difference ---
  q_range <- diff(range(data_loc$Q_missing, na.rm = TRUE))
  q_ranges1 <- c(q_ranges1, q_range)
  
  # --- Fit GAM Model for Imputation ---
  gam_model <- tryCatch({
    gam(gam_formula, data = data_loc, method = "GCV.Cp", select = TRUE)
  }, error = function(e) {
    message(paste("Error fitting model for", loc, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(gam_model)) {
    # Impute missing Q values using the fitted model.
    missing_idx <- which(is.na(data_loc$Q))
    if (length(missing_idx) > 0) {
      data_loc$Q[missing_idx] <- predict(gam_model, newdata = data_loc[missing_idx, ])
    }
    log_messages1 <- c(log_messages1, paste("Location:", loc, "- Imputation completed"))
  } else {
    log_messages1 <- c(log_messages1, paste("Skipping", loc, "- Model fitting failed"))
    next
  }
  
  # --- 5-FOLD CROSS-VALIDATION ---
  folds <- createFolds(data_loc$Q_missing, k = 5, list = TRUE, returnTrain = TRUE)
  
  rmse_folds <- future_sapply(folds, function(train_idx) {
    test_idx <- setdiff(seq_len(nrow(data_loc)), train_idx)
    train_data <- data_loc[train_idx, ]
    test_data  <- data_loc[test_idx, ]
    
    gam_model_cv <- tryCatch({
      gam(gam_formula, data = train_data, method = "GCV.Cp", select = TRUE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(gam_model_cv)) return(NA)
    
    predictions <- predict(gam_model_cv, newdata = test_data)
    sqrt(mean((test_data$Q_missing - predictions)^2, na.rm = TRUE))
  }, future.seed = 1234)
  
  rmse_loc <- mean(rmse_folds, na.rm = TRUE)
  rmse_values1 <- c(rmse_values1, rmse_loc)
  
  log_messages1 <- c(log_messages1, paste("Location:", loc, "- GAM RMSE:", round(rmse_loc, 3),
                                        "- Q Range Diff:", round(q_range, 3)))
}

# --- Final Log Output ---
cat(log_messages1, sep = "\n") 



# MODEL 2 ----
# Initialize lists to store fitted GAM models and any errors encountered
gam_models2 <- list()
gam_errors2 <- list()

# Get unique location names
unique_locations <- unique(model_data$location_name)

for (loc in unique_locations) {
  message("------------------------------------------------------")
  message(paste("Processing location:", loc))
  
  # Subset data for the current location
  loc_data <- model_data %>% filter(location_name == loc)
  
  # Clean the data using the helper function (assumes clean_location_data is defined)
  loc_data <- clean_location_data(loc_data, loc)
  if (is.null(loc_data)) {
    message(paste("Data cleaning failed for", loc, "- model fitting skipped."))
    next
  }
  
  # Identify predictor columns: exclude 'Q' and 'location_name' (keep time_num)
  predictor_cols <- setdiff(names(loc_data), c("Q", "location_name"))
  if (length(predictor_cols) == 0) {
    message(paste("Skipping", loc, "- no valid predictors remain"))
    next
  }
  
  # Dynamically build the GAM formula using all numeric predictors (including time_num)
  # Using cubic regression splines ("cr") with k = 5 basis functions.
  gam_formula_str <- paste("Q ~", paste(sprintf("s(%s, bs = 'cr', k = 5)", predictor_cols), collapse = " + "))
  gam_formula <- as.formula(gam_formula_str)
  
  # Instead of printing the full formula, log a simple confirmation
  message("GAM formula successfully built including time_num and precipitation predictors.")
  
  # Fit GAM model with error handling
  gam_fit <- tryCatch({
    gam(gam_formula, data = loc_data, method = "GCV.Cp", select = TRUE)
  }, error = function(e) {
    message(paste("Error in location", loc, ":", e$message))
    gam_errors2[[loc]] <- e$message  # store error message for future reference
    return(NULL)
  })
  
  # Store the model if fitting is successful
  if (!is.null(gam_fit)) {
    gam_models2[[loc]] <- gam_fit
    message(paste("Model successfully fitted for:", loc))
  }
}

# Summary of the first successfully fitted model
if (length(gam_models2) > 0) {
  first_model <- gam_models2[[1]]
  cat("------------------------------------------------------\n")
  cat("Summary of the first fitted GAM model (with time_num):\n")
  print(summary(first_model))
} else {
  cat("No GAM models were successfully fitted.\n")
}
# 0.572 

# MODEL 2 (first location only) ----
# Initialize storage
gam_models2_first_loc <- list()
gam_errors2_first_loc <- list()

# Get the first location name
first_loc <- unique(model_data$location_name)[1]

message("------------------------------------------------------")
message(paste("Processing first location:", first_loc))

# Subset data for the first location
loc_data <- model_data %>% filter(location_name == first_loc)

# Clean the data
loc_data <- clean_location_data(loc_data, first_loc)
if (is.null(loc_data)) {
  message(paste("Data cleaning failed for", first_loc, "- model fitting skipped."))
} else {
  
  # Identify predictors: exclude Q and location_name (keep time_num)
  predictor_cols <- setdiff(names(loc_data), c("Q", "location_name"))
  
  if (length(predictor_cols) == 0) {
    message(paste("Skipping", first_loc, "- no valid predictors remain"))
  } else {
    
    # Build GAM formula
    gam_formula_str <- paste("Q ~", paste(sprintf("s(%s, bs = 'cr', k = 5)", predictor_cols), collapse = " + "))
    gam_formula <- as.formula(gam_formula_str)
    
    message("GAM formula successfully built for first location (with time_num and precipitation predictors).")
    
    # Fit the model with error handling
    gam_fit <- tryCatch({
      gam(gam_formula, data = loc_data, method = "GCV.Cp", select = TRUE)
    }, error = function(e) {
      message(paste("Error fitting model for", first_loc, ":", e$message))
      gam_errors2_first_loc[[first_loc]] <- e$message
      return(NULL)
    })
    
    if (!is.null(gam_fit)) {
      gam_models2_first_loc[[first_loc]] <- gam_fit
      message(paste("Model successfully fitted for:", first_loc))
    }
  }
}

# View model summary
if (!is.null(gam_models2_first_loc[[first_loc]])) {
  cat("------------------------------------------------------\n")
  cat("Summary of the GAM model for first location:\n")
  print(summary(gam_models2_first_loc[[first_loc]]))
} else {
  cat("No model was successfully fitted for the first location.\n")
} # 0.572


# CV MODEL 2 ----
# --- Set Up Parallel Processing ---
plan(multisession)

# --- Initialize Storage for Results ---
rmse_values2  <- c()
q_ranges2     <- c()
log_messages2 <- c()

set.seed(1234)  # Reproducibility

# --- Loop Over Each Location ---
for (loc in unique(model_data$location_name)) {
  
  log_messages2 <- c(log_messages2, paste("Processing location:", loc))
  
  # Subset data for the current location
  data_loc <- model_data %>% filter(location_name == loc)
  
  # --- Remove Columns with Any Missing Values ---
  data_loc <- data_loc %>% select(where(~ !any(is.na(.))))
  
  # --- Remove Constant Columns ---
  constant_cols <- names(which(sapply(data_loc, function(x) length(unique(x)) <= 1)))
  if (length(constant_cols) > 0) {
    data_loc <- data_loc %>% select(-all_of(constant_cols))
  }
  
  # --- Ensure 'Q' Exists ---
  if (!"Q" %in% names(data_loc)) {
    log_messages2 <- c(log_messages2, paste("Skipping", loc, "- Q variable missing"))
    next
  }
  
  # --- Simulate 20% Missing in Q ---
  set.seed(1234)
  na_indices <- sample(seq_len(nrow(data_loc)), size = floor(0.2 * nrow(data_loc)))
  data_loc$Q_missing <- data_loc$Q
  data_loc$Q[na_indices] <- NA
  
  # --- Identify Predictors (including time_num) ---
  all_predictors <- setdiff(names(data_loc), c("Q", "Q_missing", "location_name"))
  numeric_predictors <- all_predictors[sapply(data_loc[all_predictors], is.numeric)]
  
  # --- Build GAM Formula (Model 2 style: includes time_num) ---
  if (length(numeric_predictors) == 0) {
    gam_formula <- as.formula("Q ~ 1")
    log_messages2 <- c(log_messages2, paste("Location", loc, ": No predictors; using intercept-only model"))
  } else {
    smooth_terms <- paste0("s(", numeric_predictors, ", bs = 'cr', k = 5)", collapse = " + ")
    gam_formula <- as.formula(paste("Q ~", smooth_terms))
    log_messages2 <- c(log_messages2, paste("Location", loc, ": GAM formula (with time_num) built."))
  }
  
  # --- Skip if data is too small ---
  if (nrow(data_loc) < 10) {
    log_messages2 <- c(log_messages2, paste("Skipping", loc, "- insufficient data"))
    next
  }
  
  # --- Compute Q Range ---
  q_range <- diff(range(data_loc$Q_missing, na.rm = TRUE))
  q_ranges2 <- c(q_ranges2, q_range)
  
  # --- Fit Model for Imputation ---
  gam_model <- tryCatch({
    gam(gam_formula, data = data_loc, method = "GCV.Cp", select = TRUE)
  }, error = function(e) {
    message(paste("Error fitting model for", loc, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(gam_model)) {
    missing_idx <- which(is.na(data_loc$Q))
    if (length(missing_idx) > 0) {
      data_loc$Q[missing_idx] <- predict(gam_model, newdata = data_loc[missing_idx, ])
    }
    log_messages2 <- c(log_messages2, paste("Location:", loc, "- Imputation completed"))
  } else {
    log_messages2 <- c(log_messages2, paste("Skipping", loc, "- Model fitting failed"))
    next
  }
  
  # --- 5-Fold Cross-Validation ---
  folds <- createFolds(data_loc$Q_missing, k = 5, list = TRUE, returnTrain = TRUE)
  
  rmse_folds <- future_sapply(folds, function(train_idx) {
    test_idx <- setdiff(seq_len(nrow(data_loc)), train_idx)
    train_data <- data_loc[train_idx, ]
    test_data  <- data_loc[test_idx, ]
    
    gam_model_cv <- tryCatch({
      gam(gam_formula, data = train_data, method = "GCV.Cp", select = TRUE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(gam_model_cv)) return(NA)
    
    predictions <- predict(gam_model_cv, newdata = test_data)
    sqrt(mean((test_data$Q_missing - predictions)^2, na.rm = TRUE))
  }, future.seed = 1234)
  
  rmse_loc <- mean(rmse_folds, na.rm = TRUE)
  rmse_values2 <- c(rmse_values2, rmse_loc)
  
  log_messages2 <- c(log_messages2, paste("Location:", loc, "- GAM RMSE:", round(rmse_loc, 3),
                                        "- Q Range Diff:", round(q_range, 3)))
}

# --- Final Log Output ---
cat(log_messages2, sep = "\n")


# PART 2: IMPUTATION AND REPEATING ----
# Store GAM models
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

# Remove columns made only by NA
all_q_data_imputed <- all_q_data_imputed %>% 
  select(where(~ !all(is.na(.))))

# Create a new dataset with linear imputation
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


prec_wide <- all_prec_data_imputed %>%
  select(data, location_name, prec) %>%
  pivot_wider(names_from = location_name, values_from = prec)

prec_wide <- prec_wide %>%
  rename_with(~ gsub("[-.,]", "", gsub(" ", "_", .x)))  # Replace spaces with underscores


# Merge river flow data (`Q`) with precipitation (`Prec`) data by date
model_data <- left_join(all_q_data_imputed, prec_wide, by = "data")

# Remove all rows in which Q is NA
model_data <- model_data %>% filter(!is.na(Q))
predictor_cols <- names(model_data)[6:ncol(model_data)]  # Select from column 6 onward


model_data <- model_data %>%
  mutate(time_num = as.numeric(data)) %>%  # Convert 'data' column to numeric
  select(Q, time_num, location_name, all_of(predictor_cols))  # Select relevant columns


# Select only relevant columns (Q as dependent variable, precipitation as predictors)
model_data_imputed <- model_data %>% select(Q, time_num, location_name, all_of(predictor_cols))
colSums(is.na(model_data_imputed))


# MODEL 3 ---- 
# --- GAM MODEL WITH JUST PRECIPITATION PREDICTORS ---
# (No cleaning is needed, so we directly work on model_data_imputed)
# Initialize lists to store fitted GAM models and any errors encountered
gam_models3 <- list()
gam_errors3 <- list()

# Get unique location names from the new imputed dataset
unique_locations <- unique(model_data_imputed$location_name)

for (loc in unique_locations) {
  message("------------------------------------------------------")
  message(paste("Processing location:", loc))
  
  # Subset data for the current location from the fully imputed dataset
  loc_data <- model_data_imputed %>% filter(location_name == loc)
  
  # Check if there are enough observations
  if (nrow(loc_data) < 10) {
    message(paste("Skipping", loc, "- insufficient data (<10 observations)"))
    next
  }
  
  # Identify predictor columns: exclude 'Q' and 'location_name'
  predictor_cols <- setdiff(names(loc_data), c("Q", "location_name"))
  if (length(predictor_cols) == 0) {
    message(paste("Skipping", loc, "- no valid predictors remain"))
    next
  }
  
  # Build the GAM formula dynamically using precipitation predictors
  # Here we use cubic regression splines ("cr") with k = 5 basis functions for each predictor.
  gam_formula_str <- paste("Q ~", paste(sprintf("s(%s, bs = 'cr', k = 5)", predictor_cols), collapse = " + "))
  gam_formula <- as.formula(gam_formula_str)
  
  message("GAM formula successfully built for precipitation predictors.")
  
  # Fit the GAM model with error handling using the preprocessed data
  gam_fit <- tryCatch({
    gam(gam_formula, data = loc_data, method = "GCV.Cp", select = TRUE)
  }, error = function(e) {
    message(paste("Error in location", loc, ":", e$message))
    gam_errors3[[loc]] <<- e$message  # store error message for future reference
    return(NULL)
  })
  
  if (!is.null(gam_fit)) {
    gam_models3[[loc]] <- gam_fit
    message(paste("Model successfully fitted for:", loc))
  }
}

# Summary of the first successfully fitted model
if (length(gam_models3) > 0) {
  first_model <- gam_models3[[1]]
  cat("------------------------------------------------------\n")
  cat("Summary of the first fitted GAM model:\n")
  print(summary(first_model))
} else {
  cat("No GAM models were successfully fitted.\n")
}

# JUST FOR THE FIRST LOCATION ----
# Initialize storage for the first location
gam_models3_first_loc <- list()
gam_errors3_first_loc <- list()

# Get the first location name
first_loc <- unique(model_data_imputed$location_name)[1]

message("------------------------------------------------------")
message(paste("Processing first location:", first_loc))

# Subset data for the first location from the fully imputed dataset
loc_data <- model_data_imputed %>% filter(location_name == first_loc)

# Check for minimum number of observations
if (nrow(loc_data) < 10) {
  message(paste("Skipping", first_loc, "- insufficient data (<10 observations)"))
} else {
  
  # Identify predictors: exclude Q and location_name (retain time_num if present)
  predictor_cols <- setdiff(names(loc_data), c("Q", "location_name", "time_num"))
  
  if (length(predictor_cols) == 0) {
    message(paste("Skipping", first_loc, "- no valid predictors remain"))
  } else {
    # Build the GAM formula
    gam_formula_str <- paste("Q ~", paste(sprintf("s(%s, bs = 'cr', k = 5)", predictor_cols),
                                          collapse = " + "))
    gam_formula <- as.formula(gam_formula_str)
    
    message("GAM formula successfully built for first location (with time_num and precipitation predictors).")
    
    # Fit the model with error handling
    gam_fit <- tryCatch({
      gam(gam_formula, data = loc_data, method = "GCV.Cp", select = TRUE)
    }, error = function(e) {
      message(paste("Error fitting model for", first_loc, ":", e$message))
      gam_errors3_first_loc[[first_loc]] <- e$message
      return(NULL)
    })
    
    if (!is.null(gam_fit)) {
      gam_models3_first_loc[[first_loc]] <- gam_fit
      message(paste("Model successfully fitted for:", first_loc))
    }
  }
}

# View the model summary for the first location
if (!is.null(gam_models3_first_loc[[first_loc]])) {
  cat("------------------------------------------------------\n")
  cat("Summary of the GAM model for first location:\n")
  print(summary(gam_models3_first_loc[[first_loc]]))
} else {
  cat("No model was successfully fitted for the first location.\n")
} # 0.671


# --- MODEL 3 CROSS-VALIDATION (CV) ---
# Set up parallel processing
plan(multisession)

# Initialize storage for results
rmse_values3  <- c()    # To store averaged RMSE for each location
q_ranges3     <- c()    # To store the range difference of Q per location
log_messages3 <- c()    # To store log messages during the simulation

set.seed(1234)  # For reproducibility

for (loc in unique(model_data_imputed$location_name)) {
  
  log_messages3 <- c(log_messages3, paste("Processing location:", loc))
  
  # Subset data for the current location from the new imputed dataset
  data_loc <- model_data_imputed %>% filter(location_name == loc)
  
  # Ensure that there are enough rows
  if (nrow(data_loc) < 10) {
    log_messages3 <- c(log_messages3, paste("Skipping", loc, "- insufficient data"))
    next
  }
  
  # Simulate 20% missing values in Q while keeping the original values in a backup column
  set.seed(1234)  # Reset seed for reproducibility; adjust if needed for independent missing patterns
  na_indices <- sample(seq_len(nrow(data_loc)), size = floor(0.2 * nrow(data_loc)))
  data_loc$Q_missing <- data_loc$Q  # Backup original Q values
  data_loc$Q[na_indices] <- NA      # Introduce 20% missing values
  
  # Identify numeric predictor columns: exclude Q, Q_missing, and location_name
  all_predictors <- setdiff(names(data_loc), c("Q", "Q_missing", "location_name", "time_num"))
  numeric_predictors <- all_predictors[sapply(data_loc[all_predictors], is.numeric)]
  
  # Build the GAM formula using the available numeric predictors
  if (length(numeric_predictors) == 0) {
    gam_formula <- as.formula("Q ~ 1")
    log_messages3 <- c(log_messages3, paste("Location", loc, ": No additional numeric predictors; using intercept-only model"))
  } else {
    smooth_terms <- paste(sapply(numeric_predictors, function(pred) {
      paste0("s(", pred, ", bs = 'cr', k = 5)")
    }), collapse = " + ")
    gam_formula <- as.formula(paste("Q ~", smooth_terms))
    log_messages3 <- c(log_messages3, paste("Location", loc, ": GAM formula built using precipitation predictors"))
  }
  
  # Compute the range of the original Q values
  q_range <- diff(range(data_loc$Q_missing, na.rm = TRUE))
  q_ranges3 <- c(q_ranges3, q_range)
  
  # Fit GAM model for imputation on the dataset with simulated missing values
  gam_model <- tryCatch({
    gam(gam_formula, data = data_loc, method = "GCV.Cp", select = TRUE)
  }, error = function(e) {
    message(paste("Error fitting model for", loc, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(gam_model)) {
    # Impute missing Q values using predictions from the fitted model
    missing_idx <- which(is.na(data_loc$Q))
    if (length(missing_idx) > 0) {
      data_loc$Q[missing_idx] <- predict(gam_model, newdata = data_loc[missing_idx, ])
    }
    log_messages3 <- c(log_messages3, paste("Location:", loc, "- Imputation completed"))
  } else {
    log_messages3 <- c(log_messages3, paste("Skipping", loc, "- Model fitting failed"))
    next
  }
  
  # --- 5-FOLD CROSS-VALIDATION ---
  folds <- createFolds(data_loc$Q_missing, k = 5, list = TRUE, returnTrain = TRUE)
  
  rmse_folds <- future_sapply(folds, function(train_idx) {
    test_idx <- setdiff(seq_len(nrow(data_loc)), train_idx)
    train_data <- data_loc[train_idx, ]
    test_data  <- data_loc[test_idx, ]
    
    gam_model_cv <- tryCatch({
      gam(gam_formula, data = train_data, method = "GCV.Cp", select = TRUE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(gam_model_cv)) return(NA)
    
    predictions <- predict(gam_model_cv, newdata = test_data)
    sqrt(mean((test_data$Q_missing - predictions)^2, na.rm = TRUE))
  }, future.seed = 1234)
  
  rmse_loc <- mean(rmse_folds, na.rm = TRUE)
  rmse_values3 <- c(rmse_values3, rmse_loc)
  
  log_messages3 <- c(log_messages3, paste("Location:", loc, 
                                        "- GAM RMSE:", round(rmse_loc, 3),
                                        "- Q Range Diff:", round(q_range, 3)))
}

# Print out the log messages
cat(log_messages3, sep = "\n")


# MODEL 4 ----
# just for the first location ----
# Initialize storage for the first location
gam_models4_first_loc <- list()
gam_errors4_first_loc <- list()

# Get the first location name
first_loc <- unique(model_data_imputed$location_name)[1]

message("------------------------------------------------------")
message(paste("Processing first location:", first_loc))

# Subset data for the first location from the fully imputed dataset
loc_data <- model_data_imputed %>% filter(location_name == first_loc)

# Check for minimum number of observations
if (nrow(loc_data) < 10) {
  message(paste("Skipping", first_loc, "- insufficient data (<10 observations)"))
} else {
  
  # Identify predictors: exclude Q and location_name (retain time_num if present)
  predictor_cols <- setdiff(names(loc_data), c("Q", "location_name"))
  
  if (length(predictor_cols) == 0) {
    message(paste("Skipping", first_loc, "- no valid predictors remain"))
  } else {
    # Build the GAM formula
    gam_formula_str <- paste("Q ~", paste(sprintf("s(%s, bs = 'cr', k = 5)", predictor_cols),
                                          collapse = " + "))
    gam_formula <- as.formula(gam_formula_str)
    
    message("GAM formula successfully built for first location (with time_num and precipitation predictors).")
    
    # Fit the model with error handling
    gam_fit <- tryCatch({
      gam(gam_formula, data = loc_data, method = "GCV.Cp", select = TRUE)
    }, error = function(e) {
      message(paste("Error fitting model for", first_loc, ":", e$message))
      gam_errors4_first_loc[[first_loc]] <- e$message
      return(NULL)
    })
    
    if (!is.null(gam_fit)) {
      gam_models4_first_loc[[first_loc]] <- gam_fit
      message(paste("Model successfully fitted for:", first_loc))
    }
  }
}

# View the model summary for the first location
if (!is.null(gam_models4_first_loc[[first_loc]])) {
  cat("------------------------------------------------------\n")
  cat("Summary of the GAM model for first location:\n")
  print(summary(gam_models4_first_loc[[first_loc]]))
} else {
  cat("No model was successfully fitted for the first location.\n")
} # 0.683

# ALL THE LOCATION ----
# Initialize lists to store fitted GAM models and any errors encountered
gam_models4 <- list()
gam_errors4 <- list()

# Get unique location names from the new imputed dataset
unique_locations <- unique(model_data_imputed$location_name)

for (loc in unique_locations) {
  message("------------------------------------------------------")
  message(paste("Processing location:", loc))
  
  # Subset data for the current location from the fully imputed dataset
  loc_data <- model_data_imputed %>% filter(location_name == loc)
  
  # Check that there are enough observations
  if (nrow(loc_data) < 10) {
    message(paste("Skipping", loc, "- insufficient data (<10 observations)"))
    next
  }
  
  # Identify predictor columns: exclude 'Q' and 'location_name'
  predictor_cols <- setdiff(names(loc_data), c("Q", "location_name"))
  
  if (length(predictor_cols) == 0) {
    message(paste("Skipping", loc, "- no valid predictors remain"))
    next
  }
  
  # Build the GAM formula using all numeric predictors (including time_num if present)
  # Using cubic regression splines ("cr") with k = 5 basis functions.
  gam_formula_str <- paste("Q ~", paste(sprintf("s(%s, bs = 'cr', k = 5)", predictor_cols),
                                        collapse = " + "))
  gam_formula <- as.formula(gam_formula_str)
  
  message("GAM formula successfully built including time_num and precipitation predictors.")
  
  # Fit GAM model with error handling
  gam_fit <- tryCatch({
    gam(gam_formula, data = loc_data, method = "GCV.Cp", select = TRUE)
  }, error = function(e) {
    message(paste("Error in location", loc, ":", e$message))
    gam_errors4[[loc]] <<- e$message
    return(NULL)
  })
  
  if (!is.null(gam_fit)) {
    gam_models4[[loc]] <- gam_fit
    message(paste("Model successfully fitted for:", loc))
  }
}

# Print a summary of the first successfully fitted model
if (length(gam_models4) > 0) {
  first_model <- gam_models4[[1]]
  cat("------------------------------------------------------\n")
  cat("Summary of the first fitted GAM model (with time_num):\n")
  print(summary(first_model))
} else {
  cat("No GAM models were successfully fitted.\n")
}

# CV MODEL 4 ----
# Set up parallel processing
plan(multisession)

# Initialize storage for cross-validation results
rmse_values4  <- c()    # To store averaged RMSE for each location
q_ranges4     <- c()    # To store the range difference of Q per location
log_messages4 <- c()    # To store log messages during the simulation

set.seed(1234)  # For reproducibility

for (loc in unique(model_data_imputed$location_name)) {
  
  log_messages4 <- c(log_messages4, paste("Processing location:", loc))
  
  # Subset data for the current location from the fully imputed dataset
  data_loc <- model_data_imputed %>% filter(location_name == loc)
  
  # Ensure there are enough observations
  if (nrow(data_loc) < 10) {
    log_messages4 <- c(log_messages4, paste("Skipping", loc, "- insufficient data"))
    next
  }
  
  # Simulate 20% missing values in Q while keeping the original values in Q_missing
  set.seed(1234)
  na_indices <- sample(seq_len(nrow(data_loc)), size = floor(0.2 * nrow(data_loc)))
  data_loc$Q_missing <- data_loc$Q  
  data_loc$Q[na_indices] <- NA      
  
  # Identify numeric predictor columns: exclude Q, Q_missing, and location_name
  all_predictors <- setdiff(names(data_loc), c("Q", "Q_missing", "location_name"))
  numeric_predictors <- all_predictors[sapply(data_loc[all_predictors], is.numeric)]
  
  # Build GAM formula including all available numeric predictors
  if (length(numeric_predictors) == 0) {
    gam_formula <- as.formula("Q ~ 1")
    log_messages4 <- c(log_messages4, paste("Location", loc, ": No predictors; using intercept-only model"))
  } else {
    smooth_terms <- paste(sapply(numeric_predictors, function(pred) {
      paste0("s(", pred, ", bs = 'cr', k = 5)")
    }), collapse = " + ")
    gam_formula <- as.formula(paste("Q ~", smooth_terms))
    log_messages4 <- c(log_messages4, paste("Location", loc, ": GAM formula (with time_num) built."))
  }
  
  # Compute the range of the original Q values for later reporting
  q_range <- diff(range(data_loc$Q_missing, na.rm = TRUE))
  q_ranges4 <- c(q_ranges4, q_range)
  
  # Fit the model for imputation
  gam_model <- tryCatch({
    gam(gam_formula, data = data_loc, method = "GCV.Cp", select = TRUE)
  }, error = function(e) {
    message(paste("Error fitting model for", loc, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(gam_model)) {
    # Impute missing Q values using predictions from the fitted model.
    missing_idx <- which(is.na(data_loc$Q))
    if (length(missing_idx) > 0) {
      data_loc$Q[missing_idx] <- predict(gam_model, newdata = data_loc[missing_idx, ])
    }
    log_messages4 <- c(log_messages4, paste("Location:", loc, "- Imputation completed"))
  } else {
    log_messages4 <- c(log_messages4, paste("Skipping", loc, "- Model fitting failed"))
    next
  }
  
  # --- 5-FOLD CROSS-VALIDATION ---
  folds <- createFolds(data_loc$Q_missing, k = 5, list = TRUE, returnTrain = TRUE)
  
  rmse_folds <- future_sapply(folds, function(train_idx) {
    test_idx <- setdiff(seq_len(nrow(data_loc)), train_idx)
    train_data <- data_loc[train_idx, ]
    test_data  <- data_loc[test_idx, ]
    
    gam_model_cv <- tryCatch({
      gam(gam_formula, data = train_data, method = "GCV.Cp", select = TRUE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(gam_model_cv)) return(NA)
    
    predictions <- predict(gam_model_cv, newdata = test_data)
    sqrt(mean((test_data$Q_missing - predictions)^2, na.rm = TRUE))
  }, future.seed = 1234)
  
  rmse_loc <- mean(rmse_folds, na.rm = TRUE)
  rmse_values4 <- c(rmse_values4, rmse_loc)
  
  log_messages4 <- c(log_messages4, paste("Location:", loc, "- GAM RMSE:", round(rmse_loc, 3),
                                        "- Q Range Diff:", round(q_range, 3)))
}

# Print out the log messages from cross-validation
cat(log_messages4, sep = "\n")






# ANALYSIS OF RESIDUES ----
#  MODEL 4 (first location)
if (!is.null(gam_models4_first_loc[[first_loc]])) {
  
  message("------------------------------------------------------")
  message("Residual analysis for MODEL 4 (first location)...")
  
  model4_fit <- gam_models4_first_loc[[first_loc]]
  residuals_gam <- residuals(model4_fit)
  fitted_values_gam <- fitted(model4_fit)
  
  # 1. Residuals vs Fitted
  plot(fitted_values_gam, residuals_gam,
       main = "Residuals vs Fitted Plot",
       xlab = "Fitted", ylab = "Residuals",
       pch = 19, col = "blue",
       cex.main = 1.2, cex.lab = 0.8, cex.axis = 0.8)
  abline(h = 0, col = "red", lwd = 2)
  
  # 2. Histogram
  hist(residuals_gam, breaks = 30,
       main = "Histogram of Residuals", xlab = "Residuals",
       col = "lightblue", border = "black",
       cex.main = 1.2, cex.lab = 0.8, cex.axis = 0.8)
  
  # 3. Q-Q Plot
  qqnorm(residuals_gam, main = "Q-Q Plot of Residuals",
         cex.main = 1.2, cex.lab = 0.8, cex.axis = 0.8)
  qqline(residuals_gam, col = "red", lwd = 2)
  
  # 4. Shapiro-Wilk Test
  cat("Shapiro-Wilk Test (Model 4):\n")
  print(shapiro.test(residuals_gam))
  
  # 5. ACF & PACF
  par(mfrow = c(1,2))
  acf(residuals_gam, main = "ACF", cex.main = 0.3)
  pacf(residuals_gam, main = "PACF", cex.main = 0.3)
  par(mfrow = c(1,1))
  
} else {
  message("Residual analysis skipped - Model 4 (first location) is missing.")
}

