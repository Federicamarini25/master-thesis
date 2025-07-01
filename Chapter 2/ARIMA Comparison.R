# SECOND ARIMA MODEL ----
domain <- "surfacewater"
parameter <- "Q"
from_date2 <- as.Date("2014-01-01") # 2013??
to_date2 <- as.Date("2023-12-31")

locations_df <- fetch_locations_data(domain)
if (is.null(locations_df) || !"name" %in% names(locations_df)) {
  stop("Failed to fetch location data.")
}

# Initialize a data frame to store the best ARIMA models
arima_results2 <- data.frame(
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
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date2, to_date2)
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
  
  # Fit ARIMA model
  arima_model <- auto.arima(data$Q)
  
  # Extract ARIMA order (p, d, q)
  arima_order <- arimaorder(arima_model)
  
  # Append results to the data frame
  arima_results2 <- rbind(arima_results2, data.frame(
    Location = location_name,
    p = arima_order[1],
    d = arima_order[2],
    q = arima_order[3]
  ))
}

# Print the final ARIMA model table
print(arima_results2)

# THIRD ARIMA MODEL ----
domain <- "surfacewater"
parameter <- "Q"
from_date3 <- as.Date("2017-01-01") # 2013??
to_date3 <- as.Date("2025-01-01")

locations_df <- fetch_locations_data(domain)
if (is.null(locations_df) || !"name" %in% names(locations_df)) {
  stop("Failed to fetch location data.")
}

# Initialize a data frame to store the best ARIMA models
arima_results3 <- data.frame(
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
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "d", from_date3, to_date3)
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
  
  # Fit ARIMA model
  arima_model <- auto.arima(data$Q)
  
  # Extract ARIMA order (p, d, q)
  arima_order <- arimaorder(arima_model)
  
  # Append results to the data frame
  arima_results3 <- rbind(arima_results3, data.frame(
    Location = location_name,
    p = arima_order[1],
    d = arima_order[2],
    q = arima_order[3]
  ))
}

# Print the final ARIMA model table
print(arima_results3)
