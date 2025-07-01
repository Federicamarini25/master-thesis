# Load required libraries ---
library(ggplot2)
library(dplyr)
library(animation)
library(sp)
library(spacetime)
library(gstat)
library(RColorBrewer)
library(tidyverse)
library(tidyr)
library(fields)
library(STRbook)

# --- Define Parameters ---
domain <- "surfacewater"
parameter <- "Q"
from_date <- as.Date("2020-01-01")
to_date <- as.Date("2023-12-31")

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

# Ensure all columns are of the correct type
all_q_data$coordinates.x <- as.numeric(all_q_data$coordinates.x)
all_q_data$coordinates.y <- as.numeric(all_q_data$coordinates.y)
all_q_data$data <- as.Date(all_q_data$data)  # Ensure correct column name
all_q_data <- na.omit(all_q_data)  # Remove NAs


# Create spatial scatter plot with fixed color scale range
ggplot(all_q_data, aes(x = coordinates.x, y = coordinates.y, color = Q)) +
  geom_point(size = 2) +  # Reduce point size
  scale_color_viridis_c(option = "C", trans = "sqrt", limits = c(0, 30)) + 
  theme_minimal() +
  labs(title = "River Flow",
       x = "Longitude",
       y = "Latitude",
       color = "Flow") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),  # Ensure all x labels appear
        axis.text.y = element_text(size = 7))  # Ensure all y labels appear




# TIME SERIES PLOT (pag. 26)
# Select 12 random locations
set.seed(1234)  # Set seed for reproducibility
selected_locations <- all_q_data %>%
  distinct(location_name) %>%
  sample_n(12) %>%
  pull(location_name)

# Filter data for selected locations
filtered_data <- all_q_data %>%
  filter(location_name %in% selected_locations)

# Create time-series plots
ggplot(filtered_data, aes(x = data, y = Q)) +
  geom_line() +
  facet_wrap(~ location_name, scales = "free_y") +
  theme_minimal(base_size = 8) +  # Reduce overall text size
  labs(
    title = "Time-Series Plots of River Flow",
    x = "Year",
    y = "Q (River Flow)"
  ) +
  theme(
    axis.title = element_text(size = 8),       # Smaller axis labels
    axis.text = element_text(size = 6),        # Smaller tick labels
    strip.text = element_text(size = 7),       # Smaller facet panel titles
    legend.text = element_text(size = 7),      # Smaller legend text
    legend.title = element_text(size = 8),     # Smaller legend title
    plot.caption = element_text(size = 6)      # Smaller caption text
  )


# Animation----
# Function to plot river flow (Q) for a given day
Q_t <- function(tau) {
  Q_sub <- filter(all_q_data, data == tau)  # Subset data for the specific date
  ggplot(Q_sub) +
    geom_point(aes(x = coordinates.x, y = coordinates.y, colour = Q), size = 4) +  # Plot flow data
    scale_color_gradient(name = "Q", low = "blue", high = "red") +  # Color scale
    theme_bw() +  # Black & white theme
    labs(title = paste("River Flow on", tau), x = "Longitude", y = "Latitude")
}

# Function to create animation
gen_anim <- function() {
  unique_dates <- unique(all_q_data$data)  # Get all unique dates
  for (t in unique_dates) {
    plot(Q_t(t))  # Generate plot for each date
  }
}

# Set animation options
ani.options(interval = 0.2)  # 0.2s interval between frames

# Save animation as an HTML file
saveHTML(gen_anim(),
         autoplay = FALSE, loop = FALSE, verbose = FALSE,
         outdir = ".", single.opts = "'controls': ['first', 'previous', 'play', 'next', 'last', 'loop', 'speed']",
         htmlfile = "river_flow_animation.html")  # Output filename

# Exploratory Analysis of Spatio-Temporal Data (pag 34) ----
# Empirical spatial mean 
spatial_mean_Q <- all_q_data %>%
  group_by(location_name, coordinates.x, coordinates.y) %>%
  summarise(Q_mean = mean(Q, na.rm = TRUE))  # Mean Q at each location

# Plot 1: Mean Q vs. Longitude
p1 <- ggplot(spatial_mean_Q, aes(x = coordinates.x, y = Q_mean)) +
  geom_point() +
  labs(title = "Mean River Flow vs. Longitude",
       x = "Longitude", y = "Mean River Flow") +
  theme_minimal()

# Plot 2: Mean Q vs. Latitude
p2 <- ggplot(spatial_mean_Q, aes(x = coordinates.y, y = Q_mean)) +
  geom_point() +
  labs(title = "Mean River Flow vs. Latitude",
       x = "Latitude", y = "Mean River Flow") +
  theme_minimal()

# Display plots
print(p1)
print(p2)

# Compute the empirical temporal mean: mean Q per date
temporal_mean_Q <- all_q_data %>%
  group_by(data) %>%  # Group by date
  summarise(Q_mean = mean(Q, na.rm = TRUE))  # Mean Q for each day

# Plot: Temporal Mean Q over Time
p3 <- ggplot(temporal_mean_Q, aes(x = data, y = Q_mean)) +
  geom_line(color = "blue", linewidth = 1) +  # Line plot of mean Q over time
  labs(title = "Mean River Flow Over Time",
       x = "Date", y = "Mean River Flow") +
  theme_minimal()

# Display plot
print(p3)

#Empirical covariances 
# Rename coordinates
all_q_data <- all_q_data %>%
  rename(lon = coordinates.x, lat = coordinates.y)

# Fit a linear model to remove trends (linear, not quadratic)
lm1 <- lm(Q ~ lon + data, data = all_q_data)  # Fit a linear model
all_q_data$residuals <- residuals(lm1)  # Store residuals

# Extract spatial locations of the stations
spat_df <- filter(all_q_data, data == "2020-01-01") %>%
  select(lon, lat) %>%
  arrange(lon, lat)

m <- nrow(spat_df)  # Number of stations

# Reshape data into space-wide format
X <- select(all_q_data, lon, lat, residuals, data) %>%
  spread(data, residuals) %>%
  select(-lon, -lat) %>%
  t()

# Compute empirical covariance matrices
Lag0_cov <- cov(X, use = 'complete.obs')
Lag1_cov <- cov(X[-1, ], X[-nrow(X), ], use = 'complete.obs')

# Split domain into longitude strips
spat_df$n <- 1:nrow(spat_df)
lim_lon <- range(spat_df$lon)
lon_strips <- seq(lim_lon[1], lim_lon[2], length = 5)

spat_df$lon_strip <- cut(spat_df$lon, lon_strips, labels = FALSE, include.lowest = TRUE)

# Plot empirical covariance matrices
plot_cov_strips(Lag0_cov, spat_df)  # Plot the lag-0 matrices
plot_cov_strips(Lag1_cov, spat_df)  # Plot the lag-1 matrices


# Spatio temporal Kringing (pag 172) non funxiona ----
# Create a SpatialPointsDataFrame
coordinates(all_q_data) <- ~coordinates.x + coordinates.y

# Define space and time lags
sp_lags <- 80   # Spatial binning (adjust based on dataset extent)
time_lags <- 30  # Maximum temporal lag

# Compute Spatio-Temporal Variogram
vv <- variogramST(
  Q ~ 1,            # Model formula
  data = all_q_data, # Data in spatio-temporal format
  tlags = 0:time_lags,  # Temporal lags (0 to 6 days)
  cutoff = 1000,    # Spatial cutoff distance
  width = sp_lags   # Spatial bin size
)

# Fit a Separable Variogram Model
sepVgm <- vgmST(
  stModel = "separable",
  space = vgm(10, "Exp", 400, nugget = 0.1),
  time = vgm(10, "Exp", 1, nugget = 0.1),
  sill = 20
)
sepVgm <- fit.StVariogram(vv, sepVgm)

# Fit an Alternative Metric Variogram Model
metricVgm <- vgmST(
  stModel = "metric",
  joint = vgm(100, "Exp", 400, nugget = 0.1),
  sill = 10,
  stAni = 100
)
metricVgm <- fit.StVariogram(vv, metricVgm)

# Compare Fits Using Mean Squared Error
metricMSE <- attr(metricVgm, "optim")$value
sepMSE <- attr(sepVgm, "optim")$value

# Plot the Fitted Variograms
plot(vv, list(sepVgm, metricVgm), main = "Empirical and Fitted Variograms for River Flow")

# Create Spatio-Temporal Prediction Grid
spat_pred_grid <- expand.grid(
  lon = seq(min(all_q_data$coordinates.x), max(all_q_data$coordinates.x), length = 20),
  lat = seq(min(all_q_data$coordinates.y), max(all_q_data$coordinates.y), length = 20)
) %>%
  SpatialPoints(proj4string = CRS(proj4string(all_q_data)))

# Temporal Prediction Grid (Daily Predictions)
temp_pred_grid <- seq(min(all_q_data$data), max(all_q_data$data), by = "1 day")

# Combine into STF Object
DE_pred <- STF(sp = spat_pred_grid, time = temp_pred_grid)

# Convert Data to STIDF Format
all_q_data_stidf <- as(all_q_data, "STIDF")
all_q_data_stidf <- subset(all_q_data_stidf, !is.na(all_q_data_stidf$Q))

# Perform Spatio-Temporal Kriging
pred_kriged <- krigeST(
  Q ~ 1, data = all_q_data_stidf,
  newdata = DE_pred, modelList = sepVgm,
  computeVar = TRUE
)

# Define Color Palette
color_pal <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(16))

# Plot Kriging Predictions
stplot(pred_kriged,
       main = "River Flow Predictions",
       layout = c(3, 2),
       col.regions = color_pal)

# Plot Kriging Standard Errors
pred_kriged$se <- sqrt(pred_kriged$var1.var)
stplot(pred_kriged[, "se"],
       main = "Prediction Standard Errors",
       layout = c(3, 2),
       col.regions = color_pal)

# Monthly  ----
# --- Define Parameters ---
from_date <- as.Date("2015-01-01")
to_date <- as.Date("2019-12-31")

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
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "m", from_date, to_date)
  
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

# Ensure all columns are of the correct type
all_q_data$coordinates.x <- as.numeric(all_q_data$coordinates.x)
all_q_data$coordinates.y <- as.numeric(all_q_data$coordinates.y)
all_q_data$data <- as.Date(all_q_data$data)  # Ensure correct column name
all_q_data <- na.omit(all_q_data)  # Remove NAs

# TIME SERIES PLOT (pag. 26)
# Select 12 random locations
set.seed(1234)  # Set seed for reproducibility
selected_locations <- all_q_data %>%
  distinct(location_name) %>%
  sample_n(12) %>%
  pull(location_name)

# Filter data for selected locations
filtered_data <- all_q_data %>%
  filter(location_name %in% selected_locations)

# Create time-series plots
ggplot(filtered_data, aes(x = data, y = Q)) +
  geom_line() +
  facet_wrap(~ location_name, scales = "free_y") +
  theme_minimal(base_size = 8) +  # Reduce overall text size
  labs(
    title = "Time-Series Plots of Monthly River Flow from 2015 to 2019",
    x = "Year",
    y = "Q (River Flow)"
  ) +
  theme(
    axis.title = element_text(size = 8),       # Smaller axis labels
    axis.text = element_text(size = 6),        # Smaller tick labels
    strip.text = element_text(size = 7),       # Smaller facet panel titles
    legend.text = element_text(size = 7),      # Smaller legend text
    legend.title = element_text(size = 8),     # Smaller legend title
    plot.caption = element_text(size = 6)  
  )



# Annual  ----
# --- Define Parameters ---
from_date <- as.Date("2010-01-01")
to_date <- as.Date("2024-12-31")

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
  q_data_response <- fetch_time_series_data(domain, location_code, parameter, "y", from_date, to_date)
  
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

# Ensure all columns are of the correct type
all_q_data$coordinates.x <- as.numeric(all_q_data$coordinates.x)
all_q_data$coordinates.y <- as.numeric(all_q_data$coordinates.y)
all_q_data$data <- as.Date(all_q_data$data)  # Ensure correct column name
all_q_data <- na.omit(all_q_data)  # Remove NAs

# TIME SERIES PLOT (pag. 26)
# Select 12 random locations
set.seed(1234)  # Set seed for reproducibility
selected_locations <- all_q_data %>%
  distinct(location_name) %>%
  sample_n(12) %>%
  pull(location_name)

# Filter data for selected locations
filtered_data <- all_q_data %>%
  filter(location_name %in% selected_locations)

# Create time-series plots
ggplot(filtered_data, aes(x = data, y = Q)) +
  geom_line() +
  facet_wrap(~ location_name, scales = "free_y") +
  theme_minimal(base_size = 8) +  # Reduce overall text size
  labs(
    title = "Time-Series Plots of Annual River Flow from 2010 to 2025",
    x = "Year",
    y = "Q (River Flow)"
  ) +
  theme(
    axis.title = element_text(size = 8),       # Smaller axis labels
    axis.text = element_text(size = 6),        # Smaller tick labels
    strip.text = element_text(size = 7),       # Smaller facet panel titles
    legend.text = element_text(size = 7),      # Smaller legend text
    legend.title = element_text(size = 8),     # Smaller legend title
    plot.caption = element_text(size = 6)  
  )

# HOVMOLLER DIAGRAM FOR Q ----
# Done using daily data 
# Aggregate data by location and date (if necessary)
hovmoller_data <- all_q_data %>%
  group_by(data, location_name) %>%  # Group by time and location
  summarise(Q = mean(Q, na.rm = TRUE), .groups = "drop")  # Take mean if multiple records exist

# Convert location to factor for correct ordering in plot
hovmoller_data <- hovmoller_data %>%
  mutate(location_name = factor(location_name, levels = unique(location_name)))

# Create the Hovmöller plot
ggplot(hovmoller_data, aes(x = data, y = location_name, fill = Q)) +
  geom_tile() +  # Creates a heatmap-like visualization
  scale_fill_viridis_c(option = "C", trans = "sqrt", limits = c(0, 200)) +  # Fix here: scale_fill instead of scale_color
  theme_minimal() +
  labs(title = "Hovmöller Plot of River Flow (Q)",
       x = "Date",
       y = "Location",
       fill = "Q (River Flow)") +
  theme(axis.text.y = element_text(size = 8))  # Adjust location labels




# PRECIPITATIONS ----
# --- Define Parameters ---
domain_meteo <- "meteo"
parameter_prec <- "Prec"
from_date <- as.Date("2020-01-01")
to_date <- as.Date("2023-12-31")

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

# Ensure all columns are of the correct type
all_prec_data$coordinates.x <- as.numeric(all_prec_data$coordinates.x)
all_prec_data$coordinates.y <- as.numeric(all_prec_data$coordinates.y)
all_prec_data$data <- as.Date(all_prec_data$data)  # Ensure correct column name
all_prec_data <- na.omit(all_prec_data)  # Remove NAs

# Create spatial scatter plot with fixed color scale range
ggplot(all_prec_data, aes(x = coordinates.x, y = coordinates.y, color = prec)) +
  geom_point(size = 2) +  # Reduce point size
  scale_color_viridis_c(option = "C", trans = "sqrt", limits = c(0, 30)) + 
  theme_minimal() +
  labs(title = "Precipitation amount",
       x = "Longitude",
       y = "Latitude",
       color = "Precipitation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),  # Ensure all x labels appear
        axis.text.y = element_text(size = 7))  # Ensure all y labels appear



# Average precipitation plot 
# Ensure 'data' column is in Date format
all_prec_data$data <- as.Date(all_prec_data$data)

# Compute the mean river flow across all locations for each date
mean_prec_time_series <- all_prec_data %>%
  group_by(data) %>%
  summarise(mean_prec = mean(prec, na.rm = TRUE))

# Compute the mean precipitation across all locations for each date
mean_prec_time_series <- all_prec_data %>%
  group_by(data) %>%
  summarise(mean_prec = mean(prec, na.rm = TRUE)) %>%
  filter(mean_prec <= 200)  # Remove rows where mean_prec > 200

# Plot the mean river flow time series
ggplot(mean_prec_time_series, aes(x = data, y = mean_prec)) +
  geom_line(color = "blue", linewidth = 0.3) +
  theme_bw() +
  xlab("Year") +
  ylab("Mean Precipitation Amount") +
  ggtitle("Time Series of Mean Precipitation Across All Locations") +
  theme(axis.text.x = element_text(hjust = 0.5))  # Rotate x-axis labels for readability

# HOVMOLLER DIAGRAM for Precipitation ----
from_date <- as.Date("2015-01-01")
to_date <- as.Date("2024-12-31")

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
  prec_data_response <- fetch_time_series_data(domain_meteo, location_code, parameter_prec, "m", from_date, to_date)
  
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

# Ensure all columns are of the correct type
all_prec_data$coordinates.x <- as.numeric(all_prec_data$coordinates.x)
all_prec_data$coordinates.y <- as.numeric(all_prec_data$coordinates.y)
all_prec_data$data <- as.Date(all_prec_data$data)  # Ensure correct column name
all_prec_data <- na.omit(all_prec_data)  # Remove NAs


# Aggregate data by location and date (if necessary)
hovmoller_data_prec <- all_prec_data %>%
  group_by(data, location_name) %>%  # Group by time and location
  summarise(mean_prec = mean(prec, na.rm = TRUE), .groups = "drop")  # Take mean if multiple records exist

# Convert location to factor for correct ordering in plot
hovmoller_data_prec <- hovmoller_data_prec %>%
  mutate(location_name = factor(location_name, levels = unique(location_name)))

# Create the Hovmöller plot
ggplot(hovmoller_data_prec, aes(x = data, y = location_name, fill = mean_prec)) +
  geom_tile() +  # Creates a heatmap-like visualization
  scale_fill_viridis_c(option = "C", trans = "sqrt", limits = c(0, 000)) +  # Fix here: scale_fill instead of scale_color
  theme_minimal() +
  labs(title = "Hovmöller Plot Precipitation",
       x = "Date",
       y = "Location",
       fill = "Precipitation amount") +
  theme(axis.text.y = element_text(size = 8))  # Adjust location labels



# EMPIRICAL ORTHOGONAL FUNCTION (pag 88)----
# I USE THE FIRST all_q_data
from_date <- as.Date("2020-01-01")
to_date <- as.Date("2023-12-31")

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

# Ensure all columns are of the correct type
all_q_data$coordinates.x <- as.numeric(all_q_data$coordinates.x)
all_q_data$coordinates.y <- as.numeric(all_q_data$coordinates.y)
all_q_data$data <- as.Date(all_q_data$data)  # Ensure correct column name

# Convert to time-wide format
q_wide <- all_q_data %>%
  pivot_wider(names_from = data, values_from = Q)

# Remove rows (locations) with more than 50% missing values
threshold <- ncol(q_wide) / 2  # Compute 50% threshold
q_wide_filtered <- q_wide[rowSums(is.na(q_wide)) <= threshold, ]

# Apply forward and backward filling before interpolation
q_wide_imputed <- q_wide_filtered %>%
  mutate(across(-c(location_name, coordinates.x, coordinates.y), 
                ~ na.locf(na.locf(.x, na.rm = FALSE), fromLast = TRUE), .names = "{.col}")) %>%
  mutate(across(-c(location_name, coordinates.x, coordinates.y), 
                ~ na.approx(.x, na.rm = FALSE), .names = "{.col}"))

# Check for remaining NAs
sum(is.na(q_wide_imputed))  # Should be 0 if all missing values are filled

# 2️⃣ Convert to space-wide format (transpose matrix)
Z <- as.matrix(q_wide_imputed[, -c(1:3)])  # Extract only Q values
Z <- t(Z)  # Transpose: Time -> Rows, Space -> Columns
dim(Z) # 1461   37

# 3️⃣ Compute mean and standardize data
nT <- ncol(Z)  # Number of spatial locations

# Compute row-wise mean (each time step)
spat_mean <- rowMeans(Z, na.rm = TRUE)

# Subtract mean from each row
Zspat_detrend <- sweep(Z, 1, spat_mean, "-")

Zt <- 1 / sqrt(nT - 1) * Zspat_detrend  # Standardize

# 4️⃣ Compute Singular Value Decomposition (SVD)
E <- svd(Zt)

# 5️⃣ Extract EOFs
V <- E$v  # Spatial EOFs
colnames(V) <- paste0("EOF", 1:ncol(V))  # Rename columns

# Append coordinates and location names
EOFs <- cbind(q_wide_imputed[, c("location_name", "coordinates.x", "coordinates.y")], V)

# 6️⃣ Convert Principal Components (U matrix) to long format
TS <- data.frame(E$u) %>%
  mutate(t = 1:nrow(E$u)) %>%
  pivot_longer(-t, names_to = "EOF", values_to = "PC")

# 7️⃣ Normalize the time series
TS$nPC <- TS$PC * sqrt(nT - 1)

# 8️⃣ Plot the first EOF spatial pattern
# WORKS 
ggplot(EOFs) +
  geom_point(aes(x = coordinates.x, y = coordinates.y, color = EOF1), size = 4) +  
  geom_text(aes(x = coordinates.x, y = coordinates.y, label = location_name), 
            hjust = -0.1, vjust = 0.2, size = 2, check_overlap = TRUE) +  # Adjust position & avoid clutter
  scale_color_viridis_c(name = "EOF1") +
  theme_bw() +
  xlab("Longitude") + 
  ylab("Latitude") +
  ggtitle("Spatial Pattern of EOF1")


# Second PC
ggplot(EOFs) +
  geom_point(aes(x = coordinates.x, y = coordinates.y, color = EOF2), size = 4) +  
  geom_text(aes(x = coordinates.x, y = coordinates.y, label = location_name), 
            hjust = -0.1, vjust = 0.2, size = 2, check_overlap = TRUE) +  
  scale_color_viridis_c(name = "EOF2") +
  theme_bw() +
  xlab("Longitude") + 
  ylab("Latitude") +
  ggtitle("Spatial Pattern of EOF2")


# Third PC
ggplot(EOFs) +
  geom_point(aes(x = coordinates.x, y = coordinates.y, color = EOF3), size = 4) +  
  geom_text(aes(x = coordinates.x, y = coordinates.y, label = location_name), 
            hjust = -0.1, vjust = 0.1, size = 3, check_overlap = TRUE) +  
  scale_color_viridis_c(name = "EOF3") +
  theme_bw() +
  xlab("Longitude") + 
  ylab("Latitude") +
  ggtitle("Spatial Pattern of EOF3")


# Convert principal component time series to long format
TS <- data.frame(E$u) %>%
  mutate(t = 1:nrow(E$u)) %>%
  pivot_longer(-t, names_to = "EOF", values_to = "PC")

# Normalize time series
TS$nPC <- TS$PC * sqrt(nT - 1)

# Convert time step `t` into actual date sequence
start_date <- as.Date("2020-01-01")  
TS$Date <- start_date + (TS$t - 1)  # Convert time steps to actual dates

# Plot PC1 with years on x-axis
ggplot(TS %>% filter(EOF == "X1")) +  
  geom_line(aes(x = Date, y = nPC), color = "blue") +
  theme_bw() +
  xlab("Year") +
  ylab("Normalized PC1") +
  ggtitle("Principal Component Time Series - PC1") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")  # Set x-axis to yearly ticks

# Plot PC2 (Second Principal Component)
ggplot(TS %>% filter(EOF == "X2")) +  
  geom_line(aes(x = Date, y = nPC), color = "red") +
  theme_bw() +
  xlab("Year") +
  ylab("Normalized PC2") +
  ggtitle("Principal Component Time Series - PC2") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")  # Set x-axis to yearly ticks


# Plot PC3 (Third Principal Component)
ggplot(TS %>% filter(EOF == "X3")) +  
  geom_line(aes(x = Date, y = nPC), color = "green") +
  theme_bw() +
  xlab("Year") +
  ylab("Normalized PC3") +
  ggtitle("Principal Component Time Series - PC3") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")  # Set x-axis to yearly ticks

# TS OF MEAN Q
all_q_data$data <- as.Date(all_q_data$data)

# Compute the mean river flow across all locations for each date
mean_q_time_series <- all_q_data %>%
  group_by(data) %>%
  summarise(mean_Q = mean(Q, na.rm = TRUE))

# Plot the mean river flow time series
ggplot(mean_q_time_series, aes(x = data, y = mean_Q)) +
  geom_line(color = "blue", linewidth = 1) +
  theme_bw() +
  xlab("Year") +
  ylab("Mean River Flow (Q)") +
  ggtitle("Time Series of Mean River Flow Across All Locations") +
  theme(axis.text.x = element_text(hjust = 0.5))

# Explained variance barplot  
# Compute variance explained by each EOF (PC)
variance_explained <- E$d^2 / sum(E$d^2) * 100  # Convert singular values to percentage variance explained

# Create a data frame for plotting (only first 10 PCs)
variance_df <- data.frame(
  PC = factor(1:10, levels = 1:10),  # Select only the first 10 PCs
  Variance = variance_explained[1:10]  # First 10 variance values
)

# Plot bar chart of variance explained (first 10 PCs)
ggplot(variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_bw() +
  xlab("Principal Component") +
  ylab("Variance Explained (%)") +
  ggtitle("Variance Explained by First 10 Principal Components") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Create a table with the exact variance explained for the first 10 PCs
variance_table <- data.frame(
  PC = 1:10,
  Variance_Explained = round(variance_explained[1:10], 2)  # Round to 2 decimal places
)

# Print the table
print(variance_table)