library(httr)
library(jsonlite)
library(sf)
library(leaflet)
library(ggrepel)

# RIVER SENSORS MAP ----
# --- Define Parameters ---
domain <- "surfacewater"
parameter <- "Q"
from_date <- as.Date("2020-01-01")
to_date <- as.Date("2023-12-31")

# Parameters for precipitations 
domain_meteo <- "meteo"
parameter_prec <- "Prec"

# Function to fetch water data
fetch_water_data <- function(domain, from_date, to_date, parameter = NULL) {
  base_url <- "http://www.oasi.ti.ch/web/rest/locations"
  params <- list(
    domain = domain,
    from = from_date,
    to = to_date
  )
  
  if (!is.null(parameter)) {
    params$parameter <- parameter
  }
  
  response <- GET(base_url, query = params)
  
  if (response$status_code == 200) {
    response_text <- content(response, as = "text")
    
    if (nzchar(response_text)) {  # Check if response is not empty
      data <- tryCatch({
        fromJSON(response_text, flatten = TRUE)
      }, error = function(e) {
        print(paste("Error parsing JSON:", e))
        return(NULL)
      })
      return(data)
    } else {
      print("Empty response received from API.")
      return(NULL)
    }
  } else {
    print(paste("Failed to get data. Status code:", response$status_code))
    return(NULL)
  }
}

# Fetch data
data <- fetch_water_data(domain = "surfacewater", from_date = "2022-07-19", to_date = "2024-07-18", parameter = "Q")

if (!is.null(data)) {
  # Check if coordinates exist
  if (!all(c("coordinates.x", "coordinates.y") %in% names(data))) {
    print("Error: Coordinate fields are missing.")
  } else {
    # Convert coordinates to numeric
    data$coordinates.x <- as.numeric(data$coordinates.x)
    data$coordinates.y <- as.numeric(data$coordinates.y)
    
    # Convert to sf object with correct Swiss CRS (EPSG:2056)
    data_sf <- st_as_sf(data, coords = c("coordinates.x", "coordinates.y"), crs = 2056)
    
    # Transform to WGS84 (EPSG:4326)
    data_wgs84 <- st_transform(data_sf, crs = 4326)
    
    # Extract transformed coordinates
    coords <- st_coordinates(data_wgs84)
    data_wgs84$longitude <- coords[,1]
    data_wgs84$latitude <- coords[,2]
    
    # Simple leaflet map with blue markers
    leaflet(data = data_wgs84) %>%
      addTiles() %>%
      addCircleMarkers(
        ~longitude, ~latitude,
        popup = ~name,
        label = ~name,
        color = 'blue',
        radius = 5
      )
  }
}


# METEO STATION MAP ----
# Function to fetch precipitation data
fetch_precipitation_data <- function(domain, from_date, to_date, parameter = NULL) {
  base_url <- "http://www.oasi.ti.ch/web/rest/locations"
  params <- list(
    domain = domain,
    from = from_date,
    to = to_date
  )
  
  if (!is.null(parameter)) {
    params$parameter <- parameter
  }
  
  response <- GET(base_url, query = params)
  
  if (response$status_code == 200) {
    response_text <- content(response, as = "text")
    
    if (nzchar(response_text)) {  # Check if response is not empty
      data <- tryCatch({
        fromJSON(response_text, flatten = TRUE)
      }, error = function(e) {
        print(paste("Error parsing JSON:", e))
        return(NULL)
      })
      return(data)
    } else {
      print("Empty response received from API.")
      return(NULL)
    }
  } else {
    print(paste("Failed to get data. Status code:", response$status_code))
    return(NULL)
  }
}

# Fetch precipitation data
precipitation_data <- fetch_precipitation_data(
  domain = "meteo",
  from_date = "2022-07-19",
  to_date = "2024-07-18",
  parameter = "Prec"
)

if (!is.null(precipitation_data)) {
  # Check if coordinates exist
  if (!all(c("coordinates.x", "coordinates.y") %in% names(precipitation_data))) {
    print("Error: Coordinate fields are missing.")
  } else {
    # Convert coordinates to numeric
    precipitation_data$coordinates.x <- as.numeric(precipitation_data$coordinates.x)
    precipitation_data$coordinates.y <- as.numeric(precipitation_data$coordinates.y)
    
    # Convert to sf object with correct Swiss CRS (EPSG:2056)
    prec_sf <- st_as_sf(precipitation_data, coords = c("coordinates.x", "coordinates.y"), crs = 2056)
    
    # Transform to WGS84 (EPSG:4326)
    prec_wgs84 <- st_transform(prec_sf, crs = 4326)
    
    # Extract transformed coordinates
    coords <- st_coordinates(prec_wgs84)
    prec_wgs84$longitude <- coords[,1]
    prec_wgs84$latitude <- coords[,2]
    
    # Simple leaflet map with blue markers for precipitation locations
    leaflet(data = prec_wgs84) %>%
      addTiles() %>%
      addCircleMarkers(
        ~longitude, ~latitude,
        popup = ~name,
        label = ~name,
        color = 'blue',
        radius = 5
      )
  }
}



# MAP CHAPTER 2 GAM REGRESSION ----
# Significant precipitation locations based on GAM results

significant_locations <- c("Arcisate", "Arosio", "Biasca", "Biasca_Val_Pontirone", 
                      "Bioggio", "Cabbio", "Camignolo", "Carena", "Cavargna", 
                      "Cavergno", "Chironico", "Chironico_Cala", "Colla", "Como",
                      "Diga_Luzzone", "Fusio", "Giubiasco", "Grancia", 
                      "Lago_Truzzo_S_Giacomo", "Lavertezzo__Aquino", "Locarno", 
                      "Luino", "Maggia", "Mendrisio", "Novaggio", "Novazzano", 
                      "Olivone", "Piora_CBA", "Piotta", "Preonzo__Alpe_Roscioro",
                      "Robiei", "Sonogno", "Sonvico", "Stabio")


# Special location to highlight in blue
special_location <- "Bolletta - Porto Ceresio"

# Function to fetch precipitation data
fetch_precipitation_data <- function(domain, from_date, to_date, parameter = NULL) {
  base_url <- "http://www.oasi.ti.ch/web/rest/locations"
  params <- list(
    domain = domain,
    from = from_date,
    to = to_date
  )
  
  if (!is.null(parameter)) {
    params$parameter <- parameter
  }
  
  response <- GET(base_url, query = params)
  
  if (response$status_code == 200) {
    response_text <- content(response, as = "text")
    
    if (nzchar(response_text)) {  # Check if response is not empty
      data <- tryCatch({
        fromJSON(response_text, flatten = TRUE)
      }, error = function(e) {
        print(paste("Error parsing JSON:", e))
        return(NULL)
      })
      return(data)
    } else {
      print("Empty response received from API.")
      return(NULL)
    }
  } else {
    print(paste("Failed to get data. Status code:", response$status_code))
    return(NULL)
  }
}

# Function to fetch water data
fetch_water_data <- function(domain, from_date, to_date, parameter = NULL) {
  base_url <- "http://www.oasi.ti.ch/web/rest/locations"
  params <- list(
    domain = domain,
    from = from_date,
    to = to_date
  )
  
  if (!is.null(parameter)) {
    params$parameter <- parameter
  }
  
  response <- GET(base_url, query = params)
  
  if (response$status_code == 200) {
    response_text <- content(response, as = "text")
    
    if (nzchar(response_text)) {  # Check if response is not empty
      data <- tryCatch({
        fromJSON(response_text, flatten = TRUE)
      }, error = function(e) {
        print(paste("Error parsing JSON:", e))
        return(NULL)
      })
      return(data)
    } else {
      print("Empty response received from API.")
      return(NULL)
    }
  } else {
    print(paste("Failed to get data. Status code:", response$status_code))
    return(NULL)
  }
}

# Fetch precipitation data
precipitation_data <- fetch_precipitation_data(
  domain = domain_meteo,
  from_date = from_date,
  to_date = to_date,
  parameter = parameter_prec
)

# Fetch water data for Bolletta - Porto Ceresio
data_water <- fetch_water_data(domain = "surfacewater", from_date = from_date, to_date = to_date, parameter = "Q")

if (!is.null(precipitation_data) && !is.null(data_water)) {
  # Check if coordinates exist
  if (!all(c("coordinates.x", "coordinates.y") %in% names(precipitation_data))) {
    print("Error: Coordinate fields are missing in precipitation data.")
  } else {
    # Convert coordinates to numeric for precipitation data
    precipitation_data$coordinates.x <- as.numeric(precipitation_data$coordinates.x)
    precipitation_data$coordinates.y <- as.numeric(precipitation_data$coordinates.y)
    
    # Convert to sf object with correct Swiss CRS (EPSG:2056)
    prec_sf <- st_as_sf(precipitation_data, coords = c("coordinates.x", "coordinates.y"), crs = 2056)
    
    # Transform to WGS84 (EPSG:4326)
    prec_wgs84 <- st_transform(prec_sf, crs = 4326)
    
    # Extract transformed coordinates
    coords <- st_coordinates(prec_wgs84)
    prec_wgs84$longitude <- coords[,1]
    prec_wgs84$latitude <- coords[,2]
    
    # Extract Bolletta - Porto Ceresio coordinates from water data
    bolletta_location <- data_water[data_water$name == special_location, ]
    
    if (nrow(bolletta_location) > 0) {
      bolletta_location$coordinates.x <- as.numeric(bolletta_location$coordinates.x)
      bolletta_location$coordinates.y <- as.numeric(bolletta_location$coordinates.y)
      bolletta_sf <- st_as_sf(bolletta_location, coords = c("coordinates.x", "coordinates.y"), crs = 2056)
      bolletta_wgs84 <- st_transform(bolletta_sf, crs = 4326)
      bolletta_coords <- st_coordinates(bolletta_wgs84)
      bolletta_wgs84$longitude <- bolletta_coords[,1]
      bolletta_wgs84$latitude <- bolletta_coords[,2]
      prec_wgs84 <- rbind(prec_wgs84, bolletta_wgs84)
    }
    
    # Define colors and sizes based on significance
    prec_wgs84$color <- ifelse(prec_wgs84$name == special_location, "blue", 
                               ifelse(prec_wgs84$name %in% significant_locations, "red", "black"))
    prec_wgs84$size <- ifelse(prec_wgs84$name == special_location, 10, 
                              ifelse(prec_wgs84$name %in% significant_locations, 10, 3))
    
    # Create a leaflet map
    leaflet(data = prec_wgs84) %>%
      addTiles() %>%
      addCircleMarkers(
        ~longitude, ~latitude,
        popup = ~name,
        label = ~name,
        color = ~color,
        radius = ~size,
        stroke = FALSE,
        fillOpacity = 0.8
      )
  }
}
