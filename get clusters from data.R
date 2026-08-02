#' Using USGS 15 min turbidity data
#' Analyse turbidity patterns and cluster stations based on diurnal
#' turbidity pattern
#' Authors: Galina Shinkareva, Wondwosen Seyoum


# ---- 1. Load packages ----

library(tidyverse)
library(moments)
library(gridExtra)

# ---- 2. Load and Clean Data ----
## ---- Set up path ----
current_wd <- getwd()
data_dir <- "input_data"
fig_dir <- "figures"
out_dir <- "output_data"

# Create directories if they do not exist
if (!dir.exists(fig_dir)) dir.create(fig_dir)
if (!dir.exists(out_dir)) dir.create(out_dir)

# 15 min turbidity data from USGS gauging station is used
# file_name = "USGS_15min_turbidity.rds"
file_name = "LA_turb_disch.rds"

## ---- Read data ----
data <-readRDS(file = file.path(data_dir, file_name)) %>% ungroup() # change the file name and location accordingly 

clean_data <- data %>%
  select(dateTime,
         X_63680_00000, site_no, state,
         dec_lat_va, dec_long_va)

rm(data) # remove initial data frame to free memory

# ---- 3. Get hourly data for each site ----
# Aggregate data from 15 min to hourly format

hourly_data <- clean_data %>%
  mutate(turb = if_else(X_63680_00000 < 0, 0, X_63680_00000)) %>% # Replace negative values for zeros
  mutate(hour = hour(dateTime)) %>%
  summarise(mean_turb = mean(turb, na.rm = TRUE), .by = c(site_no, hour)) %>%
  drop_na(mean_turb)

# Check for NA values
checkNA <- hourly_data %>%
  filter(is.na(mean_turb))

# ---- 4. Transform data ----

# transform into wide format
data_wide <- hourly_data %>%
  # filter(mean_turb >= 0) %>%
  pivot_wider(names_from = hour, values_from = mean_turb, names_prefix = "time_")

# ---- 5. Calculate stats for data ----

station_ids <- unique(data_wide$site_no)
feature_matrix <- data.frame(station_id = station_ids)

# Loop to go through each USGS station available
for (station in station_ids) {
  
  # Select turbidity columns
  station_data <- data_wide %>%
    filter(site_no == station) %>%
    select(starts_with("time_")) %>%
    unlist() %>%
    as.numeric()
  
  # Warn if there are any NAs
  if(any(is.na(station_data))) {
    warning("NA values found for station ", station)
  }
  
  # Replace NA with 0
  station_data[is.na(station_data)] <- 0
  
  # Skip if all values are zero
  if(all(station_data == 0)) {
    message("Skipping station ", station, " (all zeros)")
    next
  }
  
  # --- FFT Features ---
  # Fast Fourier Transform (FFT) to calculate frequencies
  fft_result <- fft(station_data)
  amplitude <- abs(fft_result)
  
  # Safe indexing
  dominant_freq <- if(length(amplitude) > 2) which.max(amplitude[2:min(12, length(amplitude))]) else NA # Exclude DC component and frequencies above 12
  dominant_amp <- if(!is.na(dominant_freq)) amplitude[dominant_freq + 1] else NA
  phase_shift <- if(!is.na(dominant_freq)) Arg(fft_result)[dominant_freq + 1] else NA
  
  feature_matrix[feature_matrix$station_id == station, c("dominant_freq", "dominant_amplitude", "phase_shift")] <-
    c(dominant_freq, dominant_amp, phase_shift)
  
  # --- Polynomial Fit ---
  time_values <- 1:24
  poly_fit <- tryCatch(lm(station_data ~ poly(time_values, 2)), error = function(e) NULL)
  if(!is.null(poly_fit)) {
    coefs <- coef(poly_fit)
    feature_matrix[feature_matrix$station_id == station, "poly_intercept"] <- coefs[1]
    feature_matrix[feature_matrix$station_id == station, "poly_linear"] <- coefs[2]
    feature_matrix[feature_matrix$station_id == station, "poly_quadratic"] <- ifelse(length(coefs) >= 3, coefs[3], NA)
    feature_matrix[feature_matrix$station_id == station, "poly_r_squared"] <- summary(poly_fit)$r.squared
  } else {
    feature_matrix[feature_matrix$station_id == station, c("poly_intercept","poly_linear","poly_quadratic","poly_r_squared")] <- NA
  }
  
  # --- Statistical Features ---
  feature_matrix[feature_matrix$station_id == station, "mean"] <- mean(station_data, na.rm = TRUE)
  feature_matrix[feature_matrix$station_id == station, "sd"] <- sd(station_data, na.rm = TRUE)
  feature_matrix[feature_matrix$station_id == station, "min"] <- min(station_data, na.rm = TRUE)
  feature_matrix[feature_matrix$station_id == station, "max"] <- max(station_data, na.rm = TRUE)
  feature_matrix[feature_matrix$station_id == station, "skewness"] <- ifelse(all(station_data==0), NA, skewness(station_data))
  feature_matrix[feature_matrix$station_id == station, "kurtosis"] <- ifelse(all(station_data==0), NA, kurtosis(station_data))
  
  # Time of max and min
  feature_matrix[feature_matrix$station_id == station, "max_time"] <- which.max(station_data)
  feature_matrix[feature_matrix$station_id == station, "min_time"] <- which.min(station_data)
}

# ---- 5. Scale the features ----
#and replace NA with column mean
scaled_features <- feature_matrix[, -1] %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  scale() # Exclude station_id

# ---- 7. Determine optimal number of clusters (using elbow method) ----
set.seed(12345)
total_deg_fr <- nrow(scaled_features) - 1 # degrees of freedom associated with the total variance
wss <- total_deg_fr * sum(apply(scaled_features, 2, var))

# Maximum number of clusters requested
max_clusters_requested <- 13

# Automatically use the smaller of:
# 1) requested maximum, or
# 2) number of observations -1 which is equal to degrees of freedom
max_clusters <- min(max_clusters_requested, total_deg_fr)

for (i in 2:max_clusters) {
  wss[i] <- sum(kmeans(scaled_features, centers = i)$withinss)
}

# Open device
plot_name <- "ElbowPlot.jpeg" 

jpeg(filename = file.path(fig_dir, plot_name), 
      width = 20, height = 15, units = "cm", res = 300)

plot(1:max_clusters, wss, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")

# Close device and write file
dev.off()

# ---- 8. Apply K-means clustering (based on elbow method) ----

k <- 4 # Choose the optimal number of clusters from the elbow plot
set.seed(12345)
kmeans_result <- kmeans(scaled_features, centers = k)

## ---- Add cluster assignments to the feature matrix ---- 

##This combines all variables used for clustering, stations id, and cluster assignments
feature_matrix$cluster <- kmeans_result$cluster
#The feature_matrix dataframe has both stations id and cluster group
#can be used for model building

## ---- Visualize clusters (example: plot mean curves for each cluster) ---- 
cluster_means <- hourly_data %>%
  left_join(feature_matrix[, c("station_id", "cluster")], by = c("site_no" = "station_id")) %>%
  summarise(mean_variable = mean(mean_turb), .by = c(cluster, hour)) 

cluster_data <- hourly_data %>%
  left_join(feature_matrix[, c("station_id", "cluster")], by = c("site_no" = "station_id")) 

cluster_data %>%
  summarise(site_count = n_distinct(site_no), .by = "cluster")
#cluster_data dataframe also has the station id and cluster assignment

plot <- cluster_means %>%

  # patterns can be renamed here:
  # mutate(pattern_new = case_when(
  #   cluster == 1 ~ "Day max",
  #   cluster == 2 ~ "Morninig max",
  #   cluster == 3 ~ "Night max",
  #   cluster == 4 ~ "Evening max"
  #   TRUE ~ "No"
  # )) %>%
  ggplot(aes(x = hour, y = mean_variable, color = factor(cluster))) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ cluster, scales = 'free_y') + #scales = 'free' for both axes labelled
  labs(title = "Mean Turbidity by Cluster", x = "Time", y = "Variable") +
  theme_minimal() +
  guides(color="none")

plot

# Create table from cluster counts
cluster_table <- as.data.frame(table(feature_matrix$cluster))
tbl <- tableGrob(cluster_table)

# Add table to plot
p <- grid.arrange(plot, tbl, ncol = 2, widths = c(3, 1))  # 3:1 width ratio
p

# ---- 9. Save results ---- 
## ---- Save plot ---- 
plot_name <- "mean_turb_by_clusters.jpeg"

ggsave(filename = file.path(fig_dir, plot_name),
       plot = p, width = 40, height = 20, units = "cm",
       dpi = 300)


cluster_table

Stations_cluster <- feature_matrix[,c(1,17)]

## ---- Save new data frame ---- 
file_name = "turb_with_clusters.csv"
write_csv(Stations_cluster, file.path(out_dir, file_name))
