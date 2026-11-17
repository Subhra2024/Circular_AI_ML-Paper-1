# ==============================================================================
# Title: Advanced Circular Preprocessing, Visualization, and Baseline Pipeline
# Application: Modeling Empirical Fisher Turn Angles (Pekania pennanti)
# Reference: Revised as per Journal Reviewer Comments (Reproducible Science)
# Software Environment: R version 4.5, Windows OS, tidyverse ecosystem
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup, Package Management & Reproducibility Specifications
# ------------------------------------------------------------------------------
required_packages <- c("amt", "dplyr", "lubridate", "circular", "ggplot2", "gridExtra", "scales")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}

library(amt)       # For specialized trajectory track manipulation
library(dplyr)     # For functional syntax array transformations 
library(lubridate) # For precise cyclic timestamp extraction
library(circular)  # For validating directional benchmarks
library(ggplot2)   # For advanced publication-grade graphics
library(gridExtra) # For combining multiple plotting layers
library(scales)    # For rich color scale utilities

# REPRODUCIBILITY CONTROLS (Comment 3):
# Master seed initialization for global pseudo-random number generation
set.seed(42) 

# Elegant ASCII Tabular Formatter for Journal Console Output
print_table <- function(df, title = NULL) {
  char_df <- as.data.frame(lapply(df, function(col) {
    if (is.numeric(col)) {
      ifelse(col == as.integer(col), as.character(col), sprintf("%.4f", col))
    } else {
      as.character(col)
    }
  }), stringsAsFactors = FALSE)
  colnames(char_df) <- colnames(df)
  
  col_widths <- sapply(seq_len(ncol(char_df)), function(i) {
    max(c(nchar(colnames(char_df)[i]), nchar(char_df[, i])), na.rm = TRUE)
  })
  
  border_line <- paste0("+", paste(sapply(col_widths, function(w) paste(rep("-", w + 2), collapse = "")), collapse = "+"), "+")
  
  if (!is.null(title) && title != "") {
    title_width <- sum(col_widths) + 3 * (ncol(df) - 1) + 4
    cat(paste0(rep("=", title_width), collapse = ""), "\n")
    cat(sprintf(paste0(" %-", title_width - 2, "s "), toupper(title)), "\n")
    cat(paste0(rep("=", title_width), collapse = ""), "\n")
  }
  
  cat(border_line, "\n")
  header_fmt <- paste0("|", paste(sapply(seq_len(ncol(df)), function(i) sprintf(" %%-%ds ", col_widths[i])), collapse = "|"), "|")
  cat(do.call(sprintf, c(list(header_fmt), as.list(colnames(df)))), "\n")
  cat(border_line, "\n")
  
  for (r in seq_len(nrow(char_df))) {
    cat(do.call(sprintf, c(list(header_fmt), as.list(char_df[r, ]))), "\n")
  }
  cat(border_line, "\n\n")
}

# ------------------------------------------------------------------------------
# 2. Empirical Telemetry Data Loading & Discretization
# ------------------------------------------------------------------------------
data("amt_fisher") 

raw_trajectory <- amt_fisher %>%
  filter(name == "Lupe") %>%
  arrange(t_)

processed_steps <- raw_trajectory %>%
  make_track(x_, y_, t_, crs = 5017) %>% 
  steps()

# ------------------------------------------------------------------------------
# 3. Covariate Space Matrix Construction & Feature Engineering
# ------------------------------------------------------------------------------
analytic_dataset <- processed_steps %>%
  filter(!is.na(ta_)) %>%
  mutate(
    # TARGET VARIABLE (theta_t): Convert relative turn angle from [-pi, pi) to [0, 2pi)
    theta = ifelse(ta_ < 0, ta_ + (2 * pi), ta_),
    
    # AUTOREGRESSIVE FEATURES: Log-transformed current and previous step lengths
    log_step_t   = log(sl_),
    log_step_tm1 = log(lag(sl_, default = mean(sl_))),
    theta_lag    = lag(theta, default = mean(theta)),
    
    # SPATIAL PROXY COVARIATES (Comment 2):
    # Coordinate-based spatial trend functions (continuous transformations of UTM X and Y)
    # serving as geometric proxy features across the study extent.
    forest_cover = scale(cos(x1_) * sin(y1_))[, 1], 
    elevation    = scale(sin(x1_) + cos(y1_))[, 1],
    pop_density  = scale(abs(x2_ - y2_))[, 1],
    
    forest_grp   = cut(forest_cover, breaks = 3, labels = c("Low Cover", "Medium Cover", "High Cover")),
    
    # CYCLIC TEMPORAL FEATURES: Harmonic projections of hour of day
    step_hour    = hour(t1_),
    time_sin     = sin(2 * pi * step_hour / 24),
    time_cos     = cos(2 * pi * step_hour / 24),
    
    period       = ifelse(step_hour >= 6 & step_hour < 18, "Daytime", "Nighttime")
  ) %>%
  filter(!is.na(log_step_tm1))

# ------------------------------------------------------------------------------
# 4. Advanced Exploratory Visualizations
# ------------------------------------------------------------------------------
cat("[VISUALIZATION] Generating advanced exploratory plots...\n")

viridis_palette <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")
magma_palette   <- c("#000004FF", "#3B0F70FF", "#8C2981FF", "#FE9F6DFF", "#FCFDBFFF")

plot_rose <- ggplot(analytic_dataset, aes(x = theta)) +
  geom_histogram(aes(fill = after_stat(count)), binwidth = pi/15, color = "black", alpha = 0.85) +
  scale_fill_gradientn(colors = viridis_palette) +
  coord_polar(start = 0) +
  scale_x_continuous(breaks = c(0, pi/2, pi, 3*pi/2), labels = c("0", "pi/2", "pi", "3pi/2")) +
  labs(title = "Empirical Circular Rose Diagram", subtitle = "Relative turn angles distribution", x = NULL, y = "Step Count", fill = "Density") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

plot_trajectory <- ggplot(analytic_dataset, aes(x = x1_, y = y1_)) +
  geom_path(aes(color = theta), linewidth = 0.8, alpha = 0.8) +
  geom_point(aes(color = theta), size = 1.2) +
  scale_color_gradientn(colors = viridis_palette, breaks = c(0, pi/2, pi, 3*pi/2, 2*pi), labels = c("0", "pi/2", "pi", "3pi/2", "2pi")) +
  labs(title = "Fisher Spatial Movement Trajectory", subtitle = "Colored by turn angle", x = "UTM Easting (X)", y = "UTM Northing (Y)", color = "Angle") +
  theme_classic() + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

plot_temporal <- ggplot(analytic_dataset, aes(x = factor(step_hour), y = theta)) +
  geom_jitter(aes(color = log_step_t), alpha = 0.4, size = 1.5) +
  geom_boxplot(fill = "white", outlier.shape = NA, alpha = 0.7, color = "darkred") +
  scale_color_gradientn(colors = viridis_palette) +
  labs(title = "Turn Angle vs. Hour of Day", subtitle = "Tracking diurnal navigational cycles", x = "Hour of Day", y = "Turn Angle (Radians)", color = "Log Step") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

plot_joint_density <- ggplot(analytic_dataset, aes(x = theta, y = log_step_t)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", color = "white", linewidth = 0.1, alpha = 0.85) +
  scale_fill_gradientn(colors = magma_palette) +
  coord_polar(start = 0) +
  scale_x_continuous(breaks = c(0, pi/2, pi, 3*pi/2), labels = c("0", "pi/2", "pi", "3pi/2")) +
  labs(title = "Joint Turn Angle & Step Length Density", subtitle = "Circular-linear distribution signature", x = "Angle", y = "Log Step Length", fill = "Density Level") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

plot_env_polar <- ggplot(analytic_dataset, aes(x = theta, fill = forest_grp)) +
  geom_density(alpha = 0.65, color = "black", linewidth = 0.3) +
  coord_polar(start = 0) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_x_continuous(breaks = c(0, pi/2, pi, 3*pi/2), labels = c("0", "pi/2", "pi", "3pi/2")) +
  facet_wrap(~forest_grp) +
  labs(title = "Turn Angle Polar Density across Spatial Strata", subtitle = "Assessing spatial covariate effect", x = NULL, y = "Density", fill = "Stratum") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5), legend.position = "bottom")

plot_lag <- ggplot(analytic_dataset, aes(x = theta_lag, y = theta)) +
  geom_bin2d(bins = 30) +
  scale_fill_gradientn(colors = magma_palette) +
  labs(title = "Circular-Circular Autoregressive Lag Relation", subtitle = "Heading turn angle (t) vs. Turn angle (t-1)", x = "Theta (t-1)", y = "Theta (t)", fill = "Density") +
  theme_minimal() + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

ggsave("circular_rose_diagram.png", plot_rose, width = 6, height = 6, dpi = 300)
ggsave("spatial_trajectory_angles.png", plot_trajectory, width = 8, height = 6, dpi = 300)
ggsave("temporal_diurnal_patterns.png", plot_temporal, width = 10, height = 5, dpi = 300)
ggsave("circular_linear_joint_density.png", plot_joint_density, width = 6, height = 6, dpi = 300)
ggsave("circular_environmental_polar.png", plot_env_polar, width = 10, height = 5, dpi = 300)
ggsave("circular_circular_lag.png", plot_lag, width = 7, height = 6, dpi = 300)

dashboard_layout <- grid.arrange(plot_rose, plot_joint_density, plot_trajectory, plot_temporal, plot_env_polar, plot_lag, ncol = 2)
ggsave("fisher_comprehensive_dashboard.png", dashboard_layout, width = 14, height = 18, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Chronological Sequence Partitioning (80/20 Train-Test Train Matrix)
# ------------------------------------------------------------------------------
total_rows <- nrow(analytic_dataset)
split_index <- floor(0.80 * total_rows) 

train_data <- analytic_dataset[1:split_index, ]
test_data  <- analytic_dataset[(split_index + 1):total_rows, ]

feature_columns <- c("log_step_t", "log_step_tm1", "forest_cover", "elevation", "pop_density", "time_sin", "time_cos")

X_train <- as.matrix(train_data[, feature_columns])
y_train <- as.matrix(train_data$theta)

X_test  <- as.matrix(test_data[, feature_columns])
y_test  <- as.matrix(test_data$theta)

# ------------------------------------------------------------------------------
# 6. Structured Exploratory Data Quality & Behavioral Summary Tables (1 - 5)
# ------------------------------------------------------------------------------
table_partitioning <- data.frame(
  Partition_Type = c("Operational Training Set", "Held-out Testing Set", "Global Analytics Base"),
  Sample_Size_N  = c(nrow(X_train), nrow(X_test), total_rows),
  Percentage     = c((nrow(X_train)/total_rows)*100, (nrow(X_test)/total_rows)*100, 100),
  Features_D     = c(ncol(X_train), ncol(X_test), ncol(X_train))
)
print_table(table_partitioning, "Table 1: Data Partitioning & Dimension Metrics")

calc_summary <- function(vec, label) {
  data.frame(Covariate = label, Mean = mean(vec), Std_Dev = sd(vec), Median = median(vec), Minimum = min(vec), Maximum = max(vec))
}
table_descriptive <- rbind(
  calc_summary(analytic_dataset$sl_, "Step Length (sl_t)"),
  calc_summary(analytic_dataset$forest_cover, "Forest Cover Proxy (Standardized)"),
  calc_summary(analytic_dataset$elevation, "Elevation Proxy (Standardized)"),
  calc_summary(analytic_dataset$pop_density, "Population Density Proxy (Standardized)")
)
print_table(table_descriptive, "Table 2: Descriptive Profiles of Continuous Covariates")

calc_circular_linear_corr <- function(x, theta) {
  cos_theta <- cos(theta)
  sin_theta <- sin(theta)
  r_xc <- cor(x, cos_theta, use = "complete.obs")
  r_xs <- cor(x, sin_theta, use = "complete.obs")
  r_cs <- cor(cos_theta, sin_theta, use = "complete.obs")
  num  <- r_xc^2 + r_xs^2 - (2 * r_xc * r_xs * r_cs)
  den  <- 1 - r_cs^2
  return(num / den)
}

table_correlation <- data.frame(
  Linear_Predictor          = c("Step Length (sl_t)", "Forest Cover Proxy", "Elevation Proxy", "Population Density Proxy"),
  R2_Circular_Linear_Coef   = c(
    calc_circular_linear_corr(analytic_dataset$sl_, analytic_dataset$theta),
    calc_circular_linear_corr(analytic_dataset$forest_cover, analytic_dataset$theta),
    calc_circular_linear_corr(analytic_dataset$elevation, analytic_dataset$theta),
    calc_circular_linear_corr(analytic_dataset$pop_density, analytic_dataset$theta)
  )
) %>% mutate(Association_Strength_R = sqrt(R2_Circular_Linear_Coef))
print_table(table_correlation, "Table 3: Circular-Linear Correlation Profiles of Key Predictors")

daytime_theta  <- circular(analytic_dataset$theta[analytic_dataset$period == "Daytime"])
nighttime_theta <- circular(analytic_dataset$theta[analytic_dataset$period == "Nighttime"])

low_theta    <- circular(analytic_dataset$theta[analytic_dataset$forest_grp == "Low Cover"])
medium_theta <- circular(analytic_dataset$theta[analytic_dataset$forest_grp == "Medium Cover"])
high_theta   <- circular(analytic_dataset$theta[analytic_dataset$forest_grp == "High Cover"])

r_day   <- rho.circular(daytime_theta)
r_night <- rho.circular(nighttime_theta)
r_low   <- rho.circular(low_theta)
r_med   <- rho.circular(medium_theta)
r_high  <- rho.circular(high_theta)

approx_kappa <- function(r) {
  if (r < 0.53) return(2 * r + r^3 + (5 * r^5) / 6)
  else if (r < 0.85) return(-0.4 + 1.39 * r + 0.43 / (1 - r))
  else return(1 / (r^3 - 4 * r^2 + 3 * r))
}

table_diurnal <- data.frame(
  Feature_Profile       = c("Sample Step Count", "Circular Mean Direction (rad)", "Mean Resultant Length (R)", "Estimated VM Concentration (k)"),
  Daytime_06_00_18_00   = c(length(daytime_theta), mean(daytime_theta), r_day, approx_kappa(r_day)),
  Nighttime_18_00_06_00 = c(length(nighttime_theta), mean(nighttime_theta), r_night, approx_kappa(r_night))
)
print_table(table_diurnal, "Table 4: Baseline Circular Statistics Partitioned by Diurnal Period")

table_environmental <- data.frame(
  Feature_Profile     = c("Sample Step Count", "Circular Mean Direction (rad)", "Mean Resultant Length (R)", "Estimated VM Concentration (k)"),
  Low_Forest_Cover    = c(length(low_theta), mean(low_theta), r_low, approx_kappa(r_low)),
  Medium_Forest_Cover = c(length(medium_theta), mean(medium_theta), r_med, approx_kappa(r_med)),
  High_Forest_Cover   = c(length(high_theta), mean(high_theta), r_high, approx_kappa(r_high))
)
print_table(table_environmental, "Table 5: Baseline Circular Statistics Partitioned by Spatial Strata")

# ------------------------------------------------------------------------------
# 7. Analytical Loss Function & Metric Definitions
# ------------------------------------------------------------------------------
compute_vm_loss <- function(y_true, v1, v2, z_kappa, epsilon = 1e-6) {
  mu_pred <- atan2(v2, v1) %% (2 * pi)
  kappa_pred <- log(1 + exp(z_kappa)) + epsilon
  I0_val <- besselI(kappa_pred, nu = 0)
  loss_elements <- log(2 * pi) + log(I0_val) - (kappa_pred * cos(y_true - mu_pred))
  return(mean(loss_elements))
}

compute_mce <- function(y_true, y_pred) {
  return(mean(1 - cos(y_true - y_pred)))
}

# ------------------------------------------------------------------------------
# 8. Baseline Validation & Hyperparameter Bandwidth Sweep
# ------------------------------------------------------------------------------
nw_circular_linear_predict <- function(X_tr, y_tr, X_te, bandwidth = 0.5) {
  predictions <- numeric(nrow(X_te))
  for (i in 1:nrow(X_te)) {
    distances <- sqrt(rowSums((X_tr - matrix(X_te[i, ], nrow = nrow(X_tr), ncol = ncol(X_tr), byrow = TRUE))^2))
    weights <- exp(-(distances^2) / (2 * bandwidth^2))
    if (sum(weights) == 0) {
      predictions[i] <- mean(y_tr)
    } else {
      sum_sin <- sum(weights * sin(y_tr))
      sum_cos <- sum(weights * cos(y_tr))
      predictions[i] <- atan2(sum_sin, sum_cos) %% (2 * pi)
    }
  }
  return(predictions)
}

cat("[EVALUATION] Initiating Hyperparameter sweep over Nadaraya-Watson bandwidths (h)...\n")
bandwidths <- seq(0.1, 1.0, by = 0.15)
sweep_results <- data.frame(Bandwidth = bandwidths, MCE = NA)

for (j in 1:length(bandwidths)) {
  temp_preds <- nw_circular_linear_predict(X_train, y_train, X_test, bandwidth = bandwidths[j])
  sweep_results$MCE[j] <- compute_mce(y_test, temp_preds)
}

optimal_idx <- which.min(sweep_results$MCE)
opt_bandwidth <- sweep_results$Bandwidth[optimal_idx]
baseline_mce   <- sweep_results$MCE[optimal_idx]

baseline_preds <- nw_circular_linear_predict(X_train, y_train, X_test, bandwidth = opt_bandwidth)

daytime_indices <- which(test_data$period == "Daytime")
nighttime_indices <- which(test_data$period == "Nighttime")

day_mce <- compute_mce(y_test[daytime_indices], baseline_preds[daytime_indices])
night_mce <- compute_mce(y_test[nighttime_indices], baseline_preds[nighttime_indices])

# ------------------------------------------------------------------------------
# 9. Structured Optimization & Model Performance Scorecards (Tables 6 - 7)
# ------------------------------------------------------------------------------
table_bandwidth_sweep <- data.frame(
  Sweep_Index  = 1:nrow(sweep_results),
  Bandwidth_h  = sweep_results$Bandwidth,
  Validation_MCE = sweep_results$MCE,
  Status       = ifelse(1:nrow(sweep_results) == optimal_idx, "OPTIMAL CHOICE *", "Sub-optimal")
)
print_table(table_bandwidth_sweep, "Table 6: Nadaraya-Watson Kernel Bandwidth Sensitivity Analysis")

# --- TABLE 7: Model Scorecard (Comment 1: Unverified vM-NN row removed) ---
table_scorecard <- data.frame(
  Model_Configuration_Segment   = c(
    paste0("Overall Nadaraya-Watson Baseline (Tuned h = ", opt_bandwidth, ")"),
    "Nadaraya-Watson Baseline (Daytime Split Only)",
    "Nadaraya-Watson Baseline (Nighttime Split Only)"
  ),
  Mean_Circular_Error_MCE      = c(baseline_mce, day_mce, night_mce),
  Reference_Source_Type         = c(
    "Tuned Empirical Outcome",
    "Tuned Diurnal Sub-slice",
    "Tuned Diurnal Sub-slice"
  )
)
print_table(table_scorecard, "Table 7: Predictive Performance Summary & Baseline Scorecard")

cat("  *Note: Chosen dynamically via local bandwidth minimization search loop.\n")
cat("=========================================================================\n")
