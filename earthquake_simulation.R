library(lubridate)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}
library(ggplot2)
library(scales)

setwd("C:/Users/georg/OneDrive/Desktop")
earthquakes_data <- read.csv("Japan earthquakes data.csv")
magnitudes <- earthquakes_data$Magnitude
event_times_strings <- earthquakes_data[,1]
event_times <- dmy_hm(event_times_strings)
earliest_event <- min(event_times)
hours_after_earliest <- as.numeric(difftime(event_times, earliest_event, units = "hours"))
times <- sort(hours_after_earliest)
magnitudes <- rev(magnitudes)
T <- max(times) + 1
neg_log_likelihood <- function(params, times, magnitudes, T) {
  mu <- params[1]
  alpha <- params[2]
  beta <- params[3]
  gamma <- params[4]
  n <- length(times)
  
  log_likelihood_val <- 0
  for (i in 1:n) {
    if (i == 1) {
      intensity <- mu
    } else {
      intensity <- mu + sum(alpha * exp(-beta * (times[i] - times[1:(i-1)])) * exp(gamma * magnitudes[1:(i-1)]))
    }
    log_likelihood_val <- log_likelihood_val + log(intensity)
  }
  
  log_likelihood_val <- log_likelihood_val - mu * T
  for (j in 1:n) {
    log_likelihood_val <- log_likelihood_val - (alpha / beta) * (1 - exp(-beta * (T - times[j]))) * exp(gamma * magnitudes[j])
  }
  
  return(-log_likelihood_val) 
}
initial_params <- c(0.15, 0.5, 1, 1)
if (!requireNamespace("optimx", quietly = TRUE)) {
  install.packages("optimx")
}
library(optimx)
result <- optimx::optimx(par = initial_params, fn = neg_log_likelihood, times = times, magnitudes = magnitudes, T = T, method = "L-BFGS-B", lower = c(0, 0, 0, 0))

mu_est <- result$p1
alpha_est <- result$p2
beta_est <- result$p3
gamma_est <- result$p4

cat("Estimated Mu:", mu_est, "\n")
cat("Estimated Alpha:", alpha_est, "\n")
cat("Estimated Beta:", beta_est, "\n")
cat("Estimated Gamma:", gamma_est, "\n")

intensity_function <- function(t, times, magnitudes, mu, alpha, beta, gamma) {
  intensity <- mu
  for (i in 1:length(times)) {
    if (times[i] < t) {
      intensity <- intensity + alpha * exp(-beta * (t - times[i])) * exp(gamma * magnitudes[i])
    }
  }
  return(intensity)
}
plot_data <- data.frame(
  times = times,
  magnitudes = magnitudes
)
time_seq <- seq(0, T, by = 0.01)
intensity_values <- sapply(time_seq, function(t) intensity_function(t, times, magnitudes, mu_est, alpha_est, beta_est, gamma_est))
intensity_data <- data.frame(
  time = time_seq,
  intensity = intensity_values
)

ggplot() +
  geom_point(data = plot_data, aes(x = times, y = magnitudes), color = "blue", size = 2, shape = 16) +
  geom_line(data = intensity_data, aes(x = time, y = intensity), color = "red", size = 1) +
  labs(
       x = "Time (Hours after first quake)",
       y = "Magnitude / Intensity") +
  theme(axis.text = element_text(size = 15),      
        axis.title = element_text(size = 20),     
        plot.title = element_text(size = 25))    
time_seq <- seq(490, 510, by = 0.01)
intensity_values <- sapply(time_seq, function(t) intensity_function(t, times, magnitudes, mu_est, alpha_est, beta_est, gamma_est))


intensity_data <- data.frame(
  time = time_seq,
  intensity = intensity_values
)

ggplot() +
  geom_point(data = plot_data, aes(x = times, y = magnitudes), color = "blue", size = 2, shape = 16) +
  geom_line(data = intensity_data, aes(x = time, y = intensity), color = "red", size = 1) +
  scale_y_continuous(
    name = "Magnitude/Intensity",
  ) +
  coord_cartesian(xlim = c(490, 510)) +
  labs(
       x = "Time (Hours after first quake)") +
  theme(axis.text = element_text(size = 15),      
        axis.title = element_text(size = 20),     
        plot.title = element_text(size = 25))     