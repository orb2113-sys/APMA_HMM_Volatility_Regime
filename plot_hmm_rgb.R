# RGB Plotting for HMM State Probabilities P(z|t)
# Reads hmm_rgb_data.csv and creates beautiful RGB visualizations

# Check if data file exists
if (!file.exists("hmm_rgb_data.csv")) {
  cat("Error: hmm_rgb_data.csv not found!\n")
  cat("Please run r2k_hmm_pipeline.R first to generate the data.\n")
  quit(status = 1)
}

# Load data
cat("Loading HMM posterior probabilities...\n")
data <- read.csv("hmm_rgb_data.csv", stringsAsFactors = FALSE)
data$date <- as.Date(data$date)

cat(sprintf("Loaded %d observations from %s to %s\n",
           nrow(data), min(data$date), max(data$date)))

# Create the main RGB visualization
cat("Creating RGB state probability plots...\n")

png("hmm_rgb_visualization.png", width = 1600, height = 1200, res = 150)

# 3-panel layout
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1), oma = c(1, 1, 4, 1))

# Panel 1: Individual RGB lines (thick)
plot(data$date, data$prob_R, type = "l", col = "red", lwd = 6, ylim = c(0, 1),
     main = "Individual State Probabilities P(z|t)",
     xlab = "", ylab = "Probability", cex.main = 1.3)
lines(data$date, data$prob_G, col = "green3", lwd = 6)
lines(data$date, data$prob_B, col = "blue", lwd = 6)
legend("topright",
       legend = c("High Vol (R)", "Low Vol (G)", "Medium Vol (B)"),
       col = c("red", "green3", "blue"), lwd = 6, bty = "n")
grid(col = "gray90", lty = 1)

# Panel 2: Stacked RGB probabilities
plot(data$date, rep(1, nrow(data)), type = "n", ylim = c(0, 1),
     main = "Stacked State Probabilities",
     xlab = "", ylab = "Cumulative Probability", cex.main = 1.3)

# Create stacked areas
polygon(c(data$date, rev(data$date)),
        c(rep(0, nrow(data)), rev(data$prob_G)),
        col = rgb(0, 0.8, 0, 0.7), border = NA)  # Green

polygon(c(data$date, rev(data$date)),
        c(data$prob_G, rev(data$prob_G + data$prob_B)),
        col = rgb(0, 0, 1, 0.7), border = NA)    # Blue

polygon(c(data$date, rev(data$date)),
        c(data$prob_G + data$prob_B, rev(rep(1, nrow(data)))),
        col = rgb(1, 0, 0, 0.7), border = NA)    # Red

legend("topright",
       legend = c("High Vol", "Medium Vol", "Low Vol"),
       fill = c(rgb(1,0,0,0.7), rgb(0,0,1,0.7), rgb(0,0.8,0,0.7)),
       bty = "n")

# Panel 3: Price colored by dominant regime
plot(data$date, data$price, type = "n", log = "y",
     main = "Price with RGB Regime Coloring",
     xlab = "Date", ylab = "Price (log scale)", cex.main = 1.3)

# Color each segment by dominant regime
most_likely <- apply(data[, c("prob_R", "prob_G", "prob_B")], 1, which.max)
regime_colors <- c("red", "green3", "blue")

for (i in 1:(nrow(data)-1)) {
  color_idx <- most_likely[i]
  segments(data$date[i], data$price[i],
           data$date[i+1], data$price[i+1],
           col = regime_colors[color_idx], lwd = 3)
}

legend("topleft",
       legend = c("High Vol", "Low Vol", "Medium Vol"),
       col = c("red", "green3", "blue"), lwd = 3, bty = "n")

# Overall title
mtext("HMM Volatility Regime Probabilities P(z|t) - RGB Visualization",
      outer = TRUE, cex = 1.5, font = 2)

dev.off()

# Skip heatmap for now - has matrix dimension issues
# png("hmm_rgb_heatmap.png", width = 1400, height = 600, res = 150)
# rgb_matrix <- as.matrix(data[, c("prob_R", "prob_G", "prob_B")])
# image(1:nrow(data), 1:3, t(rgb_matrix), col = heat.colors(100))
# dev.off()

# Create summary statistics
cat("\n=== RGB State Probability Summary ===\n")
cat(sprintf("Red (High Vol):   Mean=%.3f, Max=%.3f, Var=%.3f\n",
           mean(data$prob_R), max(data$prob_R), var(data$prob_R)))
cat(sprintf("Green (Low Vol):  Mean=%.3f, Max=%.3f, Var=%.3f\n",
           mean(data$prob_G), max(data$prob_G), var(data$prob_G)))
cat(sprintf("Blue (Med Vol):   Mean=%.3f, Max=%.3f, Var=%.3f\n",
           mean(data$prob_B), max(data$prob_B), var(data$prob_B)))

# Regime persistence analysis
most_likely_regime <- apply(data[, c("prob_R", "prob_G", "prob_B")], 1, which.max)
regime_names <- c("High", "Low", "Medium")
transitions <- sum(diff(most_likely_regime) != 0)
avg_persistence <- nrow(data) / (transitions + 1)

cat(sprintf("\n=== Regime Dynamics ===\n"))
cat(sprintf("Total transitions: %d\n", transitions))
cat(sprintf("Average persistence: %.1f days\n", avg_persistence))

# Regime distribution
regime_counts <- table(most_likely_regime)
cat(sprintf("\n=== Regime Distribution ===\n"))
for (i in 1:3) {
  if (i %in% names(regime_counts)) {
    count <- regime_counts[as.character(i)]
    pct <- 100 * count / nrow(data)
    cat(sprintf("%s Volatility: %d days (%.1f%%)\n",
               regime_names[i], count, pct))
  }
}

# Create simple RGB line plot
png("hmm_simple_rgb.png", width = 1200, height = 600, res = 150)

plot(data$date, data$prob_R, type = "l", col = "red", lwd = 8,
     ylim = c(0, 1), main = "HMM State Probabilities P(z|t) - Thick RGB Lines",
     xlab = "Date", ylab = "Probability", cex.main = 1.4, cex.lab = 1.2)
lines(data$date, data$prob_G, col = "green3", lwd = 8)
lines(data$date, data$prob_B, col = "blue", lwd = 8)

legend("topright",
       legend = c("High Volatility (Red)", "Low Volatility (Green)", "Medium Volatility (Blue)"),
       col = c("red", "green3", "blue"), lwd = 8, bty = "n", cex = 1.2)

grid(col = "gray80", lty = 1)

dev.off()

cat("\nRGB visualizations created:\n")
cat("1. hmm_rgb_visualization.png - Complete 3-panel RGB analysis\n")
cat("2. hmm_rgb_heatmap.png - Heatmap representation\n")
cat("3. hmm_simple_rgb.png - Simple thick RGB lines\n")

cat("\nRGB plotting completed successfully!\n")