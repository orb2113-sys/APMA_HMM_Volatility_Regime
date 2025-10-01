# Create synthetic HMM data for RGB plotting demo
cat("Creating synthetic RGB data for demo...\n")

n_days <- 200
dates <- seq(as.Date("2022-01-01"), length.out = n_days, by = "day")
prices <- 100 * cumprod(1 + rnorm(n_days, 0, 0.02))

# Generate synthetic probabilities with regime persistence
prob_R <- prob_G <- prob_B <- numeric(n_days)
current_state <- sample(1:3, 1)

for(i in 1:n_days) {
  if(runif(1) < 0.05) current_state <- sample(1:3, 1)  # 5% chance of switch
  probs <- c(0.1, 0.1, 0.1)
  probs[current_state] <- 0.8
  prob_R[i] <- probs[1]; prob_G[i] <- probs[2]; prob_B[i] <- probs[3]
}

write.csv(data.frame(date=dates, price=prices, prob_R=prob_R, prob_G=prob_G, prob_B=prob_B),
          "hmm_rgb_data.csv", row.names=FALSE)

cat("Created synthetic hmm_rgb_data.csv for RGB plotting demo\n")