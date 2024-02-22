# Load the required libraries
library(ggplot2)
library(openxlsx)

## KEN-SHE results
# 9v 98.8 (91.3-99.8)
# 2v 97.5 (90.0-99.4)

 ## 9-valent Study V503-001
# 9v 96.7 (80.9-99.8)

## 2-valent study, Clinical Review Memo - Cervarix, October 15, 2009
# 2v 98.67 (74.4-100)

# Set the median and confidence interval values
median_val <- 0.988 # note this is actually the mean
lower_ci <- 0.913
upper_ci <- 0.998

# Calculate the mean and standard deviation of the beta distribution
alpha_plus_beta <- ((median_val - 1/3)/(median_val + 2/3))^(-1)
beta_dist_mean <- median_val
beta_dist_std <- (upper_ci - lower_ci)/(2*1.96)

# Calculate the alpha and beta parameters of the beta distribution
beta_dist_alpha <- beta_dist_mean^2*(1-beta_dist_mean)/beta_dist_std^2 - beta_dist_mean
beta_dist_beta <- beta_dist_alpha*(1/beta_dist_mean - 1)

## Optimizing alpha and beta to fit as close as possible to KEN-SHE data 

result <- data.frame()


for (alpha in seq(25, 27, 0.01)) {
      for (beta in seq(0.3, 0.7, 0.005)) {
            lower_ci <- qbeta(0.025, alpha, beta)
            upper_ci <- qbeta(0.975, alpha, beta)
            median <- qbeta(0.5, alpha, beta)
            
            temp <- data.frame("lower" = lower_ci, "upper" = upper_ci, "med" = median, "alpha" = alpha, "beta" = beta)
            
            result <- bind_rows(result, temp)
      }
}

# filter and finding the optimal parameters 

result %>% 
      filter(lower >= 0.90, #lower <= 0.8, 
             upper >= 0.98, upper <= 1, 
             med >= 0.98) %>%  # & med < 0.99) %>% 
      arrange(lower) %>% View()

## 2-valent study, Clinical Review Memo - Cervarix, October 15, 2009
# 2v 98.67 (74.4-100)

## Setting the new alpha and beta 
beta_dist_alpha <- 26.94
beta_dist_beta <- 0.615

# Calculate the lower and upper 95% CI
lower_ci <- qbeta(0.025, beta_dist_alpha, beta_dist_beta)
upper_ci <- qbeta(0.975, beta_dist_alpha, beta_dist_beta)
median <- qbeta(0.5, beta_dist_alpha, beta_dist_beta)

# Generate random values from the beta distribution
set.seed(123)
random_vals <- rbeta(1000000, beta_dist_alpha, beta_dist_beta)


# Create a data frame with the random values
df <- data.frame(x = random_vals*100)

# Create a data frame with the random values
df <- data.frame(x = random_vals)

# Setting the alpha and beta for the 2-dose
beta_dist_alpha <- 26.94
beta_dist_beta <- 0.615

# Calculate the lower and upper 95% CI
lower_ci <- qbeta(0.025, beta_dist_alpha, beta_dist_beta)
upper_ci <- qbeta(0.975, beta_dist_alpha, beta_dist_beta)
median <- qbeta(0.5, beta_dist_alpha, beta_dist_beta)

# Generate random values from the beta distribution
set.seed(123)
random_vals <- rbeta(1000000, beta_dist_alpha, beta_dist_beta)


# Create a data frame with the random values
df <- data.frame(x = random_vals*100)

# Create a data frame with the random values
df <- data.frame(x = random_vals)

## Setting the new alpha and beta 
beta_dist_alpha <- 26.94
beta_dist_beta <- 0.615

# Calculate the lower and upper 95% CI
lower_ci <- qbeta(0.025, beta_dist_alpha, beta_dist_beta)
upper_ci <- qbeta(0.975, beta_dist_alpha, beta_dist_beta)
median <- qbeta(0.5, beta_dist_alpha, beta_dist_beta)

# Generate random values from the beta distribution
set.seed(123)
random_vals <- rbeta(1000000, beta_dist_alpha, beta_dist_beta)


# Create a data frame with the random values
df <- data.frame(x = random_vals*100)

# Create a data frame with the random values
df <- data.frame(x = random_vals)

                        # Setting the alpha and beta for the 2-dose
                        beta_dist_alpha <- 15.66
                        beta_dist_beta <- 0.82
                        
                        # Calculate the lower and upper 95% CI
                        lower_ci <- qbeta(0.025, beta_dist_alpha, beta_dist_beta)
                        upper_ci <- qbeta(0.975, beta_dist_alpha, beta_dist_beta)
                        median <- qbeta(0.5, beta_dist_alpha, beta_dist_beta)
                        
                        # Generate random values from the beta distribution
                        set.seed(123)
                        random_vals <- rbeta(1000000, beta_dist_alpha, beta_dist_beta)
                        
                        # Create a data frame with the random values
                        df_9v_2d <- data.frame(x = random_vals)

bounded <- df %>% 
      bind_cols(df_9v_2d) %>% 
      filter(`x...2` >= `x...1`) 

dose1_9v <- bounded %>% select(`x...1`) %>% rename(x=`x...1`)
dose2_9v <- bounded %>% select(`x...2`) %>% rename(x=`x...2`)

# Generate the plot using ggplot2
ggplot(df, aes(x = x)) +
      stat_density(aes(y = ..scaled..), fill = "#69b3a2", alpha = 0.6) +
      labs(title = "Beta Probability Distribution", x = "Value", y = "Probability") +
      theme_minimal() +
      xlim(0.9,1)

# generate 100 random values and save to excel
set.seed(123)
random_vals <- rbeta(100, beta_dist_alpha, beta_dist_beta)

OUT <- createWorkbook()

addWorksheet(OUT, "VaxEfficacyRandVal")

writeData(OUT, x = dose2_9v, startRow = 1, sheet="VaxEfficacyRandVal")

saveWorkbook(OUT, "/Users/clh89/MATLAB/Projects/ccTreatment_KZN/Params/VaxEfficacyRandVal_2dose_nonavalent_bounded.xlsx", overwrite = TRUE)


OUT <- createWorkbook()

addWorksheet(OUT, "VaxEfficacyRandVal")

writeData(OUT, x = dose1_9v, startRow = 1, sheet="VaxEfficacyRandVal")

saveWorkbook(OUT, "/Users/clh89/MATLAB/Projects/ccTreatment_KZN/Params/VaxEfficacyRandVal_1dose_nonavalent_bounded.xlsx", overwrite = TRUE)
