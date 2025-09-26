# Use decay model to estimate percent of carbon in fast vs. slow decaying pools
# and associated decay rates

library(dplyr)
library(ggplot2)

# colors for 5 cover crops:
covercols <- c("chocolate", "goldenrod", "forestgreen", "darkcyan", "darkorchid")
# colors for 2 moisture levels
mcols <- c("blue", "red")




# read CO2 data
dat <- read.csv("raw-data/CO2_long.csv")
dat$C_remaining <- 100*(dat$Initial-dat$CO2Cumulative)/dat$Initial

# days on which respiration was measured
days <- read.csv("raw-data/respiration-measurement-days.csv")

# set factors
dat$Residue.type <- as.factor(dat$Residue.type)
dat$Residue.type <- factor(dat$Residue.type, levels(dat$Residue.type)[c(2,5,1,4,3)])
dat$Moisture.treatment <- as.factor(dat$Moisture.treatment)
dat$Sample.number <- as.factor(dat$Sample.number)
dat$Date <- as.numeric(dat$Date)

# Plot cumulative CO2 emissions, each microcosm gets its own curve
p <- ggplot(dat, aes(x = Date, y = C_remaining, group=Sample.number, color = Residue.type)) +
  geom_line(aes(linetype=Moisture.treatment)) +
  scale_color_manual(values=covercols) +
  labs(x = "Days after glucose addition",
       y = bquote("Percent C remaining"),
       linetype = "Moisture treatment",
       color = "Residue type") +
  theme_bw()
p


# create a dataframe for coefficients 
df <- data.frame(Sample.number = levels(dat$Sample.number),
                 C1 = rep(NA, length(levels(dat$Sample.number))),
                 k1 = rep(NA, length(levels(dat$Sample.number))),
                 C2 = rep(NA, length(levels(dat$Sample.number))),
                 k2 = rep(NA, length(levels(dat$Sample.number))),
                 Residue.type = rep(NA, length(levels(dat$Sample.number))),
                 Moisture.treatment = rep(NA, length(levels(dat$Sample.number))))

dat$C_remaining_pred <- rep(NA, dim(dat)[1])

# estimate percent of C in fast pool (C1) and slow pool (C2), with mineralization rates k1 and k2 respectively
for(i in 1:length(levels(dat$Sample.number))){
  # subset to sample number i
  dat_subset <- dat[which(dat$Sample.number == df$Sample.number[i]),]
  # estimate coefficients using non-linear least squares:
  fit <- nls(
    C_remaining ~ C1 * (exp(-k1 * Date)) + C2 * (exp(-k2 * Date)),
    data = dat_subset, 
    start = list(C1 = 3, k1 = 0.01, C2 = 97, k2 = 0.001)
  )
  # add coefficients to df
  df[i , c("C1", "k1", "C2", "k2")] <- coef(fit)
  
  # assign treatments based on original dataset
  df$Residue.type[i] <- as.character(dat_subset$Residue.type[1])
  df$Moisture.treatment[i] <- as.character(dat_subset$Moisture.treatment[1])
  
  # add predicted C remaining to original dataset
  dat$C_remaining_pred[which(dat$Sample.number == df$Sample.number[i])] <- predict(fit)
}
# export averaged parameter estimates
write.csv(df, "model-output/decaymodel-twopool.csv")





# get averaged data by treatments
averaged_data <- df %>%
  group_by(Residue.type, Moisture.treatment) %>%
  summarise(
    C1 = mean(C1, na.rm = TRUE),
    C2 = mean(C2, na.rm = TRUE),
    k1 = mean(k1, na.rm = TRUE),
    k2 = mean(k2, na.rm = TRUE)
  )
# export averaged parameter estimates
write.csv(averaged_data, "model-output/decaymodel-twopool_averaged.csv")

################# plot C pools

# transform data in order to create stacked barplot of C pools
averaged_data_Cpools <- averaged_data %>% 
  select(Residue.type, Moisture.treatment, C1, C2) %>% 
  tidyr::pivot_longer(
    cols = starts_with("C"),
    names_to = "Cpool",
    values_to = "value"
  )

# stacked barchart of C pools
p <- ggplot(averaged_data_Cpools, aes(x = Residue.type, y = value, fill = Cpool)) +
  facet_wrap(~Moisture.treatment) + 
  geom_col(position = "stack", color = "black") + # geom_col for bar plots, position="stack" for stacking
  labs(
    title = "",
    x = "Residue type",
    y = "SOC pool size (% of initial)",
    fill = "C pool"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels if needed
  ) +
  scale_fill_manual(values = c("red", "orange")) 
p
ggsave(p, file="figures/C pools.png",width = 8, height = 4)


################# plot decay rates

# plot decay rates : k1
p <- ggplot(averaged_data, aes(x = Residue.type, y = k1, fill = Residue.type)) +
  geom_bar(stat = "identity") + 
  facet_grid(.~Moisture.treatment) + 
  labs(
    title = "",
    x = "Residue type",
    y = "SOC pool size (% of initial)",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels if needed
  ) +
  scale_fill_manual(values = covercols) 
p
ggsave(p, file="figures/C rates_k1.png",width = 8, height = 4)

# plot decay rates : k2
p <- ggplot(averaged_data, aes(x = Residue.type, y = k2, fill = Residue.type)) +
  geom_bar(stat = "identity") + 
  facet_grid(.~Moisture.treatment) + 
  labs(
    title = "",
    x = "Residue type",
    y = "SOC pool size (% of initial)",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels if needed
  ) +
  scale_fill_manual(values = covercols) 
p
ggsave(p, file="figures/C rates_k2.png",width = 8, height = 4)





# Summarize data to calculate treatment averages
cumulative_avg <- dat %>%
  group_by(Date, Residue.type, Moisture.treatment) %>%
  summarise(mean_C_remaining = mean(C_remaining, na.rm = TRUE),
            se_C_remaining = sd(C_remaining, na.rm = TRUE)/2,
            mean_C_remaining_pred = mean(C_remaining_pred, na.rm = TRUE),
            se_C_remaining_pred = sd(C_remaining_pred, na.rm = TRUE)/2, .groups = "drop")


# Plot cumulative CO2 emissions averages
p <- ggplot(cumulative_avg[which(cumulative_avg$Date %in% days$Days),], aes(x = Date, y = mean_C_remaining, color = Residue.type)) +
  geom_errorbar(aes(ymin = mean_C_remaining - se_C_remaining, ymax = mean_C_remaining + se_C_remaining), width=2, size=0.5) +
  geom_point(size=0.5) + 
  facet_grid(.~Moisture.treatment) +
  labs(x = "Day", y = bquote("C remaining (%"),
       linetype = "Moisture treatment", 
       color = "Residue type") +
  theme_bw() +
  # add predicted values
  geom_line(inherit.aes = F, data = cumulative_avg, aes(x = Date, y = mean_C_remaining_pred, color = Residue.type))
p

ggsave(p, file="figures/C remaining.png",width = 8, height = 4)


