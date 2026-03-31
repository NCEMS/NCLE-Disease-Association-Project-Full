cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("type: af\n")
cat("plddt_thresh: 70\n")

# Output file names
(fileName = paste("sim_mis_prop_type_af_pthr_70"))

# Libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lawstat)

# Load simulation results
cg_simulation_results <- read_excel("Data/cg_simulation_results_2026_03_09.xlsx")

# Subset only variables needed for the permutation test and compute the difference
sub_df <- cg_simulation_results[, c("gene", "misfolding_prop", "gene_ctrl", "misfolding_prop_ctrl")] %>%
  mutate(diff = misfolding_prop - misfolding_prop_ctrl)

# Compute observed test statistic
# Mean/Median paired difference in misfolding propensity
obs_diff <- mean(sub_df$diff)

obs_median_diff <- median(sub_df$diff)

# Brunner-Munzel test
bm_test <- brunner.munzel.test(sub_df$misfolding_prop, sub_df$misfolding_prop_ctrl, alternative = "greater")

# 95% Bootstrap CI of the Median (10^6 iterations)
set.seed(123)

B <- 1e6

# Bootstrap function
bootstrap_median <- function(x, B) {
  n <- length(x)
  replicate(B, median(sample(x, n, replace = TRUE)))
}

boot_x <- bootstrap_median(sub_df$misfolding_prop, B)
boot_y <- bootstrap_median(sub_df$misfolding_prop_ctrl, B)

# 95% CI
ci_x <- quantile(boot_x, c(0.025, 0.975))
ci_y <- quantile(boot_y, c(0.025, 0.975))

# Observed medians
median_x <- median(sub_df$misfolding_prop)
median_y <- median(sub_df$misfolding_prop_ctrl)

median_x
ci_x

median_y
ci_y

# Save results
results_df <- data.frame(
  BM_Statistic             = bm_test$statistic,
  BM_p_value               = bm_test$p.value,
  Median_Entangled         = median_x,
  CI_Entangled_Lower       = ci_x[1],
  CI_Entangled_Upper       = ci_x[2],
  Median_Control           = median_y,
  CI_Control_Lower         = ci_y[1],
  CI_Control_Upper         = ci_y[2]
)

write.csv(
  results_df,
  file = paste0("Results/Dataframes/", fileName, "_results.csv"),
  row.names = FALSE
)

# Plot
plot_df <- sub_df %>%
  select(misfolding_prop, misfolding_prop_ctrl) %>%
  pivot_longer(everything(), names_to = "Group", values_to = "Misfolding") %>%
  mutate(
    Group = recode(
      Group,
      misfolding_prop      = "Natively entangled\nand disease-associated",
      misfolding_prop_ctrl = "Natively non-entangled\nand non-disease-associated"
    )
  )

p <- ggplot(plot_df, aes(x = Group, y = Misfolding)) +
  geom_violin(
    fill = "grey80", color = "black",
    width = 0.9,
    trim = TRUE,
    adjust = 1.0,
    scale = "width"
  ) +
  geom_point(
    position = position_jitter(width = 0.06, height = 0),
    size = 2.5, alpha = 0.85
  ) +
  # thick median bar
  stat_summary(
    fun = median,
    geom = "errorbar",
    aes(ymin = after_stat(y), ymax = after_stat(y)),
    width = 0.25,
    linewidth = 2.2
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = "Misfolding propensity", x = NULL) +
  theme_classic(base_size = 16) +
  theme(
    axis.line  = element_line(linewidth = 1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(6, "pt"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 18),
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  )

p

# Save plot
ggsave(
  filename = paste0("Results/Plots/", fileName, "_plot.svg"),
  plot = p,
  width = 6,
  height = 6
)

ggsave(
  filename = paste0("Results/Plots/", fileName, "_plot.pdf"),
  plot = p,
  width = 6,
  height = 6
)
