setwd("C:/Users/kevin/Documents/KM/UMichigan/Submissions")

source("scripts.R")

back_pain <- readRDS("data/Back_Pain.Rds")
pain_combos <- expand.grid(unique(back_pain$x1), unique(back_pain$x2),
                           unique(back_pain$x3), unique(back_pain$x4)) |>
  rename(x1 = Var1, x2 = Var2, x3 = Var3, x4 = Var4)

#####Run the below lines once to save output###################################
#set.seed(2515) ## Seed can be anything
#pain_predictions <- runBootstrapRegression(var_index = 4, data = back_pain, 
#                                           runs = 1000)
#saveRDS(pain_predictions, "data/output/pain_wk_predictions.Rds")
###############################################################################

pain_predictions <- readRDS("data/output/pain_wk_predictions.Rds")
full_data <- getBubblePlotData(regression_data = pain_predictions,
                               data = back_pain,
                               var_index = 4)

full_data <- full_data |>
  mutate(
    x1_name = ifelse(x1 == 1, "S", "L"),
    x3_name = ifelse(x3 == 1, "AD", "PI"),
    x2_name = factor(
      case_when(
        x2 == 1 ~ "W",
        x2 == 2 ~ "S",
        x2 == 3 ~ "B"
      ),
      levels = c("W", "S", "B")
    ),
    x4_name = factor(
      case_when(
        x4 == 1 ~ "W",
        x4 == 2 ~ "S",
        x4 == 3 ~ "SI",
        x4 == 4 ~ "MODI",
        x4 == 5 ~ "MARI",
        x4 == 6 ~ "CR"
      ),
      levels = c("W", "S", "SI", "MODI",
                 "MARI", "CR")
    )
  )

full_data$combination_2_3 <- paste0(
  "(",
  full_data$x2_name,
  ",",
  full_data$x3_name,
  ")"
)
p1 <- full_data |>
  filter(x1_name == "S") |>
  ggplot(aes(x = combination_2_3, y = x4_name)) +
  geom_point(aes(
    size = ifelse(n == 0, NA, n * 100), color = x2),
    shape = 21, colour = "black",
    fill = "white", stroke = .5
  ) +
  scale_size_continuous(range = c(1, 20)) +
  geom_point(aes(
    x = combination_2_3,
    y = x4_name,
    size =
      ifelse(n == 0 | n <= 500,
             NA, .1
      )
  )) +
  geom_text(aes(label = ifelse(n == 0, "", n)),
            vjust = -1.2, size = 2
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 13),
    axis.text.y = element_text(
      angle = 45, vjust = 0.5,
      hjust = 1, size = 8
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    x = latex2exp::TeX("$(X_2,X_3)$ with $X_1=$ S (short)"),
    y = latex2exp::TeX("Predicted Back Pain ($\\hat{X}_4$)")
  ) |>
  suppressWarnings()

p2 <- full_data |>
  filter(x1_name == "L") |>
  ggplot(aes(x = combination_2_3, y = x4_name)) +
  geom_point(aes(
    size = ifelse(n == 0, NA, n * 100), color = x2),
    shape = 21, colour = "black",
    fill = "white", stroke = .5
  ) +
  scale_size_continuous(range = c(1, 20)) +
  geom_point(aes(
    x = combination_2_3,
    y = x4_name,
    size =
      ifelse(n == 0 | n <= 500,
             NA, .1
      )
  )) +
  geom_text(aes(label = ifelse(n == 0, "", n)),
            vjust = -1.2, size = 2
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 13),
    axis.text.y = element_text(
      angle = 45, vjust = 0.5,
      hjust = 1, size = 8
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    x = latex2exp::TeX("$(X_2,X_3)$ with $X_1=$ L (long)"),
    y = latex2exp::TeX("Predicted Back Pain ($\\hat{X}_4$)")
  )

gridExtra::grid.arrange(p1, p2, nrow = 2, ncol = 1)