library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)


df <- read.xlsx("mito_stats.xlsx")


df_long <- df %>%
  pivot_longer(cols = c(rrnL, rrnS, tRNA, Control_Region, Total_Length, PCG), 
               names_to = "Region", values_to = "Length")

# Replace underscores with spaces in Region labels
df_long$Region <- str_replace_all(df_long$Region, "_", " ")

# Define color mapping
df_long <- df_long %>%
  mutate(Color = case_when(
    Species == "Phyllium philippinicum" ~ "darkgreen"
  ))

ggplot(df_long, aes(x = Region, y = Length / 1000)) +
geom_boxplot(fill = "grey", color = "black") +
  geom_point(aes(color = Color), size = 3) +
  # Add labels only for Control Region and Total Length for P. philippinicum
  geom_text(
    data = df_long %>%
      filter(Species == "P. philippinicum", Region %in% c("Control Region", "Total Length")),
    aes(label = Species),
    vjust = -0.5,
    size = 3,
    color = "darkgreen"
  ) +
  scale_color_identity() +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2)
  ) +
  labs(
    y = "Length (Kb)",
    x = "Feature",
    title = "Boxplot of Length of the Mitochondrial Features in Euphasmatodea"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    panel.border = element_blank(),  # Ensures only axes, not full box
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
