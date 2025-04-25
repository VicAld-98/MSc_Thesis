#CDSPT data
codon_df <- read.csv("cdspt.csv", header = T)

library(ggplot2)
library(tidyr)

# Reshape the dataframe to long format
codon_long <- pivot_longer(codon_df, 
                           cols = -Codon, 
                           names_to = "Sample", 
                           values_to = "Value")

codon_long$Species <- gsub("_", ". ", codon_long$Sample)

# Plot with ggplot
ggplot(codon_long, aes(x = Codon, y = Value, color = Sample, group = Sample)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Codon Frequency per 1000 Codons Across Phylliidae Mitogenomes",
       x = "Codon",
       y = "Frequency per Thousand (CDspT)",
       color = "Species") +
  scale_color_discrete(
    labels = c(
      "C_oyae" = "C. oyae",
      "C_tibetense" = "C. tibetense",
      "C_westwoodii" = "C. westwoodii",
      "P_bioculatum" = "P. bioculatum",
      "P_giganteum" = "P. giganteum",
      "P_philippinicum" = "P. philippinicum"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 8),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid = element_blank()
  )
