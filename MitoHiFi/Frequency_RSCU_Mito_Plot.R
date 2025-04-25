################## Codon Frequency Plot #######################

data <- read.csv("Codon_Frequency_Mitogenome.csv")
library(ggplot2)
library(ggpattern)


# Order the data by 'Count' in descending order
data$Symbol <- factor(data$Symbol, levels = data$Symbol[order(-data$Frequency)])
library(ggplot2)
library(ggpattern)

ggplot(data, aes(x = reorder(Symbol, -Frequency), y = Frequency, pattern = Symbol)) +
  geom_bar_pattern(stat = "identity", 
                   aes(fill = ifelse(Symbol %in% c("STOP", "Ile", "Met"), "white", "grey")),  # Grey fill for others
                   color = "black",  # Black outline
                   pattern_density = 0.1, 
                   pattern_spacing = 0.02,
                   pattern_fill = "black") +  
  scale_fill_identity() +  # Ensures manual fill colors are applied
  scale_pattern_manual(values = c("STOP" = "stripe", "Ile" = "circle", "Met" = "circle", 
                                  setNames(rep("none", length(unique(data$Symbol)) - 3), 
                                           setdiff(unique(data$Symbol), c("STOP", "Ile", "Met")))), 
                       na.translate = FALSE) +
  theme_minimal() +
  labs(x = "Codon Family", y = "Frequency", title = "Codon Family Usage Within the PCGs of P. philippinicum mitogenome") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))


################ RSCU Plot #####################
library(openxlsx)
data2 <- read.xlsx("RSCU_P_philippinicum.xlsx")

ggplot(data2, aes(x = reorder(Codon, -RSCU, na.rm = TRUE), y = RSCU)) +
  geom_bar_pattern(stat = "identity", 
                   aes(fill = ifelse(Codon %in% c("ATT", "ATC", "ATA", "ATG", "TAA"), "white", "grey"),
                       pattern = ifelse(Codon %in% c("TAA", "TAG"), "stripe", 
                                        ifelse(Codon %in% c("ATT", "ATC", "ATA", "ATG"), "circle", "none"))),  
                   color = "black",  
                   pattern_density = 0.1, 
                   pattern_spacing = 0.02,
                   pattern_fill = "black") +  
  scale_fill_identity(guide = "none") +  # Removes fill legend
  scale_pattern_manual(
    values = c("stripe" = "stripe", "circle" = "circle", "none" = "none"),
    guide = "none"  # Removes pattern legend
  ) +
  theme_minimal() +
  labs(x = "Codon", y = "RSCU", title = "Codon RSCU within the PCGs of P. philippinicum mitogenome") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0.1, vjust = 0.5))


