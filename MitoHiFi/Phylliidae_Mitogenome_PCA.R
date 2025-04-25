library(openxlsx)
df <- read.xlsx("mito_stats.xlsx")
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(factoextra)

df_scaled <- scale(df[,1:6], scale = T, center = T)

#Do PCA
pca_result <- prcomp(df_scaled, center = TRUE, scale = TRUE)

#Get PCA results into a dataframe
pca_df <- as.data.frame(pca_result$x)
pca_df

#Add the "Group" column back to the df :)
pca_df$Group <- df$Group

#Plot PCA results with ggplot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "darkgreen")) +  # Customize colors
  stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed", linewidth = 1) +  # Add ellipses
  labs(title = "PCA of Length of Mitochondrial Features:Total, PCG, tRNA, rrnS, rrnL, Control Region", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        panel.border = element_blank(),legend.title = element_blank())+
  geom_text_repel(aes(label = ifelse(df$Species == "P. philippinicum", "P. philippinicum", "")),
                  size = 3, color = "black")

#Check what PC loadings are
summary(pca_result)


pca_result$rotation[, 1]
#Define PC1 loadings - probably couldve done this a better way
loadings1 <- data.frame(
  Variable = c("PCG", "rrnL", "rrnS", "tRNA", "Control_Region", "Total_Length"),
  PC1 = c(0.20552986, -0.07613023, -0.21105527, 0.25623724, -0.65267376, -0.64480755)
)


#Create loading plot for PC1
ggplot(loadings1, aes(x = Variable, y = PC1, label = Variable)) +
  geom_bar(stat = "identity", fill = "darkgray", colour = "black") +
  geom_text(aes(label = round(PC1, 3)), vjust = -0.5) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        panel.border = element_blank(),  # Ensures only axes, not full box
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)
  labs(title = "PC1 Loadings", x = "Variables", y = "Loading Value")





#PC2 Loadings
pca_result$rotation[,2]
loadings2 <- data.frame(
  Variable = c("PCG", "rrnL", "rrnS", "tRNA", "Control_Region", "Total_Length"),
  PC2 = c(-0.48770640, -0.66637362, -0.33158341, -0.444467221, -0.07629180, -0.06772958))


#Create loading plot for PC2
ggplot(loadings2, aes(x = Variable, y = PC2, label = Variable)) +
  geom_bar(stat = "identity", fill = "darkgray", color = "black") +
  geom_text(aes(label = round(PC2, 3)), vjust = -0.5) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        panel.border = element_blank(),  # Ensures only axes, not full box
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)
  labs(title = "PC2 Loadings", x = "Variables", y = "Loading Value")

#####Cumulative Proportions#####

cumulative_prop <- c(0.3553, 0.5836, 0.7851, 0.9133, 0.99775, 1.00000)
pcs <- paste("PC", 1:6)

# Create a data frame
prop_df <- data.frame(PC = pcs, Cumulative_Prop = cumulative_prop)

# Plot the bar chart
ggplot(prop_df, aes(x = PC, y = Cumulative_Prop)) +
  geom_bar(stat = "identity", fill = "darkgray", color = "black") +
  labs(title = "Cumulative Proportion of Variance by Principal Component", 
       x = "Principal Components", 
       y = "Cumulative Proportion") +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        panel.border = element_blank(),  # Ensures only axes, not full box
        axis.text.x = element_text(angle = 0, hjust = 0.5))
