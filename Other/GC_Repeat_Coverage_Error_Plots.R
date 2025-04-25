################## Error Rate #######################

data <- read.delim("HiFiASM_merqury.HiFiASM.qv", header = FALSE) #Do both before and after purge_dups
head(data)
library(ggplot2)
data$X <- 1:nrow(data)

ggplot(data, aes(x = X, y = V5)) +
  geom_bar(stat = "identity", fill = "gray", width =1) +  # Bar chart
  theme_minimal() +
  labs(title = "Error rate across contigs of HiFiASM assembly of Phyllium philippinicum genome",
       x = "Contig Index",
       y = "Error Rate") +
  scale_x_continuous(breaks = seq(0, nrow(data), by = 50)) +  # Label every 1000th contig
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




############## GC Content ###################

gcdata <- read.csv("gc_content_per_contig.csv", header = FALSE) #Name Changed
gcdata$X <- 1:nrow(gcdata)

ggplot(gcdata, aes(x = X, y = V2)) +
  geom_bar(stat = "identity", fill = "gray", width = 1) +  # Bar chart
  theme_minimal() +
  labs(title = "GC content across contigs of HiFiASM assembly of Phyllium philippinicum genome",
       x = "Contig Index",
       y = "GC (%)") +
  scale_x_continuous(breaks = seq(0, nrow(data), by = 50)) +  # Label every 1000th contig
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylim(0,100)




######################## Coverage #############################

cvdata <- read.csv("coverage_per_contig_SORTED.csv", header = FALSE) #Name Changed
cvdata$X <- 1:nrow(cvdata)

ggplot(cvdata, aes(x = X, y = V2)) +
  geom_bar(stat = "identity", fill = "gray", width = 1) +  # Bar chart
  theme_minimal() +
  labs(title = "Coverage across contigs of HiFiASM assembly of Phyllium philippinicum genome",
       x = "Contig Index",
       y = "Coverage") +
  scale_x_continuous(breaks = seq(0, nrow(data), by = 50)) +  # Label every 1000th contig
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
which.max(cvdata$V2)



################# Repeat Contig ##########################

rpdata <- read.csv("repeat_percentage.csv", header = FALSE) # Name Changed
rpdata$X <- 1:nrow(rpdata)

ggplot(rpdata, aes(x = X, y = V2)) +
  geom_bar(stat = "identity", fill = "gray", width = 1) +  # Bar chart
  theme_minimal() +
  labs(title = "Repeat content of contigs across HiFiASM assembly of Phyllium philippinicum genome",
       x = "Contig Index",
       y = "Repeat Content (%)") +
  scale_x_continuous(breaks = seq(0, nrow(rpdata), by = 50)) +  # Label every 1000th contig
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
max_value <- which.max(cvdata$V2)
print(data[454,])
##### Find all contigs with a error rate of 0
result <- rpdata[cvdata[, 2] > 99, 1]  # Filter rows where column 5 is 0 and select column 1
print(result)

