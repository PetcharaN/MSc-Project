# This script will plot the eigengene examples of each module provided by Rees et al


# Load Librarys
library(ggplot2)
library(reshape2)

# Load data 

setwd("C:/Users/asus/OneDrive - University of East Anglia/Project onedrive/Excel")
eigengene_dataset <- read.csv("eigengenes.csv")
eigengene_dataset <- melt(eigengene_dataset, id.vars = "Module",
                          variable.name = "Time", value.name = "Expression")
eigengene_dataset$Time <- as.character(eigengene_dataset$Time)
eigengene_dataset$Time <- as.numeric(substr(eigengene_dataset$Time, 2, nchar(eigengene_dataset$Time)))

# Plot the data
plot <- ggplot(eigengene_dataset, aes(x = as.numeric(Time), y = Expression, color = Module)) +
  geom_line( ) +
  labs(title = "Module Circadian Expression",
       x = "Time (hours)",
       y = "Normalised Gene Expression") +
  theme_classic() +
  scale_x_continuous(name ="Time", 
                     breaks=c(0,4,8,12,16,20))

vline_positions <- c(0,4,8,12,16,20)

plot <- plot +
  geom_vline(xintercept = vline_positions, 
             linetype = "dashed", 
             color = "grey")

plot
