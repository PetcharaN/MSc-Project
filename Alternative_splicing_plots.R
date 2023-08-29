# This script will plot the expression of different alternative splicing events in wheat circadian clock genes


# Load Librarys
library(ggplot2)
library(reshape2)

# Load data 
# Save dataset as csv first
setwd("C:/Users/asus/OneDrive - University of East Anglia/Project onedrive/Excel")
expression_dataset <- read.csv("splicing_events.csv")

# Reformat dataset to long 
expression_dataset <- melt(expression_dataset, id.vars = "Bin",
                           variable.name = "Time", value.name = "Expression")
expression_dataset$Time <- as.character(expression_dataset$Time)
expression_dataset$Time <- as.numeric(substr(expression_dataset$Time, 2, nchar(expression_dataset$Time)))

# MYB gene splicing events 
MYB <- c("TraesCS7A02G299400:E006",
         "TraesCS7B02G188000:E006",
         "TraesCS7D02G295400:E007",
         "TraesCS7A02G553800:E004",
         "TraesCS7B02G478200:E002",
         "TraesCS6A02G258600:E006",
         "TraesCS6A02G258600:E008",
         "TraesCS6A02G258600:E009",
         "TraesCS6A02G258600:E010",
         "TraesCS6A02G258600:E016",
         "TraesCS6B02G266800:E005",
         "TraesCS6B02G266800:E010",
         "TraesCS6D02G239800:E008",
         "TraesCS4A02G474100:E005",
         "TraesCS4A02G474100:E011",
         "TraesCS7D02G014900:E008",
         "TraesCS7A02G470700:E008",
         "TraesCS7D02G458000:E009")

MYB_filtered_data <- expression_dataset[expression_dataset$Bin%in%MYB,]

# Plot the data
plot <- ggplot(MYB_filtered_data, aes(x = as.numeric(Time), y = Expression, color = Bin)) +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS7A02G299400:E006",], col = "#02D09F", linetype= "longdash") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS7B02G188000:E006",], col = "#02D09F", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS7D02G295400:E007",], col = "#02D09F", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS7A02G553800:E004",], col = "#02D09F", linetype= "longdash") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS7B02G478200:E002",], col = "#02D09F", linetype= "dotted") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS6A02G258600:E006",], col = "#0066CC", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS6A02G258600:E008",], col = "#0066CC", linetype= "dotdash") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS6A02G258600:E009",], col = "#0066CC", linetype= "dotted") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS6A02G258600:E010",], col = "#0066CC", linetype= "dotdash") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS6A02G258600:E016",], col = "#0066CC", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS6B02G266800:E005",], col = "#006600", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS6B02G266800:E010",], col = "#006600", linetype= "dotdash") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS6D02G239800:E008",], col = "#02D09F", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS4A02G474100:E005",], col = "#9933CC", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS4A02G474100:E011",], col = "#9933CC", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS7D02G014900:E008",], col = "#02D09F", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS7A02G470700:E008",], col = "#02D09F", linetype= "solid") +
  geom_line( data = MYB_filtered_data[MYB_filtered_data$Bin == "TraesCS7D02G458000:E009",], col = "#02D09F", linetype= "solid") +
   
  labs(title = "MYB Circadian Splicing Activity",
       x = "Time (hours)",
       y = "Relative Alternative Splicing Isoform Expression") +
  theme_classic() +
scale_x_continuous(name ="Time", 
                   breaks=c(0,4,8,12,16,20))


vline_positions <- c(0,4,8,12,16,20)
plot <- plot +
  geom_vline(xintercept = vline_positions, 
             linetype = "dashed", 
             color = "grey")

plot


#LNK1 plot

# LNK1 gene splicing events 
LNK1 <- c("TraesCS4B02G153700:E005",
         "TraesCS4B02G153700:E010",
         "TraesCS4B02G153700:E013")

LNK1_filtered_data <- expression_dataset[expression_dataset$Bin%in%LNK1,]

# Plot the data
plot <- ggplot(LNK1_filtered_data, aes(x = as.numeric(Time), y = Expression, color = Bin)) +
  geom_line( data = LNK1_filtered_data[LNK1_filtered_data$Bin == "TraesCS4B02G153700:E005",], col = "#C65911", linetype= "solid") +
  geom_line( data = LNK1_filtered_data[LNK1_filtered_data$Bin == "TraesCS4B02G153700:E010",], col = "#C65911", linetype= "longdash") +
  geom_line( data = LNK1_filtered_data[LNK1_filtered_data$Bin == "TraesCS4B02G153700:E013",], col = "#C65911", linetype= "longdash") +

  
  labs(title = "LNK1 Circadian Splicing Activity",
       x = "Time (hours)",
       y = "Relative Alternative Splicing Isoform Expression") +
  theme_classic() +
  scale_x_continuous(name ="Time", 
                     breaks=c(0,4,8,12,16,20))


vline_positions <- c(0,4,8,12,16,20)
plot <- plot +
  geom_vline(xintercept = vline_positions, 
             linetype = "dashed", 
             color = "grey")

plot

#PRR plot

# PRR gene splicing events 
PRR <- c("TraesCS2D02G079600:E007",
          "TraesCS4D02G199600:E006")

PRR_filtered_data <- expression_dataset[expression_dataset$Bin%in%PRR,]

# Plot the data
plot <- ggplot(PRR_filtered_data, aes(x = as.numeric(Time), y = Expression, color = Bin)) +
  geom_line( data = PRR_filtered_data[PRR_filtered_data$Bin == "TraesCS2D02G079600:E007",], col = "#CE0259", linetype= "solid") +
  geom_line( data = PRR_filtered_data[PRR_filtered_data$Bin == "TraesCS4D02G199600:E006",], col = "#D774E8", linetype= "solid") +
 
  
  labs(title = "PRR Circadian Splicing Activity",
       x = "Time (hours)",
       y = "Relative Alternative Splicing Isoform Expression") +
  theme_classic() +
  scale_x_continuous(name ="Time", 
                     breaks=c(0,4,8,12,16,20))


vline_positions <- c(0,4,8,12,16,20)
plot <- plot +
  geom_vline(xintercept = vline_positions, 
             linetype = "dashed", 
             color = "grey")

plot



