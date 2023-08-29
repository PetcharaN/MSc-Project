# This script will plot the expression of different alternative splicing events in triadic wheat circadian clock genes

setwd("C:/Users/asus/OneDrive - University of East Anglia/Project onedrive/Excel")

#LHY assessment

# Reformat dataset to long 

LHY_dataset <- read.csv("LHY_triad_expression.csv")
LHY_dataset <- melt(LHY_dataset, id.vars = "bin",
                           variable.name = "Time", value.name = "Expression")
LHY_dataset$Time <- as.character(LHY_dataset$Time)
LHY_dataset$Time <- as.numeric(substr(LHY_dataset$Time, 2, nchar(LHY_dataset$Time)))

#Plot LHY triad

plot <- ggplot(LHY_dataset, aes(x = as.numeric(Time), y = Expression, color = bin)) +
  geom_line() +
  labs(title = "TaLHY triad alternative splicing events",
       x = "Time",
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

#RVE assessment


RVE_dataset <- read.csv("RVE_triad_expression.csv")
RVE_dataset <- melt(RVE_dataset, id.vars = "bin",
                    variable.name = "Time", value.name = "Expression")
RVE_dataset$Time <- as.character(RVE_dataset$Time)
RVE_dataset$Time <- as.numeric(substr(RVE_dataset$Time, 2, nchar(RVE_dataset$Time)))

#Plot LHY triad

plot <- ggplot(RVE_dataset, aes(x = as.numeric(Time), y = Expression, color = bin)) +
  geom_line( data = RVE_dataset[RVE_dataset$bin == "TraesCS6A02G258600:E006",], col = "#0066CC", linetype= "solid") +
  geom_line( data = RVE_dataset[RVE_dataset$bin == "TraesCS6A02G258600:E008",], col = "#0066CC", linetype= "dotdash") +
  geom_line( data = RVE_dataset[RVE_dataset$bin == "TraesCS6A02G258600:E009",], col = "#0066CC", linetype= "dotted") +
  geom_line( data = RVE_dataset[RVE_dataset$bin == "TraesCS6A02G258600:E010",], col = "#0066CC", linetype= "dotdash") +
  geom_line( data = RVE_dataset[RVE_dataset$bin == "TraesCS6A02G258600:E016",], col = "#0066CC", linetype= "solid") +
  geom_line( data = RVE_dataset[RVE_dataset$bin == "TraesCS6B02G266800:E005",], col = "#9933CC", linetype= "solid") +
  geom_line( data = RVE_dataset[RVE_dataset$bin == "TraesCS6B02G266800:E010",], col = "#9933CC", linetype= "dotdash") +
  geom_line( data = RVE_dataset[RVE_dataset$bin == "TraesCS6D02G239800:E008",], col = "#02D09F", linetype= "solid") +
  labs(title = "TaRVE triad alternative splicing events",
       x = "Time",
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
