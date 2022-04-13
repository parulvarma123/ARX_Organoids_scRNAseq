library(ggplot2)

#Stacked bar plot for percentage of cells for all datasets
#Proportion of cell type in each dataset was calculated using Idents function from Seurat
#Used excel sheet to arrange the CellTypes and OrganoidTypes (Ctrl, Mut) to read into the data file

#########################################################
## For hCS30DIV

#read the percent cells file for hCS30DIV
data <- read.delim(file = "~/hCS30DPercentCells.txt")
                   
#OrganoidType is either Ctrl or Mut
ggplot(data= data, aes(x=OrganoidType, y= Percentage, fill= CellTypes))+geom_col()
                   
#Preferred order of cell types
data$CellTypes <- factor(data$CellTypes, levels = c( "Radial Glia","Cycling Progenitors","Intermediate Progenitors","Deep Layer Neurons","Early Cajal-Retzius Neurons", "Immature Interneurons","Astrocytes", "Other"))
                   
#Make the plot with custom colors matching the colors of the UMAP
hCS30DCellPer <- ggplot(data = data, aes(x=OrganoidType, y= Percentage,fill= factor(CellTypes)), bg= "white")+geom_col() +scale_fill_manual(values = c("hotpink1","deepskyblue", "violetred4","gold1", "orangered", "olivedrab3", "purple3","gray"))
                   
hCS30DCellPer+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(),axis.line = element_line(colour = "black"))
                   
#Save the plot
ggsave("~/hCS30DCellPercent.jpeg", width = 4, height = 4, units = c("in"), dpi = 300, bg= "white")
dev.off()

##########################################################
##For hCS120DIV

#Read the percent cells file for hCS120DIV
data1 <- read.delim(file = "~/hCS120DPercentCells.txt")
                    
#OrganoidType is either Ctrl or Mut
ggplot(data= data1, aes(x=OrganoidType, y= Percentage,fill= CellTypes))+geom_col()
                    
#Preferred order of cell types
data1$CellTypes <- factor(data1$CellTypes, levels = c( "Radial Glia","oRG", "Cycling Progenitors","Intermediate Progenitors","Deep Layer Neurons","Upper Layer Neurons", "Interneuron Progenitor", "Interneurons", "Other"))
                    
#Make the plot with custom colors matching the colors of the UMAP
hCS120DCellPer <- ggplot(data = data1, aes(x=OrganoidType, y= Percentage,fill= factor(CellTypes)), bg= "white")+geom_col()+scale_fill_manual(values = c( "hotpink1","navyblue","deepskyblue","violetred4", "gold1","orangered", "chartreuse", "green4","gray"))
                    
hCS120DCellPer+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(),axis.line = element_line(colour = "black"))
                    
#Save the plot
ggsave("~/hCS120DCellPercentStackedBarchart.jpeg", width = 4, height = 4, units = c("in"), dpi = 300, bg= "white")
dev.off()

###########################################################
##For hSS30DIV

#Read the percent cells file for hSS30DIV
data2 <- read.delim(file = "~/hSS30DPercentCells.txt")
                    
#OrganoidType is either Ctrl or Mut
ggplot(data= data2, aes(x=OrganoidType, y= Percentage,fill= CellTypes))+geom_col()
                    
#Preferred order of cell types
data2$CellTypes <- factor(data2$CellTypes, levels = c( "Radial Glia", "Cycling Progenitors", "Ventral Intermediate Progenitor","Immature Interneurons", "Astrocytes", "Other","Unknown"))
                    
#Make the plot with custom colors matching the colors of the UMAP
hSS30DCellPer <- ggplot(data = data2, aes(x=OrganoidType, y= Percentage,fill= factor(CellTypes)), bg= "white")+geom_col() +scale_fill_manual(values = c("hotpink1", "deepskyblue", "greenyellow","olivedrab3","purple3", "gray28", "gray"))
                    
hSS30DCellPer+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(),axis.line = element_line(colour = "black"))
                    
#Save the plot
ggsave("~/hSS30DCellPercent.jpeg", width = 4, height = 4, units = c("in"), dpi = 300, bg= "white")
dev.off()

###########################################################
##For hSS120DIV

#Read the percent cells file for hSS120DIV
data3 <- read.delim(file = "~/hSS120DPercentCells.txt")
                    
#OrganoidType is either Ctrl or Mut
ggplot(data= data3, aes(x=OrganoidType, y= Percentage,fill= CellTypes))+geom_col()
                    
#Preferred order of cell types
data3$CellTypes <- factor(data3$CellTypes, levels = c( "Radial Glial", "Cycling Progenitors","Ventral Intermediate Progenitor", "Immature Interneurons", "SST Interneurons", "Other","Unknown"))
                    
#Make the plot with custom colors matching the colors of the UMAP
hSS120DCellPer <- ggplot(data = data3, aes(x=OrganoidType, y= Percentage,fill= factor(CellTypes)), bg= "white")+geom_col() +scale_fill_manual(values = c("hotpink", "deepskyblue", "darkred","olivedrab3", "green4", "gray28", "gray"))
                    
hSS120DCellPer+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(),axis.line = element_line(colour = "black"))
                    
#Save the plot
ggsave("~/hSS120DCellPercent.jpeg", width = 4, height = 4, units = c("in"), dpi = 300, bg= "white")
dev.off()