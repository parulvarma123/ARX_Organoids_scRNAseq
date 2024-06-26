library(ggplot2)

#Stacked bar plot for percentage of cells for all datasets
#Proportion of cell type in each dataset was calculated using Idents function from Seurat
#Used excel sheet to arrange the CellTypes and OrganoidTypes (Ctrl, Mut) to read into the data file

#########################################################
## For hCO30DIV

#read the percent cells file for hCO30DIV
data <- read.delim(file = "~/hCS30DPercentCells.txt")
                   
#OrganoidType is either Ctrl or Mut
ggplot(data= data, aes(x=OrganoidType, y= Percentage, fill= CellTypes))+geom_col()
                   
#Preferred order of cell types
data$CellTypes <- factor(data$CellTypes, levels = c("Cajal-Retzius Neurons","Deep Layer Neurons","Intermediate Progenitors","Cycling Progenitors", "Radial Glia"))
ggplot(data= data, aes(x=OrganoidType, y= Percentage, fill= CellTypes))+geom_col()

#Make the plot with custom colors matching the colors of the UMAP
hCS30DCellPer <- ggplot(data = data, aes(x=OrganoidType, y= Percentage,fill= factor(CellTypes)), bg= "white")+geom_col() +scale_fill_manual(values = c("orangered", "gold1",  "violetred4", "deepskyblue", "hotpink1"))

hCS30DCellPer+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(),axis.line = element_line(colour = "black"))
                   
#Save the plot
ggsave("~/hCS30DCellPercent.jpeg", width = 4, height = 4, units = c("in"), dpi = 300, bg= "white")
dev.off()

##########################################################
##For hCO120DIV

#Read the percent cells file for hCO120DIV
data1 <- read.delim(file = "~/hCS120DPercentCells.txt")
                    
#OrganoidType is either Ctrl or Mut
ggplot(data= data1, aes(x=OrganoidType, y= Percentage,fill= CellTypes))+geom_col()
                    
#Preferred order of cell types
data2$CellTypes <- factor(data2$CellTypes, levels = c("Upper Layer Neurons","Deep Layer Neurons","Intermediate Progenitors", "Cycling Progenitors", "Outer Radial Glia", "Radial Glia"))

#now check the plot
ggplot(data= data2, aes(x=OrganoidType, y= Percentage,fill= CellTypes))+geom_col()


#Make the plot with custom colors matching the colors of the UMAP
hCS120DCellPer <- ggplot(data = data2, aes(x=OrganoidType, y= Percentage,fill= factor(CellTypes)), bg= "white")+geom_col()+scale_fill_manual(values = c("orangered", "gold1","violetred4","deepskyblue","navyblue", "hotpink1"))

hCS120DCellPer
hCS120DCellPer+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(),axis.line = element_line(colour = "black"))
                    
#Save the plot
ggsave("~/hCS120DCellPercentStackedBarchart.jpeg", width = 4, height = 4, units = c("in"), dpi = 300, bg= "white")
dev.off()

###########################################################
##For hGEO30DIV

#Read the percent cells file for hGEO30DIV
data2 <- read.delim(file = "~/hSS30DPercentCells.txt")
                    
#OrganoidType is either Ctrl or Mut
ggplot(data= data2, aes(x=OrganoidType, y= Percentage,fill= CellTypes))+geom_col()
                    
#Preferred order of cell types
data4$CellTypes <- factor(data4$CellTypes, levels = c("Inhibitory Neurons", "Ventral Intermediate Progenitors","Cycling Progenitors", "Radial Glial Cells" ))

#Make the plot with custom colors matching the colors of the UMAP
hSS30DCellPer <- ggplot(data = data4, aes(x=OrganoidType, y= Percentage,fill= factor(CellTypes)), bg= "white")+geom_col() +scale_fill_manual(values = c("olivedrab3", "violetred4", "deepskyblue", "hotpink1"))

hSS30DCellPer+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(),axis.line = element_line(colour = "black"))
                    
#Save the plot
ggsave("~/hSS30DCellPercent.jpeg", width = 4, height = 4, units = c("in"), dpi = 300, bg= "white")
dev.off()

###########################################################
##For hGEO120DIV

#Read the percent cells file for hGEO120DIV
data3 <- read.delim(file = "~/hSS120DPercentCells.txt")
                    
#OrganoidType is either Ctrl or Mut
ggplot(data= data3, aes(x=OrganoidType, y= Percentage,fill= CellTypes))+geom_col()
                    
#Preferred order of cell types
data6$CellTypes <- factor(data6$CellTypes, levels = c("Inhibitory Neurons", "Ventral Intermediate Progenitors","Cycling Progenitors", "Radial Glial Cells" ))

#Make the plot with custom colors matching the colors of the UMAP
hSS120DCellPer <- ggplot(data = data6, aes(x=OrganoidType, y= Percentage,fill= factor(CellTypes)), bg= "white")+geom_col() +scale_fill_manual(values = c("olivedrab3", "violetred4", "deepskyblue", "hotpink1" ))

hSS120DCellPer+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(),axis.line = element_line(colour = "black"))
                  
#Save the plot
ggsave("~/hSS120DCellPercent.jpeg", width = 4, height = 4, units = c("in"), dpi = 300, bg= "white")
dev.off()
