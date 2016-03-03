setwd("~/prog/r/spatial_organization_measurement_2-03-2016")

FISHdata <- read.csv2('3DFISH.csv', header=TRUE, sep=',', dec='.')

FISHdata$X_volume <- gsub(',', '', FISHdata$X_volume)
FISHdata$X_volume <- as.numeric(FISHdata$X_volume)
FISHdata$X_surface <- gsub(',', '', FISHdata$X_surface)
FISHdata$X_surface <- as.numeric(FISHdata$X_surface)
FISHdata$Nuc_volume <- gsub(',', '', FISHdata$Nuc_volume)
FISHdata$Nuc_volume <- as.numeric(FISHdata$Nuc_volume)

FISHdata$XN <- FISHdata$X_volume/FISHdata$Nuc_volume
FISHdata <- na.omit(FISHdata)
aggregate(XN ~ Sp+Tissue, FISHdata, mean)
numFISHdata <- FISHdata
numFISHdata$Sp <- numFISHdata$Tissue <- NULL

FISHdata$Sp_Tissue <- paste(FISHdata$Sp,FISHdata$Tissue, ' ')

