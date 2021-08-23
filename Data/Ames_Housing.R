setwd("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit")
source("./R-Code/Sim_functions.R")

library(dplyr)

data <- read.csv("./R-Code/DataSets/Ames_Housing.csv", sep=",")[,-1]
data <- data[sample(nrow(data)),]

# eliminate single records
data[,19] <- bin(data[,19], 1868, 1908, 10, bin_min = T)

# Indicate categorical variables
cat <- c(1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 28, 29, 30, 31, 32, 33, 35, 39, 40, 41,42, 53, 55, 57, 58, 59, 62, 63, 64, 71, 72, 73, 75, 76, 77, 78)

# and  transform them to factors
for (i in cat){
  data[,i] <- as.factor(data[,i])
}

# Data Cleaning
data[,19] <- recode_factor(data[,19], "1917"="1918", "1913"="1912", "1911"="1912", "2010"="2009")
data[,1] <- recode_factor(data[,1], "One_and_Half_Story_PUD_All_Ages"="One_and_Half_Story_Unfinished_All_Ages")
data <- cut_factor(data,2,2)
data <- cut_factor(data,9,2)
data <- cut_factor(data, 12,2)
data <- cut_factor(data, 14,2)
data <- cut_factor(data, 22,1)
data <- cut_factor(data,23,1)
data <- cut_factor(data,24,1)
data <- cut_factor(data,25,25)
data <- cut_factor(data,39,2)
data <- cut_factor(data,42,1)
data <- cut_factor(data, 53, 70)
data <- cut_factor(data,73,4)
data[,77] <- recode_factor(data[,77], "VWD"="Oth")

summarise_cat(data, cat) 

# Transform categorical variables with only two levels to binary variables
data[, which(sapply(data, nlevels) == 2)]  <- sapply(data[, which(sapply(data, nlevels) == 2)],as.numeric)-1

# Set column names
Y <- ncol(data)-2
colnames(data) <- paste(replace(replace(sapply(data, FUN=class), which(sapply(data, FUN=class) %in% c('factor', 'character')), 'C'), which(!(sapply(data, FUN=class) %in% c('factor', 'character'))), 'X'), formatC(1:ncol(data), width=2, flag=0), sep = "")
colnames(data)[Y] <- 'Y'

# Create entropy boxplot 
cat <- which(grepl("C", colnames(data)))
entropy_cat(data, cat)

# Split data set
strat_categ <- which.max(sapply(data, nlevels))

#Stratified split on the Variable with most levels
library(caret)
train.index <- createDataPartition(data[,strat_categ], p = 2/3, list = FALSE)
train_data <- data[ train.index,]
test_data  <- data[-train.index,]

