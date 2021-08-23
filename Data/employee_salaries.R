setwd("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit")
source("./R-Code/Sim_functions.R")

require(dplyr)

# Read the data
data <- read.csv("./R-Code/DataSets/employee-salaries.csv", sep=",")
data <- data[sample(nrow(data)),]


# remove obsolete variables
data <- data[, -c(1, 2, 5, 6, 7, 12, 13)]
data <- data[-which(is.na(data[,1])),]

# Define categorical variables
categorical <- c(seq(1, ncol(data), 1)[-c(2,7)])

# Change year of entry into seniority
data[,7] <- max(data[,7]) - data[,7]

# Transform categorical variables to factors
for(c in categorical){
  data[,c] <- as.factor(data[,c])
}

# Put all levels with 2 or less occurences into the level "Others"
cutoff <- 2
data <- others_factor(data, categorical, cutoff)

# Transform categorical variables with only 2 factors to binary variables
data[, which(sapply(data, nlevels) == 2)]  <- sapply(data[, which(sapply(data, nlevels) == 2)],as.numeric)-1

# set column names
Y <- 2
colnames(data) <- paste(replace(replace(sapply(data, FUN=class), which(sapply(data, FUN=class) %in% c('factor', 'character')), 'C'), which(!(sapply(data, FUN=class) %in% c('factor', 'character'))), 'X'), formatC(1:ncol(data), width=2, flag=0), sep = "")
colnames(data)[Y] <- 'Y'

# create entropy boxplot
categorical <- which(grepl("C", colnames(data)))
entropy_cat(data, categorical)

# Split data set
strat_categ <- which.max(sapply(data, nlevels))

#Stratified split on the Variable with most levels
library(caret)
train.index <- createDataPartition(data[,strat_categ], p = 2/3, list = FALSE)
train_data <- data[ train.index,]
test_data  <- data[-train.index,]


