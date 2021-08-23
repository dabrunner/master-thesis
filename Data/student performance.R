setwd("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit")
source("./R-Code/Sim_functions.R")

# Read the data
data <- read.table("./R-Code/DataSets/student-por.csv",sep=";",header=TRUE)

# remove variables
data <- data[, -c(32,31)]

# define categorical variables
categorical <- c(1,2,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)

# and transform them to factors
for(c in categorical){
  data[,c] <- as.factor(data[,c])
}

summarise_cat(data, categorical)

# Transform categorical variables with only two occurences to numerical/binary variables
data[, which(sapply(data, nlevels) == 2)]  <- sapply(data[, which(sapply(data, nlevels) == 2)],as.numeric)-1

# Set column names
Y <- 31
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

