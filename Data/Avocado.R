setwd("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit")
source("./R-Code/Sim_functions.R")

# Read data
data <- read.csv("./R-Code/DataSets/avocado.csv", sep=",")[,-c(1,2)]
data <- data[sample(nrow(data)),]

# Define categorical variables
categorical <- c(10,11, 12)

# and transform them to factors
for(c in categorical){
  data[,c] <- as.factor(data[,c])
}

# Transform factors with only two variables to binary variables
data[, which(sapply(data, nlevels) == 2)]  <- sapply(data[, which(sapply(data, nlevels) == 2)],as.numeric)-1

# Create entropy plot
entropy_cat(data, categorical[c(2,3)])
summarise_cat(data, categorical)

# Set colnames
Y <- 1
colnames(data) <- paste(replace(replace(sapply(data, FUN=class), which(sapply(data, FUN=class) %in% c('factor', 'character')), 'C'), which(!(sapply(data, FUN=class) %in% c('factor', 'character'))), 'X'), formatC(1:ncol(data), width=2, flag=0), sep = "")
colnames(data)[Y] <- 'Y'

# Split data set
strat_categ <- which.max(sapply(data, nlevels))

#Stratified split on the Variable with most levels
library(caret)
train.index <- createDataPartition(data[,strat_categ], p = 2/3, list = FALSE)
train_data <- data[ train.index,]
test_data  <- data[-train.index,]

