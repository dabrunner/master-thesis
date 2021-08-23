setwd("/Users/danielbrunner/Documents/Ausbildung/HSG/Masterstudium/Masterarbeit")
# 
# n <- 150 #1000
# p <- 3
# k <- 2
# ngl <- 2 # 100, sample(c(50, 100), 1)
# pl <- .9
# type <- "latent" # "nonlinear"#"latent"#"hermite"#"latent"#"global"
# 
# data <- create_data(n, p, k, ngl = ngl, pl = pl, type = type)
# train_test <- createDataPartition(factor(data$g), p = 0.8)
# 
# train <- list(
#   x = data$x[train_test$Resample1, ],
#   g = data$g[train_test$Resample1],
#   y = data$y[train_test$Resample1]
# )
#  l <- data$l[train_test$Resample1]

#saveRDS(train, "./data_figure2.rds")
train <- readRDS("./data_figure2.rds")
 
X <- train$x
G <- train$g
G2 <- as.numeric(G)
G2 <- as.data.frame(G2)
G2[G2==1] <- "A"
G2[G2==2] <- "C"
G2[G2==3] <- "B"
G2[G2==4] <- "D"
G <- as.factor(G2$G2)

X <- X[,-3]


# means encoding
p <- dim(X)[2]
CM <- as.matrix(aggregate(X, list(G), mean)[, 2:(p + 1)])

plot(X[, 1], X[,2], xlim=c(-4,2), ylim=c(-3,3))
plot(CM[,1], CM[,2],  xlim=c(-4,2), ylim=c(-3,3))


# data for gg plot
df1 <- data.frame(X, G)
colnames(df1) <- c("X1", "X2", "Category")

df2 <- data.frame(CM, factor(c("A", "B", "C", "D")))
colnames(df2) <-  c("X1", "X2", "Category")

plot1 <- ggplot(df1, aes(x=X1, y=X2, shape = Category, colour = Category, fill = Category))+
  geom_point(size = 4)+
  scale_shape_manual(values=c(21,21, 24, 24)) +
  scale_fill_manual(values =alpha(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), .5)) +
  xlim(-4,2.5)+ylim(-3.53,3.5) +
  scale_color_manual(values=c("black", "black", "black", "black")) +
  labs(x=expression(X[1]), y=expression(X[2]), colour = "Category") +
  ggtitle(expression(paste("(", X[1], ", " , X[2], ")", sep = ""))) +
  theme_minimal(base_family = "Palatino") +
  theme(legend.position = "none", plot.title=element_text(hjust=0.5))

plot2 <- ggplot(df2, aes(x=X1, y=X2, shape = Category, colour = Category, fill = Category))+
  geom_point(size = 4)+
  scale_shape_manual(values=c(21,21, 24, 24)) +
  scale_fill_manual(values =alpha(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), .5)) +
  scale_color_manual(values=c("black", "black", "black", "black")) +
  xlim(-4,2.5)+ylim(-3.53,3.5) +
  labs(x=expression(X[1]), y=expression(X[2]), colour = "Category") +
  ggtitle(expression(paste("(", "E[", X[1], '|', G, "], E[" , X[2], '|', G, "])", sep = ""))) +
  theme_minimal(base_family = "Palatino") +
  theme(plot.title=element_text(hjust=0.5))


require(gridExtra)
require(cowplot)
plot_grid(plot1, plot2, align = "h", nrow = 1, rel_widths = c(0.44, 0.56))


ggsave("plot_MA.jpg", units = c("cm"), height = 1.1*7.44, width = 1.1*17.8, dpi = 150)



