library(mclust)
setwd("~/R/mclust/misc/www/images/")

###################################################

X <- iris[,1:4]
irisMclust <- Mclust(X)

png("fig01a.png", width = 480, height = 350)
par(mfrow = c(1,1), mar=c(5, 8, 3, 2))
plot(irisMclust, what = "BIC")
dev.off()

png("fig01b.png", width = 400, height = 350)
par(mfrow = c(1,1))
plot(irisMclust, what = "classification")
dev.off()

system("convert fig01a.png fig01b.png +append fig01.png")

###################################################

precipMclust <- Mclust(precip)
summary(precipMclust, parameters = TRUE)

png("fig02.png", width = 960, height = 350)
par(mfrow = c(1,2), mar=c(5, 10, 2, 2))
plot(precipMclust, what = "classification")
par(mar=c(5, 6, 2, 6))
plot(precipMclust, what = "uncertainty")
dev.off()

###################################################

densWaiting <- with(faithful, densityMclust(waiting))
faithfulDens <- densityMclust(faithful)

png("fig03.png", width = 960, height = 350)
par(mfrow = c(1,2), mar=c(5, 10, 2, 2))
plot(densWaiting, what = "density", col = "dodgerblue3", lwd = 2,
     data = faithful$waiting, breaks = 15)
par(mar=c(5, 6, 2, 6))
plot(faithfulDens, what = "density", type = "level")
points(faithful, cex = 0.5)
dev.off()

png("fig04.png", width = 960, height = 350)
par(mfrow = c(1,2), mar=c(5, 10, 2, 2))
plot(faithfulDens, what = "density", type = "contour")
par(mar=c(5, 6, 2, 6))
plot(faithfulDens, what = "density", type = "persp", col = grey(0.8))
dev.off()

###################################################

X = iris[,1:4]
Class = iris[,5]
irisMCLUSTDA <- MclustDA(X, Class, modelType = "EDDA")
predClass <- predict(irisMCLUSTDA)$classification
err <- classError(predClass, irisMCLUSTDA$class)

png("fig05.png", width = 960, height = 350)
par(mfrow = c(1,2), mar=c(5, 10, 2, 2))
plot(irisMCLUSTDA, what = "scatterplot", dim = c(1,3))
par(mar=c(5, 6, 2, 6))
clPairs(irisMCLUSTDA$data[,c(1,3)], predClass)
points(irisMCLUSTDA$data[err$misclassified,c(1,3)], pch = 4)
dev.off()

###################################################

DR = MclustDR(irisMCLUSTDA)
summary(DR)

png("fig06.png", width = 960, height = 350)
par(mfrow = c(1,2), mar=c(5, 10, 2, 2))
plot(DR, what = "contour")
par(mar=c(5, 6, 2, 6))
plot(DR, what = "boundaries", ngrid = 200)
dev.off()




