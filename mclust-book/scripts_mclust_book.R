# Chapter 1 ###################################################################

## install.packages("mclust")
library("mclust")


# Chapter 2 ###################################################################

## ----------------------------------------------------------------------------
data(Snapper, package = "FSAdata")
x <- Snapper[,1]

mod1 <- densityMclust(x, G = 1, plot = FALSE)
mod2 <- densityMclust(x, G = 2, plot = FALSE)
mod3 <- densityMclust(x, G = 3, plot = FALSE)
mod4 <- densityMclust(x, G = 4, plot = FALSE)

par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 2, 0, 0))
x0 <- extendrange(x, f = 0.1)
x0 <- seq(x0[1], x0[2], length = 1000)
#
cdens <- predict(mod1, newdata = x0, what = "cdens")
cdens <- sweep(cdens, 2, mod1$parameters$pro, "*")
plot(mod1, x, what = "density", lwd = 2, breaks = 20, ylim = c(0, 0.4))
matplot(x0, cdens, type = "l", lty = 2, col = 1, add = TRUE)
#
cdens <- predict(mod2, newdata = x0, what = "cdens")
cdens <- sweep(cdens, 2, mod2$parameters$pro, "*")
plot(mod2, x, what = "density", lwd = 2, breaks = 20, ylim = c(0, 0.4))
matplot(x0, cdens, type = "l", lty = 2, col = 1, add = TRUE)
#
cdens <- predict(mod3, newdata = x0, what = "cdens")
cdens <- sweep(cdens, 2, mod3$parameters$pro, "*")
plot(mod3, x, what = "density", lwd = 2, breaks = 20, ylim = c(0, 0.4))
matplot(x0, cdens, type = "l", lty = 2, col = 1, add = TRUE)
#
cdens <- predict(mod4, newdata = x0, what = "cdens")
cdens <- sweep(cdens, 2, mod4$parameters$pro, "*")
plot(mod4, x, what = "density", lwd = 2, breaks = 20, ylim = c(0, 0.4))
matplot(x0, cdens, type = "l", lty = 2, col = 1, add = TRUE)
#
mtext("Fish length (in)", side = 1, line = 1, outer = TRUE)
mtext("Density", side = 2, line = 1, outer = TRUE)


# Chapter 3 ###################################################################

## ----------------------------------------------------------------------------
data("diabetes", package = "rrcov")
X <- diabetes[, 1:5]
Class <- diabetes$group
table(Class)

## ----------------------------------------------------------------------------
clp <- clPairs(X, Class, lower.panel = NULL)
clPairsLegend(0.1, 0.3, class = clp$class, col = clp$col, pch = clp$pch)

## ----------------------------------------------------------------------------
mod <- Mclust(X, G = 3, modelNames = "VVV")

## ----------------------------------------------------------------------------
summary(mod)

## ----------------------------------------------------------------------------
summary(mod, parameters = TRUE)

## ----------------------------------------------------------------------------
table(Class, Cluster = mod$classification)
adjustedRandIndex(Class, mod$classification)

## ----------------------------------------------------------------------------
plot(mod, what = "classification")

## ----------------------------------------------------------------------------
plot(mod, what = "classification", fillEllipses = TRUE)

## ----------------------------------------------------------------------------
plot(mod, what = "classification", dimens = c(3, 4), fillEllipses = TRUE)

## ----------------------------------------------------------------------------
plot(mod, dimens = c(3, 4), what = "uncertainty")

## ----------------------------------------------------------------------------
data("thyroid", package = "mclust")
X <- data.matrix(thyroid[, 2:6])
Class <- thyroid$Diagnosis
clp <- clPairs(X, Class, lower.panel = NULL,
               symbols = c(0, 1, 2), 
               colors = c("gray50", "black", "red3")) 
clPairsLegend(0.1, 0.3, title = "Thyroid diagnosis:", class = clp$class, 
              col = clp$col, pch = clp$pch)

## ----------------------------------------------------------------------------
mod <- Mclust(X)

## ----------------------------------------------------------------------------
mod$BIC

## ----------------------------------------------------------------------------
summary(mod$BIC, k = 5)

## ----------------------------------------------------------------------------
plot(mod, what = "BIC", 
     legendArgs = list("bottomright", ncol = 5))

## ----------------------------------------------------------------------------
summary(mod, parameters = TRUE)

## ----------------------------------------------------------------------------
plot(mod, what = "classification")

## ----------------------------------------------------------------------------
table(Class, Cluster = mod$classification)
adjustedRandIndex(Class, mod$classification)

## ----------------------------------------------------------------------------
z  <- mod$z               # posterior conditional probabilities
cl <- mod$classification  # MAP clustering
G  <- mod$G               # number of clusters
sclass <- 10 # class separation
sedge <- 3   # edge spacing
L <- nrow(z) + G*(sclass+2*sedge)
plot(1:L, runif(L), ylim = c(0, 1), type = "n", axes = FALSE, 
     ylab = "Posterior conditional probabilities", xlab = "")
axis(2)
col <- mclust.options("classPlotColors")
l <- sclass
for (k in 1:G)
{
 i <- which(cl == k)
 ord <- i[order(z[i, k], decreasing = TRUE)]
 for (j in 1:G)
    points((l+sedge)+1:length(i), z[ord, j], 
           pch = as.character(j), col = col[j])
 rect(l, 0, l+2*sedge+length(i), 1, 
      border = col[k], col = col[k], lwd = 2, density = 0)
 l <- l + 2*sedge + length(i) + sclass
}

## ----------------------------------------------------------------------------
data("wine", package = "gclus")
Class <- factor(wine$Class, levels = 1:3,
                labels = c("Barolo", "Grignolino", "Barbera"))
X <- data.matrix(wine[, -1])

## ----------------------------------------------------------------------------
mod <- Mclust(X)
summary(mod$BIC, k = 3)

## ----------------------------------------------------------------------------
plot(mod, what = "BIC", 
     ylim = range(mod$BIC[, -(1:2)], na.rm = TRUE),
     legendArgs = list(x = "bottomleft"))

## ----------------------------------------------------------------------------
summary(mod)
table(Class, mod$classification)
adjustedRandIndex(Class, mod$classification)

## ----------------------------------------------------------------------------
norm01 <- function(x) (x - min(x))/(max(x) - min(x))
M <- apply(t(mod$parameters$mean), 2, norm01)
heatmap(M, Rowv = NA, scale = "none", margins = c(8, 2),
        labRow = paste("Cluster", 1:mod$G), cexRow = 1.2)

## ----------------------------------------------------------------------------
plot(mod, what = "classification", dimens = c(1, 2, 7))

## ----------------------------------------------------------------------------
data("faithful", package = "datasets")
plot(faithful)

## ----------------------------------------------------------------------------
BIC <- mclustBIC(faithful)
BIC
plot(BIC)

## ----------------------------------------------------------------------------
mod1 <- Mclust(faithful, x = BIC)

## ----------------------------------------------------------------------------
ICL <- mclustICL(faithful)
ICL

## ----------------------------------------------------------------------------
plot(ICL)

## ----------------------------------------------------------------------------
mod2 <- Mclust(faithful, G = 2, modelNames = "VVE")

## ----------------------------------------------------------------------------
plot(mod1, what = "classification", fillEllipses = TRUE)
plot(mod2, what = "classification", fillEllipses = TRUE)

## ----------------------------------------------------------------------------
LRT <- mclustBootstrapLRT(faithful, modelName = "VVV")
LRT

## ----------------------------------------------------------------------------
plot(LRT, G = 1)
plot(LRT, G = 2)

## ----------------------------------------------------------------------------
data("hemophilia", package = "rrcov")
X <- hemophilia[, 1:2]
Class <- as.factor(hemophilia$gr) 
clp <- clPairs(X, Class, symbols = c(16, 0), colors = "black")
clPairsLegend(0.8, 0.2, class = clp$class, col = clp$col, pch = clp$pch)

## ----------------------------------------------------------------------------
mod <- Mclust(X, G = 2, modelName = "VVV")
summary(mod, parameters = TRUE)

## ----------------------------------------------------------------------------
plot(mod, what = "classification", fillEllipses = TRUE)

## ----------------------------------------------------------------------------
clp <- clPairs(X, Class, symbols = c(16, 0), colors = "black")
clPairsLegend(0.8, 0.2, class = clp$class, col = clp$col, pch = clp$pch)
plot(mod, what = "classification", fillEllipses = TRUE)

## ----------------------------------------------------------------------------
boot <- MclustBootstrap(mod, nboot = 999, type = "bs")

## ----------------------------------------------------------------------------
summary(boot, what = "se")

## ----------------------------------------------------------------------------
summary(boot, what = "ci", conf.level = 0.9)

## ----------------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(boot, what = "pro")

## ----------------------------------------------------------------------------
par(mfcol = c(2, 2))
plot(boot, what = "mean")

## ----------------------------------------------------------------------------
wlboot <- MclustBootstrap(mod, nboot = 999, type = "wlbs")
summary(wlboot, what = "se")

## ----------------------------------------------------------------------------
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
boot.ci <- summary(boot, what = "ci")
wlboot.ci <- summary(wlboot, what = "ci")
for (j in 1:mod$d)
{ 
  plot(1:mod$G, mod$parameters$mean[j, ], col = 1:mod$G, pch = 15,
       ylab = colnames(X)[j], xlab = "Mixture component",
       ylim = range(boot.ci$mean, wlboot.ci$mean), 
       xlim = c(.5, mod$G+.5), xaxt = "n")
  points(1:mod$G+0.2, mod$parameters$mean[j, ], col = 1:mod$G, pch = 15)
  axis(side = 1, at = 1:mod$G)
  with(boot.ci, 
       errorBars(1:G, mean[1, j, ], mean[2, j, ], col = 1:G))
  with(wlboot.ci, 
       errorBars(1:G+0.2, mean[1, j, ], mean[2, j, ], col = 1:G, lty = 2))
}

## ----------------------------------------------------------------------------
data("precip", package = "datasets")
names(precip)[c(45,52)] = c("Bismarck", "Pittsburgh")
dotchart(sort(precip), cex = 0.6, pch = 19,
         xlab = "Average annual rainfall (in inches)")

## ----------------------------------------------------------------------------
mod <- Mclust(precip)
summary(mod, parameters = TRUE)
plot(mod, what = "BIC", legendArgs = list(x = "bottomleft"))

## ----------------------------------------------------------------------------
par(mfrow=c(1, 2))
plot(mod, what = "classification")
plot(mod, what = "uncertainty")

## ----------------------------------------------------------------------------
x <- data.frame(precip, clusters = mod$classification)
rownames(x) <- make.unique(names(precip)) # correct duplicated names 
x <- x[order(x$precip), ]
dotchart(x$precip, labels = rownames(x),
         groups = factor(x$clusters, levels = 2:1, 
                         labels = c("Cluster 2", "Cluster 1")),
         cex = 0.6, pch = 19, 
         color = mclust.options("classPlotColors")[x$clusters],
         xlab = "Average annual rainfall (in inches)")

## ----------------------------------------------------------------------------
data("EuroUnemployment", package = "mclust")
summary(EuroUnemployment)

## ----------------------------------------------------------------------------
HC_EII <- hc(EuroUnemployment, modelName = "EII")
HC_VVV <- hc(EuroUnemployment, modelName = "VVV")

## ----------------------------------------------------------------------------
plot(HC_EII, what = "merge", labels = TRUE, hang = 0.02)
plot(HC_VVV, what = "merge", labels = TRUE, hang = 0.02)

## ----------------------------------------------------------------------------
plot(HC_EII, what = "loglik")
plot(HC_VVV, what = "loglik")

## ----------------------------------------------------------------------------
data("HRstars", package = "GDAdata")
set.seed(0)
initial <- kmeans(HRstars[, -1], centers = 100, nstart=10)$cluster
HC_VVV <- hc(HRstars[, -1], modelName = "VVV", 
             partition = initial, use = "VARS")
HC_VVV

## ----------------------------------------------------------------------------
data("flea", package = "tourr")
X <- data.matrix(flea[, 1:6])
Class <- factor(flea$species, 
                labels = c("Concinna", "Heikertingeri", "Heptapotamica")) 
table(Class)

## ----------------------------------------------------------------------------
col <- mclust.options("classPlotColors")
clp <- clPairs(X, Class, lower.panel = NULL, gap = 0,
               symbols = c(16, 15, 17), 
               colors = adjustcolor(col, alpha.f = 0.5))
clPairsLegend(x = 0.1, y = 0.3, class = clp$class, 
              col = col, pch = clp$pch,
              title = "Flea beatle species")

## ----------------------------------------------------------------------------
# set the default for the current session
mclust.options("hcUse" = "VARS")
mod1 <- Mclust(X)
# or specify the initialization method only for this model
# mod1 <- Mclust(X, initialization = list(hcPairs = hc(X, use = "VARS")))
summary(mod1)
table(Class, mod1$classification)
adjustedRandIndex(Class, mod1$classification)

## ----------------------------------------------------------------------------
mod2 <- Mclust(X[, 6:1])
summary(mod2)
table(Class, mod2$classification)
adjustedRandIndex(Class, mod2$classification)

## ----------------------------------------------------------------------------
mod3 <- Mclust(X, initialization = list(hcPairs = hcRandomPairs(X, seed = 1)))
summary(mod3)
table(Class, mod3$classification)
adjustedRandIndex(Class, mod3$classification)

## ----------------------------------------------------------------------------
set.seed(20190603)
BIC <- NULL
for (i in 1:50)
{
  # get BIC table from initial random start
  BIC0 <- mclustBIC(X, verbose = FALSE,
                    initialization = list(hcPairs = hcRandomPairs(X)))
  # update BIC table by merging best BIC values for each
  # G and modelNames
  BIC  <- mclustBICupdate(BIC, BIC0)
}
summary(BIC, k = 5)

## ----------------------------------------------------------------------------
mod4 <- Mclust(X, x = BIC)
summary(mod4)
table(Class, mod4$classification)
adjustedRandIndex(Class, mod4$classification)

## ----------------------------------------------------------------------------
mclust.options("hcUse" = "SVD")  # restore the default
mod5 <- Mclust(X) # X is the unscaled flea data
# or specify only for this model fit
# mod5 <- Mclust(X, initialization = list(hcPairs = hc(X, use = "SVD")))
summary(mod5)
table(Class, mod5$classification)
adjustedRandIndex(Class, mod5$classification)  

## ----------------------------------------------------------------------------
data("iris", package = "datasets")
str(iris)
ms <- mstep(iris[, 1:4], modelName = "VVV",
            z = unmap(iris$Species))
str(ms, 1)
es <- estep(iris[, 1:4], modelName = "VVV", 
            parameters = ms$parameters)
str(es, 1)


# Chapter 4 ###################################################################

## ----------------------------------------------------------------------------
data("wdbc", package = "mclust")
X <- wdbc[, c("Texture_mean", "Area_extreme", "Smoothness_extreme")]
Class <- wdbc[, "Diagnosis"]

## ----------------------------------------------------------------------------
set.seed(123)
train <- sample(1:nrow(X), size = round(nrow(X)*2/3), replace = FALSE)
X_train <- X[train, ]
Class_train <- Class[train]
tab <- table(Class_train)
cbind(Counts = tab, "%" = prop.table(tab)*100)
X_test <- X[-train, ]
Class_test <- Class[-train]
tab <- table(Class_test)
cbind(Counts = tab, "%" = prop.table(tab)*100)

## ----------------------------------------------------------------------------
clp <- clPairs(X_train, Class_train, lower.panel = NULL)
clPairsLegend(0.1, 0.3, col = clp$col, pch = clp$pch, 
              class = ifelse(clp$class == "B", "Benign", "Malign"),
              title = "Breast cancer diagnosis:")

## ----------------------------------------------------------------------------
mod1 <- MclustDA(X_train, Class_train, modelType = "EDDA")
summary(mod1)

## ----------------------------------------------------------------------------
summary(mod1, parameters = TRUE)

## ----------------------------------------------------------------------------
summary(mod1, newdata = X_test, newclass = Class_test)

## ----------------------------------------------------------------------------
plot(mod1, what = "scatterplot")

## ----------------------------------------------------------------------------
plot(mod1, what = "error")

## ----------------------------------------------------------------------------
mod2 <- MclustDA(X_train, Class_train)
summary(mod2, newdata = X_test, newclass = Class_test)

## ----------------------------------------------------------------------------
plot(mod2, what = "scatterplot")

## ----------------------------------------------------------------------------
plot(mod2, what = "scatterplot", dimens = c(1, 2))

## ----------------------------------------------------------------------------
mod3 <- MclustDA(X_train, Class_train, G = 2, modelNames = "EEE")
summary(mod3, newdata = X_test, newclass = Class_test)

## ----------------------------------------------------------------------------
set.seed(20190520)
cv1 <- cvMclustDA(mod1)
str(cv1)
cv2 <- cvMclustDA(mod2)
str(cv2)

## ----------------------------------------------------------------------------
unlist(cv1[c("ce", "se.ce", "brier", "se.brier")])
unlist(cv2[c("ce", "se.ce", "brier", "se.brier")])

## ----------------------------------------------------------------------------
set.seed(2)
models <- mclust.options("emModelNames")
tab_CE <- tab_Brier <- 
  matrix(as.double(NA), nrow = length(models)+1, ncol = 5)
rownames(tab_CE) <- rownames(tab_Brier) <- 
  c(paste0("EDDA[", models, "]"), "MCLUSTDA")
colnames(tab_CE) <- colnames(tab_Brier) <- 
  c("Train", "10-fold CV", "se(CV)", "lower", "upper")
for (i in seq(models))
{
  mod <- MclustDA(X, Class, modelType = "EDDA", 
                  modelNames = models[i], verbose = FALSE)
  pred <- predict(mod, X)
  cv <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
  #
  tab_CE[i, 1] <- classError(pred$classification, Class)$errorRate
  tab_CE[i, 2] <- cv$ce
  tab_CE[i, 3] <- cv$se.ce
  tab_CE[i, 4] <- cv$ce - cv$se.ce
  tab_CE[i, 5] <- cv$ce + cv$se.ce
  #
  tab_Brier[i, 1] <- BrierScore(pred$z, Class)
  tab_Brier[i, 2] <- cv$brier
  tab_Brier[i, 3] <- cv$se.brier
  tab_Brier[i, 4] <- cv$brier - cv$se.brier
  tab_Brier[i, 5] <- cv$brier + cv$se.brier
}
i <- length(models)+1
mod <- MclustDA(X, Class, modelType = "MclustDA", verbose = FALSE)
pred <- predict(mod, X)
cv <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
#
tab_CE[i, 1] <- classError(pred$classification, Class)$errorRate
tab_CE[i, 2] <- cv$ce
tab_CE[i, 3] <- cv$se.ce
tab_CE[i, 4] <- cv$ce - cv$se.ce
tab_CE[i, 5] <- cv$ce + cv$se.ce
#
tab_Brier[i, 1] <- BrierScore(pred$z, Class)
tab_Brier[i, 2] <- cv$brier
tab_Brier[i, 3] <- cv$se.brier
tab_Brier[i, 4] <- cv$brier - cv$se.brier
tab_Brier[i, 5] <- cv$brier + cv$se.brier

## ----------------------------------------------------------------------------
tab_CE

## ----------------------------------------------------------------------------
library("ggplot2")
df <- data.frame(rownames(tab_CE), tab_CE)
colnames(df) <- c("model", "train", "cv", "se", "lower", "upper")
df$model <- factor(df$model, levels = rev(df$model))
ggplot(df, aes(x = model, y = cv, ymin = lower, ymax = upper)) +
  geom_point(aes(shape = "s1", color = "c1")) + 
  geom_errorbar(width = 0.5, col = "dodgerblue3") + 
  geom_point(aes(y = train, shape = "s2", color = "c2")) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.01), lim = c(0, NA)) +
  scale_color_manual(name = "", 
                     breaks = c("c1", "c2"),
                     values = c("dodgerblue3", "black"),
                     labels = c("CV", "Train")) +
  scale_shape_manual(name = "", 
                     breaks = c("s1", "s2"),
                     values = c(19, 0),
                     labels = c("CV", "Train")) +
  ylab("Classification error") + xlab("") + coord_flip() +
  theme(legend.position = "top")

## ----------------------------------------------------------------------------
tab_Brier

## ----------------------------------------------------------------------------
df <- data.frame(rownames(tab_Brier), tab_Brier)
colnames(df) <- c("model", "train", "cv", "se", "lower", "upper")
df$model <- factor(df$model, levels = rev(df$model))
ggplot(df, aes(x = model, y = cv, ymin = lower, ymax = upper)) +
  geom_point(aes(shape = "s1", color = "c1")) + 
  geom_errorbar(width = 0.5, col = "dodgerblue3") + 
  geom_point(aes(y = train, shape = "s2", color = "c2")) +
  scale_y_continuous(breaks = seq(0, 0.2, by = 0.01), lim = c(0, NA)) +
  scale_color_manual(name = "", 
                     breaks = c("c1", "c2"),
                     values = c("dodgerblue3", "black"),
                     labels = c("CV", "Train")) +
  scale_shape_manual(name = "", 
                     breaks = c("s1", "s2"),
                     values = c(19, 0),
                     labels = c("CV", "Train")) +
  ylab("Brier score") + xlab("") + coord_flip() +
  theme(legend.position = "top")

## ----------------------------------------------------------------------------
# confusion matrix
(tab <- table(Predict = cv1$classification, Class = Class_train))
tab[2, 2]/sum(tab[, 2])  # sensitivity
tab[1, 1]/sum(tab[, 1])  # specificity

## ----------------------------------------------------------------------------
threshold <- seq(0, 1, by = 0.01)
sensitivity <- specificity <- rep(NA, length(threshold))
for(i in 1:length(threshold))
{
  pred <- factor(ifelse(cv1$z[, "M"] > threshold[i], "M", "B"),
                 levels = c("B", "M"))
  tab <- table(pred, Class_train)
  sensitivity[i] <- tab[2, 2]/sum(tab[, 2])
  specificity[i] <- tab[1, 1]/sum(tab[, 1])
}

## ----------------------------------------------------------------------------
auc_approx <- function(tpr, fpr)
{
  x <- 1 - fpr
  y <- tpr
  dx <- c(diff(x), 0)
  dy <- c(diff(y), 0)
  sum(y * dx) + sum(dy * dx)/2
}
auc_approx(tpr = sensitivity, fpr = 1 - specificity)

## ----------------------------------------------------------------------------
threshold <- seq(0, 1, by = 0.01)
sensitivity2 <- specificity2 <- rep(NA, length(threshold))
for(i in 1:length(threshold))
{
  pred <- factor(ifelse(cv2$z[, "M"] > threshold[i], "M", "B"),
                 levels = c("B", "M"))
  tab <- table(pred, Class_train)
  sensitivity2[i] <- tab[2, 2]/sum(tab[, 2])
  specificity2[i] <- tab[1, 1]/sum(tab[, 1])
}

## ----------------------------------------------------------------------------
plot(1-specificity, sensitivity, type = "l", lwd = 2)
abline(h = c(0, 1), v = c(0, 1), lty = 3)
abline(a = 0, b = 1, lty = 2)
plot(1-specificity2, sensitivity2, type = "l", lwd = 2,
     xlab = "1-specificity", ylab = "sensitivity") 
abline(h = c(0, 1), v = c(0, 1), lty = 3)
abline(a = 0, b = 1, lty = 2)

## ----------------------------------------------------------------------------
J <- sensitivity + specificity - 1 
threshold[which.max(J)]     # optimal threshold
sensitivity[which.max(J)]   # sensitivity at optimal threshold
specificity[which.max(J)]   # specificity at optimal threshold

## ----------------------------------------------------------------------------
data("bankruptcy", package = "MixGHD")
X <- bankruptcy[, -1]
Class <- factor(bankruptcy$Y, levels = c(1:0), 
                labels = c("solvent", "bankrupt"))
cl <- clPairs(X, Class)
legend("bottomright", legend = cl$class, 
       pch = cl$pch, col = cl $col, inset = 0.02)

## ----------------------------------------------------------------------------
mod <- MclustDA(X, Class, modelType = "EDDA")
summary(mod)

## ----------------------------------------------------------------------------
plot(mod, what = "scatterplot")
plot(mod, what = "error")

## ----------------------------------------------------------------------------
(C <- matrix(c(0, 1, 10, 0), nrow = 2, ncol = 2, byrow = TRUE))

## ----------------------------------------------------------------------------
rowSums(C)

## ----------------------------------------------------------------------------
pred <- predict(mod)
(tab <- table(Class, Predicted = pred$classification))
sum(tab * C)

## ----------------------------------------------------------------------------
pred <- predict(mod, prop = mod$prop*rowSums(C))
(tab <- table(Class, Predicted = pred$classification))
sum(tab * C)

## ----------------------------------------------------------------------------
set.seed(20190515)
# generate training data from a balanced case-control sample
n_train <- 1000
class_train <- factor(sample(0:1, size = n_train, prob = c(0.5, 0.5), 
                             replace = TRUE))
x_train <- ifelse(class_train == 1, rnorm(n_train, mean = 3, sd = 1), 
                                    rnorm(n_train, mean = 0, sd = 1))

hist(x_train[class_train == 0], breaks = 11, xlim = range(x_train), 
     main = "", xlab = "x", 
     col = adjustcolor("dodgerblue2", alpha.f = 0.5), border = "white")
hist(x_train[class_train == 1], breaks = 11, add = TRUE,
     col = adjustcolor("red3", alpha.f = 0.5), border = "white")
box()

# generate test data from mixture f(x) = 0.9 * N(0,1) + 0.1 * N(3,1)
n <- 10000
mixpro <- c(0.9, 0.1)
class_test <- factor(sample(0:1, size = n, prob = mixpro, 
                            replace = TRUE))
x_test <- ifelse(class_test == 1, rnorm(n, mean = 3, sd = 1), 
                                  rnorm(n, mean = 0, sd = 1))
hist(x_test[class_test == 0], breaks = 15, xlim = range(x_test), 
     main = "", xlab = "x", 
     col = adjustcolor("dodgerblue2", alpha.f = 0.5), border = "white")
hist(x_test[class_test == 1], breaks = 11, add = TRUE,
     col = adjustcolor("red3", alpha.f = 0.5), border = "white")
box()

## ----------------------------------------------------------------------------
mod <- MclustDA(x_train, class_train)
summary(mod, parameters = TRUE)

## ----------------------------------------------------------------------------
pred <- predict(mod, newdata = x_test)
classError(pred$classification, class_test)$error
BrierScore(pred$z, class_test)

## ----------------------------------------------------------------------------
priorProp <- seq(0.01, 0.99, by = 0.01)
CE <- BS <- rep(as.double(NA), length(priorProp))
for (i in seq(priorProp))
{
  pred <- predict(mod, newdata = x_test, 
                  prop = c(1-priorProp[i], priorProp[i]))
  CE[i] <- classError(pred$classification, class = class_test)$error
  BS[i] <- BrierScore(pred$z, class_test)
}

## ----------------------------------------------------------------------------
matplot(priorProp, cbind(CE, BS), type = "l", lty = 1, lwd = 2, xaxt = "n",
        xlab = "Class prior probability", ylab = "", ylim = c(0, max(CE, BS)), 
        col = c("red3", "dodgerblue3"),
        panel.first = 
          { abline(h = seq(0, 1, by = 0.05), col = "grey", lty = 3)
            abline(v = seq(0, 1, by = 0.05), col = "grey", lty = 3) 
          })
axis(side = 1, at = seq(0, 1, by = 0.1))
abline(v = mod$prop[2],             # training proportions
       lty = 2, lwd = 2)            
abline(v = mean(class_test == 1),   # test proportions (usually unknown)
       lty = 3, lwd = 2)   
legend("topleft", legend = c("ClassError", "BrierScore"),
       col = c("red3", "dodgerblue3"), lty = 1, lwd = 2, inset = 0.02)

## ----------------------------------------------------------------------------
(priorProbs <- classPriorProbs(mod, x_test))

## ----------------------------------------------------------------------------
pred <- predict(mod, newdata = x_test, prop = priorProbs)
classError(pred$classification, class = class_test)$error
BrierScore(pred$z, class_test)

## ----------------------------------------------------------------------------
(prior_test <- prop.table(table(class_test)))
pred <- predict(mod, newdata = x_test, prop = prior_test)
classError(pred$classification, class = class_test)$error
BrierScore(pred$z, class_test)

## ----------------------------------------------------------------------------
data("wdbc", package = "mclust")
x <- with(wdbc, 
    0.2322*Texture_mean + 0.01117*Area_extreme + 68.37*Smoothness_extreme)
Class <- wdbc[, "Diagnosis"]
mod <- MclustDA(x, Class, modelType = "MclustDA")
summary(mod)

## ----------------------------------------------------------------------------
(prop <- mod$prop)
col <- mclust.options("classPlotColors")
x0 <- seq(0, max(x)*1.1, length = 1000)
par1 <- mod$models[["B"]]$parameters
f1 <- dens(par1$variance$modelName, data = x0, parameters = par1)
par2 <- mod$models[["M"]]$parameters
f2 <- dens(par2$variance$modelName, data = x0, parameters = par2)
matplot(x0, cbind(prop[1]*f1, prop[2]*f2), type = "l", lty = 1, 
        col = col, ylab = "Class density", xlab = "x")
legend("topright", title = "Diagnosis:", legend = names(prop), 
       col = col, lty = 1, inset = 0.02)

## ----------------------------------------------------------------------------
set.seed(20190520)
cv <- cvMclustDA(mod)  # by default: prop = mod$prop
unlist(cv[c("ce", "se.ce")])

## ----------------------------------------------------------------------------
set.seed(20190520)
cv <- cvMclustDA(mod, prop = c(0.5, 0.5))
unlist(cv[c("ce", "se.ce")])

## ----------------------------------------------------------------------------
x0 <- seq(min(x), max(x), length.out = 1000)
pred <- predict(mod, newdata = x0) 
(threshold1 <- approx(pred$z[, 2], x0, xout = 0.5)$y)
pred <- predict(mod, newdata = x0, prop = c(0.5, 0.5))
(threshold2 <- approx(pred$z[, 2], x0, xout = 0.5)$y)

## ----------------------------------------------------------------------------
par(mfrow=c(1, 2))
plot(mod, what = "scatterplot", main = TRUE)
abline(v = threshold1, lty = 2)
plot(mod, what = "classification", main = TRUE)
abline(v = threshold1, lty = 2)

## ----------------------------------------------------------------------------
set.seed(1)
threshold <- seq(0.1, 0.9, by = 0.05)
ngrid <- length(threshold)
cv <- data.frame(threshold, error = numeric(ngrid))
cverr <- cvMclustDA(mod, verbose = FALSE)
for (i in seq(threshold))
{
  cv$error[i] <- classError(ifelse(cverr$z[, 2] > threshold[i], "M", "B"),
                            Class)$errorRate
}  
min(cv$error)
threshold[which.min(cv$error)]

ggplot(cv, aes(x = threshold, y = error)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  ylab("CV misclassification error") + 
  xlab("Probability threshold of malignant (M) tumor class")

## ----------------------------------------------------------------------------
set.seed(1)
priorProb <- seq(0.1, 0.9, by = 0.05)
ngrid <- length(priorProb)
cv_error2 <- data.frame(priorProb, 
                        cv = numeric(ngrid), 
                        lower = numeric(ngrid), 
                        upper = numeric(ngrid))
for (i in seq(priorProb))
{
  cv <- cvMclustDA(mod, prop = c(1-priorProb[i], priorProb[i]), 
                   verbose = FALSE)
  cv_error2$cv[i]    <- cv$ce
  cv_error2$lower[i] <- cv$ce - cv$se.ce
  cv_error2$upper[i] <- cv$ce + cv$se.ce
}  
min(cv_error2$cv)
priorProb[which.min(cv_error2$cv)]

ggplot(cv_error2, aes(x = priorProb, y = cv)) +
  geom_point() + 
  geom_linerange(aes(ymin = lower, ymax = upper)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  ylab("CV misclassification error") + 
  xlab("Malignant (M) tumor class prior probability")

## ----------------------------------------------------------------------------
n <- 200
pars <- list(pro = c(0.5, 0.5),
             mean = matrix(c(-1, 1), nrow = 2, ncol = 2, byrow = TRUE),
             variance = mclustVariance("EII", d = 2, G = 2))
pars$variance$sigmasq <- 1
data <- sim("EII", parameters = pars, n = n, seed = 12)
class <- data[, 1]
X <- data[, -1]
clPairs(X, class, symbols = c(1, 2))

# Randomly remove labels
cl <- class; cl[sample(1:n, size = 195)] <- NA
# table(cl, useNA = "ifany")
clPairs(X, ifelse(is.na(cl), 0, class),
        symbols = c(0, 16, 17), colors = c("grey", 4, 2))

## ----------------------------------------------------------------------------
mod      <- MclustDA(X, class, modelType = "EDDA") 
mod_EDDA <- MclustDA(X[!is.na(cl), ], cl[!is.na(cl)], modelType = "EDDA", modelNames = "EII")
mod_SSC  <- MclustSSC(X, cl)

ngrid = 100
xgrid = seq(-4, 4, length.out = ngrid)
ygrid = seq(-5, 5, length.out = ngrid)
xygrid = expand.grid(xgrid, ygrid)

pred_mod  <- predict(mod, newdata = xygrid)
pred_EDDA <- predict(mod_EDDA, newdata = xygrid)
pred_SSC  <- predict(mod_SSC, newdata = xygrid)

col = mclust.options("classPlotColors")[class]
pch = class
pch[!is.na(cl)] = ifelse(cl[!is.na(cl)] == 1, 19, 17)

plot(X, pch = pch, col = col)
contour(xgrid, ygrid, matrix(pred_mod$z[, 1], ngrid, ngrid), 
        add = TRUE, levels = 0.5, drawlabels = FALSE, lwd = 2, col = "cyan2")
contour(xgrid, ygrid, matrix(pred_EDDA$z[, 1], ngrid, ngrid), 
        add = TRUE, levels = 0.5, drawlabels = FALSE, lty = 2, lwd = 2)
contour(xgrid, ygrid, matrix(pred_SSC$z[, 1], ngrid, ngrid), 
        add = TRUE, levels = 0.5, drawlabels = FALSE, lty = 3, lwd = 2)

## ----------------------------------------------------------------------------
data("olive", package = "pgmm")
X <- olive[, 3:10]
class <- factor(olive$Region, levels = 1:3, 
                labels = c("South", "Sardinia", "North"))
table(class)

## ----------------------------------------------------------------------------
mod_EDDA_full <- MclustDA(X, class, modelType = "EDDA")
pred_EDDA_full <- predict(mod_EDDA_full, newdata = X)
classError(pred_EDDA_full$classification, class)$errorRate
BrierScore(pred_EDDA_full$z, class)

## ----------------------------------------------------------------------------
set.seed(20200807)
pct_labeled_data <- 10
n <- nrow(X)
cl <- class
is.na(cl) <- sample(1:n, round(n*(1-pct_labeled_data/100)))
table(cl, useNA = "ifany")

## ----------------------------------------------------------------------------
mod_SSC <- MclustSSC(X, cl)
plot(mod_SSC, what = "BIC")
mod_SSC$BIC
pickBIC(mod_SSC$BIC, 5) - max(mod_SSC$BIC)  # BIC diff for the top-5 models

## ----------------------------------------------------------------------------
summary(mod_SSC)

## ----------------------------------------------------------------------------
pred_SSC <- predict(mod_SSC, newdata = X[is.na(cl), ])
table(Predicted = pred_SSC$classification, Class = class[is.na(cl)])
classError(pred_SSC$classification, class[is.na(cl)])$errorRate
BrierScore(pred_SSC$z, class[is.na(cl)])

## ----------------------------------------------------------------------------
set.seed(20201013)
pct_labeled_data <- c(5, seq(10, 90, by = 10), 95)
BS <- matrix(as.double(NA), nrow = length(pct_labeled_data), ncol = 2,
             dimnames = list(pct_labeled_data, c("EDDA", "SSC")))
for (i in seq(pct_labeled_data))
{
  cl <- class
  labeled <- sample(1:n, round(n*pct_labeled_data[i]/100))
  cl[-labeled] <- NA
  # Classification on labeled data
  mod_EDDA  <- MclustDA(X[labeled, ], cl[labeled], 
                        modelType = "EDDA")
  # prediction for the unlabeled data
  pred_EDDA <- predict(mod_EDDA, newdata = X[-labeled, ])
  BS[i, 1]  <- BrierScore(pred_EDDA$z, class[-labeled])
  # Semi-supervised classification
  mod_SSC  <- MclustSSC(X, cl)
  # prediction for the unlabeled data
  pred_SSC <- predict(mod_SSC, newdata = X[-labeled, ])
  BS[i, 2] <- BrierScore(pred_SSC$z, class[-labeled])
}
BS

matplot(pct_labeled_data, BS, type = "b", 
        lty = 1, pch = c(19, 15), col = c(2, 4), xaxt = "n",
        xlab = "Percentage of labeled data", ylab = "Brier score")
axis(side = 1, at = pct_labeled_data)
abline(h = BrierScore(pred_EDDA_full$z, class), lty = 2)
legend("topright", pch = c(19, 15), col = c(2, 4), lty = 1, 
       legend = c("EDDA", "SSC"), inset = 0.02)


# Chapter 5 ###################################################################

## ----------------------------------------------------------------------------
data("stamps1", package = "multimode")
str(stamps1)
Thickness <- stamps1$thickness

## ----------------------------------------------------------------------------
dens <- densityMclust(Thickness)

## ----------------------------------------------------------------------------
summary(dens, parameters = TRUE)

## ----------------------------------------------------------------------------
br <- seq(min(Thickness), max(Thickness), length.out = 21)
plot(dens, what = "density", data = Thickness, breaks = br)

## ----------------------------------------------------------------------------
with(dens$parameters, 
     data.frame(mean = mean,
                sd = sqrt(variance$sigmasq),
                CoefVar = sqrt(variance$sigmasq)/mean*100))

## ----------------------------------------------------------------------------
x <- c(0.07, 0.08, 0.1, 0.12)
predict(dens, newdata = x, what = "dens")
predict(dens, newdata = x, logarithm = TRUE)
predict(dens, newdata = x, what = "cdens")

## ----------------------------------------------------------------------------
predict(dens, newdata = x, what = "z")

## ----------------------------------------------------------------------------
Year <- stamps1$year
table(Year)
h1 <- hist(Thickness[Year == "1872"], breaks = br, plot = FALSE)
h1$density <- h1$density*prop.table(table(Year))[1]
h2 <- hist(Thickness[Year == "1873-1874"], breaks = br, plot = FALSE)
h2$density <- h2$density*prop.table(table(Year))[2]
x <- seq(min(Thickness)-diff(range(Thickness))/10, 
         max(Thickness)+diff(range(Thickness))/10, length = 200)
cdens <- predict(dens, x, what = "cdens")
cdens <- sweep(cdens, 2, dens$parameters$pro, "*")
col <- adjustcolor(mclust.options("classPlotColors")[1:2], alpha = 0.3)
ylim <- range(h1$density, h2$density, cdens)
plot(h1, xlab = "Thickness", freq = FALSE, main = "", border = "white", 
col = col[1], xlim = range(x), ylim = ylim)
plot(h2, add = TRUE, freq = FALSE, border = "white", col = col[2])
matplot(x, cdens, type = "l", lwd = 1, lty = 1, col = 1, add = TRUE)
box()
legend("topright", legend = levels(Year), col = col, pch = 15, inset = 0.02,
       title = "Overprinted years:", title.adj = 0.2)

## ----------------------------------------------------------------------------
data("acidity", package = "BNPdensity")

## ----------------------------------------------------------------------------
summary(mclustBIC(acidity), k = 5)

## ----------------------------------------------------------------------------
dens_E2 <- densityMclust(acidity, G = 2, modelNames = "E")
summary(dens_E2, parameters = TRUE)

## ----------------------------------------------------------------------------
mclustBootstrapLRT(acidity, modelName = "V")

## ----------------------------------------------------------------------------
dens_V3 <- densityMclust(acidity, G = 3, modelNames = "V")
summary(dens_V3, parameters = TRUE)

## ----------------------------------------------------------------------------
plot(dens_E2, what = "density",
     ylim = c(0, max(dens_E2$density, dens_V3$density)))
rug(acidity)
plot(dens_V3, what = "density", 
     ylim = c(0, max(dens_E2$density, dens_V3$density)))
rug(acidity)

## ----------------------------------------------------------------------------
plot(dens_E2, what = "diagnostic", type = "cdf")
plot(dens_E2, what = "diagnostic", type = "qq")

## ----------------------------------------------------------------------------
plot(dens_V3, what = "diagnostic", type = "cdf")
plot(dens_V3, what = "diagnostic", type = "qq")

## ----------------------------------------------------------------------------
data("faithful", package = "datasets")
plot(faithful)
dens <- densityMclust(faithful)
summary(dens, parameters = TRUE)

## ----------------------------------------------------------------------------
plot(faithful)

## ----------------------------------------------------------------------------
plot(dens, what = "dens")

## ----------------------------------------------------------------------------
plot(dens, what = "density", type = "image")

## ----------------------------------------------------------------------------
plot(dens, what = "density", type = "persp")

## ----------------------------------------------------------------------------
plot(dens, what = "density")

## ----------------------------------------------------------------------------
plot(dens, what = "density", type = "image")
plot(dens, what = "density", type = "persp")

## ----------------------------------------------------------------------------
plot(dens, what = "density", type = "hdr")

## ----------------------------------------------------------------------------
data("aircraft", package = "sm")
X <- log(subset(aircraft, subset = (Period == 3), select = 3:8))
PCA <- prcomp(X, scale = TRUE)
summary(PCA)
PCA$rotation        # loadings
Z <- PCA$x[, 1:3]   # PCA projection

## ----------------------------------------------------------------------------
BIC <- mclustBIC(Z)
summary(BIC, k = 5)

## ----------------------------------------------------------------------------
densAircraft <- densityMclust(Z, G = 3, modelNames = "VVE", plot = FALSE)
plot(densAircraft, what = "density", type = "hdr", 
     data = Z, points.cex = 0.5)

## ----------------------------------------------------------------------------
library("mclustAddons")

## ----------------------------------------------------------------------------
data("suicide", package = "mclustAddons")
dens <- densityMclust(suicide)
rug(suicide)            # add data points at the bottom of the graph
abline(v = 0, lty = 3)  # draw a vertical line at the natural boundary


## ----------------------------------------------------------------------------
bdens <- densityMclustBounded(suicide, lbound = 0)
summary(bdens, parameters = TRUE)
plot(bdens, what = "density")
rug(suicide)            # add data points at the bottom of the graph
abline(v = 0, lty = 3)  # draw a vertical line at the natural boundary

## ----------------------------------------------------------------------------
plot(dens, what = "density")
rug(suicide)
abline(v = 0, lty = 3)
#
plot(bdens, what = "density")
rug(suicide)
abline(v = 0, lty = 3)

## ----------------------------------------------------------------------------
data("racial", package = "mclustAddons")
bdens <- densityMclustBounded(racial$PropWhite, lbound = 0, ubound = 1)
plot(bdens, what = "density", 
     lwd = 2, col = "dodgerblue2",
     data = racial$PropWhite, breaks = 15,
     xlab = "Proportion of white students enrolled in schools")
rug(racial$PropWhite)        # add data points at the bottom of the graph
abline(v = c(0, 1), lty = 3)  # draw a vertical line at the natural boundary

## ----------------------------------------------------------------------------
f <- function(x) 
  0.7*dnorm(x, mean = 0, sd = 1) + 0.3*dnorm(x, mean = 4, sd = 1)
curve(f, from = -4, to = 8)

## ----------------------------------------------------------------------------
set.seed(20190525)
par <- list(pro = c(0.7, 0.3), mean = c(0, 4), 
            variance = mclustVariance("E", G = 2))
par$variance$sigmasq <- c(1, 1)
x <- sim(modelName = "E", parameters = par, n = 1e4)[, -1]

## ----------------------------------------------------------------------------
par(mfrow = c(2, 2))
prob <- c(0.25, 0.5, 0.75, 0.95)
(hdr <- hdrlevels(f(x), prob))
for (j in seq(prob))
{
  curve(f, from = -4, to = 8)
  mtext(side = 3, paste0(prob[j]*100, "% HDR"), adj = 0)
  abline(h = hdr[j], lty = 2)
  rug(x, col = "lightgrey")
  rug(x[f(x) >= hdr[j]])
}

## ----------------------------------------------------------------------------
dens <- densityMclust(x, plot = FALSE)
hdrlevels(predict(dens, x), prob)


# Chapter 6 ###################################################################

## ----------------------------------------------------------------------------
data("Snapper", package = "FSAdata")
x <- Snapper[,1]
mod <- Mclust(x, G = 4, modelNames = "V")
summary(mod, parameters = TRUE)

## ----------------------------------------------------------------------------
mclust1Dplot(x, what = "classification", 
             parameters = mod$parameters, z = mod$z,
             xlab = "Fish length")
mclust1Dplot(x, what = "uncertainty",
             parameters = mod$parameters, z = mod$z,
             xlab = "Fish length")

## ----------------------------------------------------------------------------
mod <- Mclust(faithful)
mclust2Dplot(data = faithful, what = "classification", 
             parameters = mod$parameters, z = mod$z)
mclust2Dplot(data = faithful, what = "uncertainty",
             parameters = mod$parameters, z = mod$z)

## ----------------------------------------------------------------------------
surfacePlot(data = faithful, parameters = mod$parameters,
            what = "density", type = "contour",
            transformation = "log")

## ----------------------------------------------------------------------------
surfacePlot(data = faithful, parameters = mod$parameters,
            what = "density", type = "image")

## ----------------------------------------------------------------------------
surfacePlot(data = faithful, parameters = mod$parameters,
            what = "density", type = "persp")

## ----------------------------------------------------------------------------
surfacePlot(data = faithful, parameters = mod$parameters,
            what = "uncertainty", type = "image",
            transformation = "sqrt")

## ----------------------------------------------------------------------------
data("iris", package = "datasets")
mod <- Mclust(iris[,1:4], G = 3)
summary(mod)

## ----------------------------------------------------------------------------
coordProj(data = iris[,1:4], dimens = c(2,4), what = "classification",
          parameters = mod$parameters, z = mod$z)
coordProj(data = iris[,1:4], dimens = c(2,4), what = "uncertainty",
          parameters = mod$parameters, z = mod$z)
coordProj(data = iris[,1:4], dimens = c(2,4), what = "error",
          parameters = mod$parameters, z = mod$z, truth = iris$Species)

## ----------------------------------------------------------------------------
coordProj(data = iris[,1:4], dimens = c(2,4), what = "classification",
          parameters = mod$parameters, z = mod$z)

## ----------------------------------------------------------------------------
coordProj(data = iris[,1:4], dimens = c(2,4), what = "uncertainty",
          parameters = mod$parameters, z = mod$z)

## ----------------------------------------------------------------------------
coordProj(data = iris[,1:4], dimens = c(2,4), what = "error",
          parameters = mod$parameters, z = mod$z, truth = iris$Species)

## ----------------------------------------------------------------------------
set.seed(1)
nrow(iris)
Q <- randomOrthogonalMatrix(ncol(iris[,1:4]), 2)
dim(Q)
QTQ <- crossprod(Q)           # equivalently t(Q) %*% Q
zapsmall(QTQ)                 # 2 x 2 identity matrix
# projection of iris data onto coordinates Q
irisProj <- as.matrix(iris[,1:4]) %*% Q  

## ----------------------------------------------------------------------------
randProj(data = iris[,1:4], seeds = c(1,13,79,201), 
         what = "classification",
         parameters = mod$parameters, z = mod$z)

## ----------------------------------------------------------------------------
randProj(data = iris[,1:4], seed = 1, what = "classification",
         parameters = mod$parameters, z = mod$z)

## ----------------------------------------------------------------------------
randProj(data = iris[,1:4], seed = 13, what = "classification",
         parameters = mod$parameters, z = mod$z)

## ----------------------------------------------------------------------------
randProj(data = iris[,1:4], seed = 79, what = "classification",
         parameters = mod$parameters, z = mod$z)

## ----------------------------------------------------------------------------
randProj(data = iris[,1:4], seed = 201, what = "classification",
         parameters = mod$parameters, z = mod$z)

## ----------------------------------------------------------------------------
plot(crimcoords(iris[,1:4], mod$classification))

## ----------------------------------------------------------------------------
data("thyroid", package = "mclust")
CRIMCOORDS <- crimcoords(thyroid[,-1], thyroid$Diagnosis)
summary(CRIMCOORDS)

## ----------------------------------------------------------------------------
plot(CRIMCOORDS)
points(CRIMCOORDS$means %*% CRIMCOORDS$basis, pch = 3, cex = 1.5, lwd = 2)
legend("topright", legend = levels(thyroid$Diagnosis), inset = 0.02,
       col = mclust.options("classPlotColors")[1:3],
       pch = mclust.options("classPlotSymbols")[1:3])

## ----------------------------------------------------------------------------
data("wine", package = "gclus")
Class <- factor(wine$Class, levels = 1:3,
                labels = c("Barolo", "Grignolino", "Barbera"))
X <- data.matrix(wine[,-1])
mod <- Mclust(X, G = 3, modelNames = "VVE")
table(Class, Cluster = mod$classification)

## ----------------------------------------------------------------------------
drmod <- MclustDR(mod, lambda = 1)
summary(drmod)

## ----------------------------------------------------------------------------
plot(drmod, what = "contour")

## ----------------------------------------------------------------------------
plot(drmod, what = "boundaries")

## ----------------------------------------------------------------------------
mod <- Mclust(iris[, 1:4], G = 3)
drmod <- MclustDR(mod, lambda = .5)
summary(drmod)

## ----------------------------------------------------------------------------
plot(drmod, what = "evalues")

## ----------------------------------------------------------------------------
plot(drmod, what = "pairs")

## ----------------------------------------------------------------------------
sdrmod <- MclustDRsubsel(drmod, verbose = TRUE)
summary(sdrmod)

## ----------------------------------------------------------------------------
zapsmall(cor(drmod$dir, sdrmod$dir)^2)

## ----------------------------------------------------------------------------
plot(sdrmod, what = "contour", nlevels = 7)

## ----------------------------------------------------------------------------
plot(sdrmod, what = "classification")

## ----------------------------------------------------------------------------
plot(sdrmod, what = "boundaries")

## ----------------------------------------------------------------------------
plot(sdrmod, what = "density")

## ----------------------------------------------------------------------------
data("banknote", package = "mclust")
mod <- MclustDA(data = banknote[, -1], class = banknote$Status)
summary(mod)

## ----------------------------------------------------------------------------
drmod <- MclustDR(mod, lambda = .5)
summary(drmod)

## ----------------------------------------------------------------------------
plot(drmod, what = "evalues")

## ----------------------------------------------------------------------------
plot(drmod, what = "pairs", lower.panel = NULL)
clPairsLegend(0.1, 0.4, class = levels(drmod$classification), 
              col = mclust.options("classPlotColors")[1:2],
              pch = mclust.options("classPlotSymbols")[1:2],
              title = "Swiss banknote data")

## ----------------------------------------------------------------------------
sdrmod <- MclustDRsubsel(drmod, verbose = TRUE)
summary(sdrmod)

## ----------------------------------------------------------------------------
zapsmall(cor(drmod$dir, sdrmod$dir)^2)

## ----------------------------------------------------------------------------
plot(sdrmod, what = "contour", nlevels = 15)
plot(sdrmod, what = "classification")

## ----------------------------------------------------------------------------
mod <- Mclust(faithful)
DF <- data.frame(mod$data, cluster = factor(mod$classification))
library("ggplot2")
ggplot(DF, aes(x = eruptions, y = waiting, 
               colour = cluster, shape = cluster)) +
  geom_point()

## ----------------------------------------------------------------------------
library("tidyr")
DF <- data.frame(mod$BIC[], G = 1:nrow(mod$BIC))
DF <- pivot_longer(DF, cols = 1:14, names_to = "Model", values_to = "BIC")
DF$Model <- factor(DF$Model, levels = mclust.options("emModelNames"))
ggplot(DF, aes(x = G, y = BIC, colour = Model, shape = Model)) +
  geom_point() + 
  geom_line() +
  scale_shape_manual(values = mclust.options("bicPlotSymbols")) +
  scale_color_manual(values = mclust.options("bicPlotColors")) +
  scale_x_continuous(breaks = unique(DF$G)) +
  xlab("Number of mixture components") +
  guides(shape = guide_legend(ncol=2))

## ----------------------------------------------------------------------------
mod <- Mclust(iris[, 1:4], G = 3)

## ----------------------------------------------------------------------------
means <- data.frame(Profile = 1:mod$G, t(mod$parameters$mean))
means <- pivot_longer(means, cols = -1, 
                      names_to = "Variable",
                      values_to = "Mean")
means$Profile  <- factor(means$Profile)
means$Variable <- factor(means$Variable, 
                         levels = rownames(mod$parameters$mean))
means

ggplot(means, aes(Variable, Mean, group = Profile, 
                  shape = Profile, color = Profile)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = NULL, y = "Latent profiles means") +
  scale_color_manual(values = mclust.options("classPlotColors")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "bottom")

## ----------------------------------------------------------------------------
damod <- MclustDA(iris[, 1:4], iris$Species)
drmod <- MclustDR(damod)
DF1 <- data.frame(drmod$dir[, 1:2], class = damod$class)
DF2 <- do.call("rbind", by(DF1, DF1[, 3], 
                           function(x) x[chull(x), ]))
ggplot() + 
  geom_point(data = DF1, 
             aes(x = Dir1, y = Dir2, color = class, shape = class)) + 
  geom_polygon(data = DF2, 
               aes(x = Dir1, y = Dir2, fill = class), 
               alpha = 0.3) +
  scale_color_manual(values = mclust.options("classPlotColors")) +
  scale_fill_manual(values = mclust.options("classPlotColors")) +
  scale_shape_manual(values = mclust.options("classPlotSymbols"))

## ----------------------------------------------------------------------------
mod <- densityMclust(faithful$waiting, plot = FALSE)
x <- extendrange(faithful$waiting, f = 0.1)
x <- seq(x[1], x[2], length.out = 101)
pred <- data.frame(x, density = predict(mod, newdata = x))
ggplot(faithful, aes(waiting)) +
  geom_histogram(aes(y = stat(density)), bins = 15, 
                 fill = "slategray3", colour = "grey92") +
  geom_line(data = pred, aes(x, density))

## ----------------------------------------------------------------------------
data("hemophilia", package = "rrcov")
X <- hemophilia[, 1:2]
mod <- Mclust(X, G = 2, modelName = "VVV")
boot <- MclustBootstrap(mod, nboot = 999, type = "bs")
DF <- data.frame(mixcomp = rep(1:boot$G, each = boot$nboot), 
                 pro = as.vector(boot$pro))
ggplot(DF, aes(x = pro)) + 
  geom_histogram(aes(y = stat(density)), bins = 15, 
                 fill = "slategray3", colour = "grey92") +
  facet_grid(~ mixcomp) +
  xlab("Mixing proportions") +
  ylab("Density of bootstrap distribution") 

## ----------------------------------------------------------------------------
DF0 <- data.frame(rbind(cbind(mixcomp = 1, boot$mean[, , 1]), 
                        cbind(mixcomp = 2, boot$mean[, , 2])))
DF0 <- tidyr::pivot_longer(DF0, cols = 2:3, names_to = "variable", values_to = "mean")
DF <- rbind(
  data.frame("mixcomp"  = 1,
             "variable" = rep(colnames(boot$mean[, , 1]), 
                              each = dim(boot$mean)[1]),
             "mean"     = as.vector(boot$mean[, , 1])),
  data.frame("mixcomp"  = 2,
             "variable" = rep(colnames(boot$mean[, , 2]), 
                              each = dim(boot$mean)[1]),
             "mean"     = as.vector(boot$mean[, , 2])))
ggplot(DF, aes(x = mean)) +
   geom_histogram(aes(y = stat(density)), bins = 15,
                  fill = "slategray3", colour = "grey92") +
   facet_grid(mixcomp ~ variable, scales = "free_x") +
   xlab("Means of mixture") +
   ylab("Density of bootstrap distribution")

## ----------------------------------------------------------------------------
mclust.options("bicPlotColors")
mclust.options("classPlotColors")

## ----------------------------------------------------------------------------
palette.colors(palette = "Okabe-Ito")

## ----------------------------------------------------------------------------
# get and save default palettes
bicPlotColors <- mclust.options("bicPlotColors")
classPlotColors <- mclust.options("classPlotColors")
# set Okabe-Ito palette for use in mclust
bicPlotColors_Okabe_Ito <-
   palette.colors(palette = "Okabe-Ito")[c(9, 1, 2:8, 2:6, 9, 1)]
names(bicPlotColors_Okabe_Ito) <- names(bicPlotColors)
classPlotColorsWong <- palette.colors(palette = "Okabe-Ito")[-1]
mclust.options("bicPlotColors" = bicPlotColors_Okabe_Ito)
mclust.options("classPlotColors" = classPlotColorsWong)

## ----------------------------------------------------------------------------
mod <- Mclust(iris[,3:4])
plot(mod, what = "BIC")
plot(mod, what = "classification")

## ----------------------------------------------------------------------------
mclust.options("bicPlotColors" = bicPlotColors)
mclust.options("classPlotColors" = classPlotColors)


# Chapter 7 ###################################################################

## ----------------------------------------------------------------------------
data("chevron", package = "mclust")
summary(chevron)
noise <- with(chevron, class == "noise")
X <- chevron[,2:3]
plot(X, cex = 0.5)
plot(X, cex = 0.5, col = ifelse(noise, "grey", "black"), 
     pch = ifelse(noise, 3, 1))

## ----------------------------------------------------------------------------
library("prabclus")
nnc <- NNclean(X, k = 5)
table(nnc$z)
clPairs(X, nnc$z, colors = c("darkgrey", "black"), symbols = c(3, 1))

## ----------------------------------------------------------------------------
modNoise <- Mclust(X, initialization = list(noise = (nnc$z == 0)))
summary(modNoise$BIC)

## ----------------------------------------------------------------------------
summary(modNoise, parameters = TRUE)
addmargins(table(chevron$class, modNoise$classification), 2)
plot(modNoise, what = "classification")

## ----------------------------------------------------------------------------
clPairs(X, nnc$z, colors = c("darkgrey", "black"), symbols = c(3, 1))

## ----------------------------------------------------------------------------
plot(modNoise, what = "classification")

## ----------------------------------------------------------------------------
hypvol(X)

## ----------------------------------------------------------------------------
library("cluster")
ehull <- ellipsoidhull(as.matrix(X))
volume(ehull)
modNoise.ehull <- Mclust(X, Vinv = 1/volume(ehull),
                         initialization = list(noise = (nnc$z == 0)))
summary(modNoise.ehull)
tab <- table(chevron$class, modNoise.ehull$classification)
addmargins(tab, 2)

## ----------------------------------------------------------------------------
library("geometry")
chull <- convhulln(X, options = "FA")
chull$vol

## ----------------------------------------------------------------------------
modNoise.chull <- Mclust(X, Vinv = 1/chull$vol,
                         initialization = list(noise = (nnc$z == 0)))
summary(modNoise.chull)
tab <- table(chevron$class, modNoise.chull$classification)
addmargins(tab, 2)

## ----------------------------------------------------------------------------
library("covRobust")
nnve <- cov.nnve(X, k = 5)
table(nnve$classification)

## ----------------------------------------------------------------------------
modNoise.nnve <- Mclust(X, initialization = 
                        list(noise = (nnve$classification == 0)))
summary(modNoise.nnve$BIC)
summary(modNoise.nnve)
addmargins(table(chevron$class, modNoise.nnve$classification), 2)

## ----------------------------------------------------------------------------
set.seed(123)
(Sigma <- array(c(2, 0.5, 0.5, 0.5, 1, 0, 0, 0.1, 2, -0.5, -0.5, 0.5), 
                dim = c(2, 2, 3)))
var.decomp <- sigma2decomp(Sigma)
str(var.decomp)
par <- list(pro = c(1/3, 1/3, 1/3), 
            mean = cbind(c(0, 3), c(3, 0), c(-3, 0)),
            variance = var.decomp)
data <- sim(par$variance$modelName, parameters = par, n = 200)
noise <- matrix(runif(100, -10, 10), nrow = 50, ncol = 2)
X <- rbind(data[, 2:3], noise)
cluster <- c(data[, 1], rep(0, 50))
clPairs(X, ifelse(cluster == 0, 0, 1), 
        colors = "black", symbols = c(16, 1), cex = c(0.5, 1))

## ----------------------------------------------------------------------------
nnc <- NNclean(X, k = 5)
modNoise <- Mclust(X, initialization = list(noise = (nnc$z == 0)))
summary(modNoise$BIC)
summary(modNoise, parameters = TRUE)

## ----------------------------------------------------------------------------
clPairs(X, ifelse(cluster == 0, 0, 1), 
        colors = "black", symbols = c(16, 1), cex = c(0.5, 1))

## ----------------------------------------------------------------------------
plot(modNoise, what = "classification")

## ----------------------------------------------------------------------------
plot(modNoise, what = "classification")
table(cluster, Classification = modNoise$classification)

## ----------------------------------------------------------------------------
data("galaxies", package = "MASS")
# now fix a typographical error in the data
# see help("galaxies", package = "MASS")
galaxies[78] <- 26960 
galaxies <- galaxies / 1000

## ----------------------------------------------------------------------------
set.seed(20190410)
BIC <- NULL
for(i in 1:50)
{
  BIC0 <- mclustBIC(galaxies, verbose = FALSE,
                    initialization = list(hcPairs = hcRandomPairs(galaxies)))
  BIC  <- mclustBICupdate(BIC, BIC0)
}
summary(BIC, k = 5)
plot(BIC)

mod <- densityMclust(galaxies, x = BIC)
summary(mod, parameters = TRUE)

plot(mod, what = "density", data = galaxies, breaks = 11)
rug(galaxies)

## ----------------------------------------------------------------------------
plot(BIC)

## ----------------------------------------------------------------------------
plot(mod, what = "density", data = galaxies, breaks = 11)
rug(galaxies)

## ----------------------------------------------------------------------------
set.seed(20190411)
BICp <- NULL
for(i in 1:50)
{
  BIC0p <- mclustBIC(galaxies, verbose = FALSE, 
                     prior = priorControl(),
            initialization = list(hcPairs = hcRandomPairs(galaxies)))
  BICp  <- mclustBICupdate(BICp, BIC0p)
}
summary(BICp, k = 5)
plot(BICp)

modp <- densityMclust(galaxies, x = BICp)
summary(modp, parameters = TRUE)

plot(modp, what = "density", data = galaxies, breaks = 11)
rug(galaxies)

## ----------------------------------------------------------------------------
plot(BICp)

## ----------------------------------------------------------------------------
plot(modp, what = "density", data = galaxies, breaks = 11)
rug(galaxies)

## ----------------------------------------------------------------------------
defaultPrior(galaxies, G = 3, modelName = "V")

## ----------------------------------------------------------------------------
data("olive", package = "pgmm")
# recode of labels for Region and Area
Regions <- c("South", "Sardinia", "North")
Areas <- c("North Apulia", "Calabria", "South Apulia", "Sicily", 
           "Inland Sardinia", "Coastal Sardinia", "East Liguria", 
           "West Liguria", "Umbria")
olive$Region <- factor(olive$Region, levels = 1:3, labels = Regions)
olive$Area <- factor(olive$Area, levels = 1:9, labels = Areas)
with(olive, table(Area, Region))

## ----------------------------------------------------------------------------
X <- scale(olive[, 3:10])

BIC <- mclustBIC(X, G = 1:15)
summary(BIC)
plot(BIC, legendArgs = list(x = "bottomright", ncol = 5))

BICp <- mclustBIC(X, G = 1:15, prior = priorControl())
summary(BICp)
plot(BICp, legendArgs = list(x = "bottomright", ncol = 5))

## ----------------------------------------------------------------------------
plot(BIC, legendArgs = list(x = "bottomright", ncol = 5))

## ----------------------------------------------------------------------------
plot(BICp, legendArgs = list(x = "bottomright", ncol = 5))

## ----------------------------------------------------------------------------
mod <- Mclust(X, x = BIC)
summary(mod)

table(Region = olive$Region, cluster = mod$classification)
adjustedRandIndex(olive$Region, mod$classification)

table(Area = olive$Area, cluster = mod$classification)
adjustedRandIndex(olive$Area, mod$classification)

## ----------------------------------------------------------------------------
modp <- Mclust(X, x = BICp)
summary(modp)

table(Region = olive$Region, cluster = modp$classification)
adjustedRandIndex(olive$Region, modp$classification)

table(Area = olive$Area, cluster = modp$classification)
adjustedRandIndex(olive$Area, modp$classification)

## ----------------------------------------------------------------------------
data("Baudry_etal_2010_JCGS_examples", package = "mclust")
plot(ex4.1)

## ----------------------------------------------------------------------------
mod_ex4.1 <- Mclust(ex4.1)
summary(mod_ex4.1)

## ----------------------------------------------------------------------------
CLUSTCOMBI <- clustCombi(mod_ex4.1)
summary(CLUSTCOMBI)

## ----------------------------------------------------------------------------
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
plot(CLUSTCOMBI, ex4.1, what = "classification")

## ----------------------------------------------------------------------------
plot(CLUSTCOMBI, what = "tree")
plot(CLUSTCOMBI, what = "tree", type = "rectangle", yaxis = "step")

## ----------------------------------------------------------------------------
optimClust <- clustCombiOptim(CLUSTCOMBI, reg = 2, plot = TRUE)
str(optimClust)

## ----------------------------------------------------------------------------
optimClust <- clustCombiOptim(CLUSTCOMBI, reg = 2, plot = TRUE)

## ----------------------------------------------------------------------------
data("faithful", package = "datasets")
mod <- Mclust(faithful)
clPairs(faithful, mod$classification)
plot(as.densityMclust(mod), what = "density", add = TRUE)

## ----------------------------------------------------------------------------
GMMHD <- gmmhd(mod)
summary(GMMHD)

## ----------------------------------------------------------------------------
plot(GMMHD, what = "mode")
plot(GMMHD, what = "cores")
plot(GMMHD, what = "clusters")

## ----------------------------------------------------------------------------
clPairs(faithful, mod$classification)
plot(as.densityMclust(mod), what = "density", add = TRUE)

## ----------------------------------------------------------------------------
plot(GMMHD, what = "mode")

## ----------------------------------------------------------------------------
plot(GMMHD, what = "cores")

## ----------------------------------------------------------------------------
plot(GMMHD, what = "clusters")

## ----------------------------------------------------------------------------
data("yeast", package = "MixSAL")
dup <- duplicated(yeast)
sum(dup) # count replicated observations
mod <- Mclust(yeast[!dup, -1])
summary(mod)
adjustedRandIndex(yeast$Site[!dup], mod$classification)

## ----------------------------------------------------------------------------
GMMHD <- gmmhd(mod)
summary(GMMHD)

## ----------------------------------------------------------------------------
table(yeast$Site[!dup], cluster = GMMHD$cluster)
adjustedRandIndex(yeast$Site[!dup], GMMHD$cluster)

## ----------------------------------------------------------------------------
plot(GMMHD, what = "cores", 
     col = c("grey50", "red3", "dodgerblue2"), 
     pch = c(0, 17, 16))
surfacePlot(GMMHD$x, parameters = GMMHD$MclustDR$parameters,
            what = "density", type = "contour", nlevels = 15,
            col = "black", drawlabels = FALSE, add = TRUE)
plot(GMMHD, what = "clusters", 
     col = c("red3", "dodgerblue2"), pch = c(17, 16))

## ----------------------------------------------------------------------------
data("Baudry_etal_2010_JCGS_examples", package = "mclust")
mod <- Mclust(ex4.1)
GMMHD <- gmmhd(mod)
summary(GMMHD)
optimClust <- clustCombiOptim(clustCombi(mod), reg = 2)
table(GMMHD = GMMHD$cluster, CLUSTCOMBI = optimClust$cluster)

## ----------------------------------------------------------------------------
mod <- Mclust(faithful)
sim0 <- sim(modelName = mod$modelName, 
            parameters = mod$parameters,
            n = nrow(faithful), seed = 0)
sim1 <- sim(modelName = mod$modelName, 
            parameters = mod$parameters,
            n = nrow(faithful), seed = 1)

## ----------------------------------------------------------------------------
xlim <- range(c(faithful[, 1], sim0[, 2], sim1[, 2]))
ylim <- range(c(faithful[, 2], sim0[, 3], sim1[, 3]))
mclust2Dplot(data = sim0[, -1], parameters = mod$parameters,
             classification = sim0[, 1], xlim = xlim, ylim = ylim)
mclust2Dplot(data = sim1[, -1], parameters = mod$parameters,
             classification = sim1[, 1], xlim = xlim, ylim = ylim)

## ----------------------------------------------------------------------------
head(sim0)

## ----------------------------------------------------------------------------
par <- list(
  pro = c(0.5, 0.3, 0.2),
  mean = matrix(c(0, 0, 3, 3, -4, 1), nrow = 2, ncol = 3),
  variance = sigma2decomp(matrix(c(1, 0.6, 0.6, 1.5), nrow = 2, ncol = 2), 
                          G = 3))
str(par)

## ----------------------------------------------------------------------------
sim <- sim(modelName = "EEI", parameters = par, n = 10000, seed = 123)
cluster <- sim[, 1]
x <- sim[, 2:3]

## ----------------------------------------------------------------------------
data("stlouis", package = "mix")
x <- data.frame(stlouis[, -(1:3)], row.names = NULL)
table(complete.cases(x))
apply(x, 2, function(x) prop.table(table(complete.cases(x))))

## ----------------------------------------------------------------------------
library("ggplot2")
df <- data.frame(obs = rep(1:nrow(x), times = ncol(x)),
                 var = rep(colnames(x), each = nrow(x)),
                 missing = as.vector(is.na(x)))
ggplot(data = df, aes(x = var, y = obs)) +
  geom_tile(aes(fill = missing)) +
  scale_fill_manual(values = c("lightgrey", "black")) +
  labs(x = "Variables", y = "Observations") +
  theme_minimal() +
  theme(axis.text.x = element_text(margin = margin(b = 10))) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank())

## ----------------------------------------------------------------------------
ggplot(data = df, aes(x = var, y = obs)) +
  geom_tile(aes(fill = missing)) +
  scale_fill_manual(values = c("lightgrey", "black")) +
  labs(x = "Variables", y = "Observations") +
  theme_minimal() +
  theme(axis.text.x = element_text(margin = margin(b = 10))) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank())

## ----------------------------------------------------------------------------
ximp <- imputeData(x, seed = 123)

## ----------------------------------------------------------------------------
imputePairs(x, ximp)

## ----------------------------------------------------------------------------
imputePairs(x, ximp)

## ----------------------------------------------------------------------------
x <- as.matrix(iris[, 1:4])
mod <- Mclust(x, G = 3, modelNames = "EEE")
table(iris$Species, mod$classification)
adjustedRandIndex(iris$Species, mod$classification)
mod$parameters$mean  # component means

## ----------------------------------------------------------------------------
set.seed(20171211)
isMissing <- sample(c(TRUE, FALSE), size = prod(dim(x)), 
                    replace = TRUE, prob = c(0.1, 0.9))
x[isMissing] <- NA
table(cmpObs <- complete.cases(x))

## ----------------------------------------------------------------------------
set.seed(20171212)
nImp  <- 100
muImp <- array(NA, c(ncol(x), 3, nImp))
clImp <- array(NA, c(nrow(x), nImp))
for(i in 1:nImp)
{ 
  xImp <- imputeData(x, verbose = FALSE)
  modImp <- Mclust(xImp, G = 3, modelNames = "EEE", verbose = FALSE)
  if (i == 1) clImp[, i] <- modImp$classification
  mcl <- matchCluster(clImp[, 1], modImp$classification)
  clImp[, i]  <- mcl$cluster
  muImp[, , i] <- modImp$parameters$mean[, mcl$ord]
}

# majority rule
cl <- apply(clImp, 1, function(x) majorityVote(x)$majority)
table(iris$Species, cl)
adjustedRandIndex(iris$Species, cl)

# pooled estimate of cluster means
apply(muImp, 1:2, mean)

