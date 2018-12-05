install.packages("hexSticker", dependencies = TRUE)
# devtools::install_github("GuangchuangYu/hexSticker")

library(hexSticker)
library(mclust)

mod <- Mclust(faithful, G = 2)

sticker(
  expression(
  { 
    par(mar = c(0,0,0,0)+0.2, bg = "transparent", cex = 0.3)
    plot(mod, what = "classification", main = "", xlab = "", ylab = "", 
         axes = FALSE, lwd = 2)
    # box()
  }),
  s_x=1, s_y=.8, s_width=1, s_height=1,
  package="mclust", p_x=1, p_y=1.5, p_size=10,
  h_color="steelblue", h_fill="grey98", p_color="steelblue",
  filename="~/R/mclust/misc/mclust_sticker.png"
)


set.seed(20181125)
m <- "EVE"
G <- 3
var <- mclustVariance(m, d = 2, G = G)
var$scale <- 1
var$shape <- apply(matrix(c(3,1,1.1,1,1,2),2,3), 2, 
                   function(a) sigma2decomp(diag(a))$shape)
var$orientation <- apply(matrix(c(1,1,-1,1),2,2), 2, normalize)
par <- list(pro = rep(1/G, G),
            mean = matrix(c(-1,-1, -3,3, 3,3), 2, 3),
            variance = var)
x <- sim(m, parameters = par, n = 200)[,-1]
plot(x)
mod <- Mclust(x, modelNames = m, G = G)
summary(mod)
plot(mod, what = "classification", asp = 1)


sticker(
  expression(
  { 
    par(mar = c(0,0,0,0), bg = "transparent", cex = 0.3)
    plot(mod, what = "classification", main = "", xlab = "", ylab = "", 
         axes = FALSE, lwd = 1)
    # box()
  }),
  s_x=1, s_y=.7, s_width=1.2, s_height=1.2,
  package="mclust", p_x=1, p_y=1.5, p_size=9,
  h_color="steelblue", h_fill="grey98", p_color="steelblue",
  filename="~/R/mclust/misc/mclust_logo.png"
)
