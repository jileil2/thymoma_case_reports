##################### set working directory
setwd("/Users/jileilin/Desktop/Research/thymomal case reports/")

##################### load packages
library('grDevices')

#################### read data
data <- read.csv('data4.csv')[, -1]


#################### figure 4
# ordinary least squares
lmod <- lm(change ~ log(interval), data = data)

# prepare figures
subset <- with(data, which(is.na(change) == FALSE & is.na(interval) == FALSE))

# confidence shade
confshade <- function(x, ylo, yhi, col = 8) {
  #
  # Draw a [Shade]d [Conf]idence band.
  n <- length(x)
  for (i in 1:(n-1)) {
    polygon(c(x[i], x[i + 1], x[i + 1], x[i]),
            c(ylo[i], ylo[i + 1], yhi[i + 1], yhi[i]),
            col = col,
            border = F)
  }
}

# plot
plot(log(interval[subset]), fitted.values(lmod), type = 'l',
     main = 'Tumor Size against Interval',
     ylab = 'Reduction of tumor size (%)',
     xlab = 'log[Interval (months)]',
     ylim = c(-100, 0), lwd = 3,
     col = 'magenta')
box(lwd = 4)
points(y = change[subset], x = log(interval[subset]),
       col = 'magenta')
text(x = 2.8, y = -26, 'Adjusted R^2 = 0.12')
text(x = 2.49, y = -18, 'Slope = -6.54')
text(x = 2.50, y = -10, 'p-value = 0.01')
abline(v = -0.772, lty = 2)
abline(h = -50, lty = 2)

# prediction interval
pred_conf <- predict(lmod, 
                     interval = "confidence", 
                     level = 0.95)

# turn into data.frame
pred_conf_df <- as.data.frame(pred_conf)
confshade(log(interval[subset]),
          ylo = pred_conf_df$lwr,
          yhi = pred_conf_df$upr,
          col = adjustcolor("magenta", alpha.f = 0.3))
