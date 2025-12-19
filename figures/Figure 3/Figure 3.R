##################### set working directory
##################### change working directory to the places where the data live
setwd("/Users/jileilin/Desktop/Research/thymomal case reports/")

##################### load packages
library("strucchange")
library('nnet')
library('quantreg')
library('multcomp')
library('forestplot')
library('ggplot2')

#################### read data
data <- read.csv('data2.csv')[, -1]
n <- nrow(data)
death.ind <- ifelse(data$outcome == 'Death', 1, 0)
subset <- which(data$non.derm.IBM == '-')

#################### foresplot
# logistic regression
mod <- glm(death.ind ~ Sex + AChR + ct + Striated + rs + Age, 
           family = 'binomial', data = data, 
           subset = subset, control = glm.control(maxit = 50))
# odds ratio
odds <- coef(summary(mod))[2:6, 1]
odds[1] <- -odds[1]
ses <- coef(summary(mod))[2:6, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
# p-value
pvs <- coef(summary(mod))[2:6, 4]

# groupwise proportion
covariates <- data[, c('Sex', 'AChR', 'ct', 'Striated', 'rs')]
ppp <- c()
for (j in 1:5) {
  tab <- with(data, table(covariates[, j][subset], death.ind[subset]))
  tab <- cbind(tab, tab / rowSums(tab))
  if (j == 1) {
    ppp <- rbind(ppp, c(paste0(tab[1, 2], '/', sum(tab[1, 1:2]), ', ', round(tab[2, 4] * 100)),
                        paste0(tab[2, 2], '/', sum(tab[2, 1:2]), ', ', round(tab[1, 4] * 100))))
  } else {
    ppp <- rbind(ppp, c(paste0(tab[2, 2], '/', sum(tab[2, 1:2]), ', ', round(tab[2, 4] * 100)),
                        paste0(tab[1, 2], '/', sum(tab[1, 1:2]), ', ', round(tab[1, 4] * 100))))
  }
}

# prepare data
base_data <- tibble::tibble(mean  = odds[, 1],
                            lower = odds[, 2],
                            upper = odds[, 3],
                            OR = round(mean, 2),
                            pvs = round(pvs, 4),
                            propp = ppp[, 1],
                            propn = ppp[, 2],
                            variable = c('Female', 'AChR', 
                                         'Cardiovascular intervention',
                                         'StrAbs', 'Respiratory support'))
base_data |>
  forestplot(labeltext = c(variable, OR, pvs),
             align = c("l", "c", "c"),      
             xlim = c(0, 100),
             title = expression(bold("Entire Cohort")),   
             xlab = expression(bold("Ratio of death and hospice care rate")),
             xlog = TRUE,
             boxsize = 0.35,
             lwd.ci = 5,                       # Thicker CI lines
             xticks = c(0.1, 0.2, 0.5, 1, 5, 10, 25),
             zero = 1,       
             lty.zero = 2, 
             col = fpColors(zero = "black"),  # Location of vertical line
             lwd.zero = 4,   
             txt_gp = fpTxtGp(
               xlab = gpar(fontsize = 12),
               ticks = gpar(fontsize = 12),    # <--- Larger axis tick labels
               label = gpar(fontsize = 12)
             ) # Enlarge axis tick labels
  ) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_header(variable = c("", "Variable"),
                OR = c("", "Adj-OR"),
                pvs = c("", "p-value")) |>
  fp_set_zebra_style("#EFEFEF")

