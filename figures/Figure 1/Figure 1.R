##################### set working directory
##################### change working directory to the places where the data live
setwd("/Users/jileilin/Desktop/Research/thymomal case reports/")

##################### load packages
library("strucchange")
library('nnet')
library('quantreg')
library('multcomp')
library('forestplot')

#################### read data
data <- read.csv('data1.csv')[, -1]


#################### forest plot (Figure 3A)
thymoma.ind <- ifelse(data$thymoma == '+', 1, 0)
# covariates
death <- ifelse(data$outcome == 'Death', 1, 0)
covariates <- with(data, data.frame(Sex, EOMG, Oplus, 
                                    Striated, GG, Cardiac, death))
covariates$GG[which(covariates$GG == 'nobiopsy')] <- NA
# subset
subset <- which(data$year >= 1993)

# plot
mean <- c()
lower <- c()
upper <- c()
pvs <- c()
for (j in 1:7) {
  
  # response
  if (j == 1) {
    y.j <- ifelse(covariates[, j] == names(table(covariates[, j]))[1], 1, 0)
  } else {
    y.j <- ifelse(covariates[, j] == names(table(covariates[, j]))[2], 1, 0)
  } 
  
  # logistic regression
  if (j >= 7) {
    mod.j <- glm(y.j ~ thymoma.ind + Age.TP + Sex, family = 'binomial',
                 subset = subset, data = data, control = glm.control(maxit = 50))
  } else if (j == 1) {
    mod.j <- glm(y.j ~ thymoma.ind + Age.TP, family = 'binomial',
                 data = data, control = glm.control(maxit = 50))
  } else {
    mod.j <- glm(y.j ~ thymoma.ind + Age.TP + Sex, family = 'binomial',
                 data = data, control = glm.control(maxit = 50))
  }
  
  # confidence intervals
  cis.j <- apply(confint(mod.j), 2, exp)[2, ]
  mean[j] <- exp(coef(mod.j))[-1][1]
  lower[j] <- cis.j[1]
  upper[j] <- cis.j[2]
  
  # p-values
  pvs[j] <- coef(summary(mod.j))[2, 4]
}

# p-value
pvs <- round(pvs, 4)


# groupwise proportion
ppp <- c()
for (j in 1:7) {
  if (j == 7) {
    tab <- with(data, table(thymoma.ind[subset], covariates[, j][subset]))
  } else {
    tab <- with(data, table(thymoma.ind, covariates[, j]))
  }
  tab <- cbind(tab, tab / rowSums(tab))
  if (j == 1) {
    ppp <- rbind(ppp, c(paste0(tab[2, 1], '/', sum(tab[2, 1:2]), ', ', round(tab[2, 4] * 100)),
                        paste0(tab[1, 1], '/', sum(tab[1, 1:2]), ', ', round(tab[1, 4] * 100))))
  } else {
    ppp <- rbind(ppp, c(paste0(tab[2, 2], '/', sum(tab[2, 1:2]), ', ', round(tab[2, 4] * 100)),
                        paste0(tab[1, 2], '/', sum(tab[1, 1:2]), ', ', round(tab[1, 4] * 100))))
  }
}

# forest plot
base_data <- tibble::tibble(mean  = mean, lower = lower, upper = upper,
                            OR = round(mean, 2), pvs = pvs,
                            propp = ppp[, 1], propn = ppp[, 2],
                            variable = c('Female', 'LOMG', 
                                         'MG-like Ocular symptoms', 
                                         'StrAb', 'Giant cells/granulomas', 
                                         'Cardiac involvement', 'Death'))
base_data |>
  forestplot(labeltext = c(variable, propp, propn, OR, pvs),
             align = c("l", "c", "c", "c", "c"),     # <-- Align columns
             xlim = c(0, 25),
             title = expression(bold("\\large A. Thymoma")),
             xlab = expression(bold("\\Large Association with thymoma versus no thymoma")),
             xlog = TRUE,
             boxsize = 0.35,
             lwd.ci = 5,                       # Thicker CI lines
             xticks = c(0.2, 0.5, 1, 2, 10, 40),
             zero = 1,       
             lty.zero = 2, 
             col = fpColors(zero = "black"),  # Location of vertical line
             lwd.zero = 4,                     # <--- Thicker vertical line
             txt_gp = fpTxtGp(
               xlab = gpar(fontsize = 12),
               ticks = gpar(fontsize = 12),    # <--- Larger axis tick labels
               label = gpar(fontsize = 12)
             ) # Enlarge axis tick labels
  ) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_header(variable = c("", "Variable\n"),
                propp = c("", "Prevalence \n Thymoma+"),
                propn = c("", "\n Thymoma-"),
                OR = c("", "Adj-OR\n"),
                pvs = c("", "\\emph{p}-value\n")) |>
  fp_set_zebra_style("#EFEFEF")


#################### forest plot (Figure 3B)
# subset
subset <- which(data$myocarditis_subset == 1)
subset2 <- which(data$year >= 1993 & data$myocarditis_subset == 1)
myocarditis.ind <- ifelse(data$myocarditis == 'Myocarditis', 1, 0)

# covariates
death <- ifelse(data$outcome == 'Death', 1, 0)
covariates <- with(data, data.frame(scale(Age), Sex, EOMG, thymoma, GG, 
                                    cort, imm, rs,
                                    death))
covariates$GG[which(covariates$GG == 'nobiopsy')] <- NA

# logistic regression
mean <- c()
lower <- c()
upper <- c()
pvs <- c()
for (j in 1:9) {
  
  # response
  if (j == 2) {
    y.j <- ifelse(covariates[, j] == names(table(covariates[, j]))[1], 1, 0)
  } else if (j == 1) {
    y.j <- covariates[, j]
  } else {
    y.j <- ifelse(covariates[, j] == names(table(covariates[, j]))[2], 1, 0)
  } 
  
  
  # logistic regression
  if (j >= 6) {
    mod.j <- glm(y.j ~ myocarditis.ind + Age.TP + Sex, family = 'binomial',
                 subset = subset2, data = data, control = glm.control(maxit = 50))
  } else if (j == 1) {
    mod.j <- mod.j <- lm(y.j ~ myocarditis.ind + Sex, data = data,
                         subset = subset)
  } else if (j == 2) {
    mod.j <- glm(y.j ~ myocarditis.ind + Age.TP, family = 'binomial',
                 subset = subset,
                 data = data, control = glm.control(maxit = 50))
  } else {
    mod.j <- glm(y.j ~ myocarditis.ind + Age.TP + Sex, family = 'binomial',
                 subset = subset,
                 data = data, control = glm.control(maxit = 50))
  }
  
  cis.j <- apply(confint(mod.j), 2, exp)[2, ]
  
  # confidence intervals
  mean[j] <- exp(coef(mod.j))[-1][1]
  lower[j] <- cis.j[1]
  upper[j] <- cis.j[2]
  
  # p-value
  pvs[j] <- coef(summary(mod.j))[2, 4]
}
pvs <- round(pvs, 4)


# groupwise proportion
ppp <- c()
for (j in 1:9) {
  if (j >= 6) {
    tab <- with(data, table(myocarditis.ind[subset2], covariates[, j][subset2]))
    tab <- cbind(tab, tab / rowSums(tab))
  } else if (j == 1) {
    med <- round(with(data, tapply(Age[subset], myocarditis.ind[subset], mean))[c(2, 1)])
    iqrs <- with(data, tapply(Age[subset], myocarditis.ind[subset], function(x) paste0(round(quantile(x, .25)), 
                                                                                       ' to ', round(quantile(x, .75)))))
    tab <- paste0(med, ' (', iqrs, ')')
    
  } else {
    tab <- with(data, table(myocarditis.ind[subset], covariates[, j][subset]))
    tab <- cbind(tab, tab / rowSums(tab))
  }
  if (j == 1) {
    ppp <- rbind(ppp, tab)
  } else if (j == 2) {
    ppp <- rbind(ppp, c(paste0(tab[2, 1], '/', sum(tab[2, 1:2]), ', ', round(tab[2, 3] * 100)),
                        paste0(tab[1, 1], '/', sum(tab[1, 1:2]), ', ', round(tab[1, 3] * 100))))
  } else {
    ppp <- rbind(ppp, c(paste0(tab[2, 2], '/', sum(tab[2, 1:2]), ', ', round(tab[2, 4] * 100)),
                        paste0(tab[1, 2], '/', sum(tab[1, 1:2]), ', ', round(tab[1, 4] * 100))))
  }
}

# forest plot
base_data <- tibble::tibble(mean  = mean,
                            lower = lower,
                            upper = upper,
                            OR = round(mean, 2),
                            propp = ppp[, 1],
                            propn = ppp[, 2],
                            pvs = pvs,
                            variable = c('Age', 'Female', 'LOMG', 'Thymoma', 
                                         'Giant cells/granulomas',  'Corticosteroids',
                                         'Immunomodulators',
                                         'Respiratory support',
                                         'Death'))

base_data |>
  forestplot(labeltext = c(variable, propp, propn, OR, pvs),
             align = c("l", "c", "c", "c", "c"),     # <-- Align columns
             xlim = c(0, 50),
             title = expression(bold("\\large B. Myocarditis")),
             xlab = expression(bold("\\Large Association with myocarditis versus no myocarditis")),
             xlog = TRUE,
             boxsize = .35,
             lwd.ci = 5,                       # Thicker CI lines
             xticks = c(0.25, 0.5, 1, 2, 5, 50),
             zero = 1,       
             lty.zero = 2, 
             col = fpColors(zero = "black"),  # Location of vertical line
             lwd.zero = 4,                     # <--- Thicker vertical line
             txt_gp = fpTxtGp(
               xlab = gpar(fontsize = 12),
               ticks = gpar(fontsize = 12),    # <--- Larger axis tick labels
               label = gpar(fontsize = 12)
             ) # Enlarge axis tick labels
  ) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_header(variable = c("", "Variable\n"),
                propp = c("", "Prevalence \n Myocarditis+"),
                propn = c("", "\n Myocarditis-"),
                OR = c("", "Adj-OR\n"),
                pvs = c("", "\\emph{p}-value\n")) |>
  fp_set_zebra_style("#EFEFEF")
