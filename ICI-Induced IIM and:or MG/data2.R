##################### set working directory
setwd("/Users/jileilin/Desktop/Research/thymomal case reports/")

##################### load packages
library("strucchange")
library('nnet')
library('quantreg')
library('multcomp')
library('DescTools')

#################### read data
data <- read.csv('data2.csv')[, -1]

# adjust contrast level
cancer.detailed <- data$cancer.detailed
ll <- c()
for (i in 1:nrow(data)) {
  if (is.na(cancer.detailed[i]) == TRUE) {
    ll[i] <- NA
  } else if (cancer.detailed[i] == 'thymic') {
    ll[i] <- 1
  } else if (cancer.detailed[i] == 'lung') {
    ll[i] <- 2
  } else if (cancer.detailed[i] == 'melanoma') {
    ll[i] <- 3
  } else if (cancer.detailed[i] == 'gastrointestinal') {
    ll[i] <- 4
  } else if (cancer.detailed[i] == 'genitourinary') {
    ll[i] <- 5
  } else if (cancer.detailed[i] == 'others') {
    ll[i] <- 6
  }
}
cancer.detailed <- reorder(cancer.detailed, ll)
contrasts(cancer.detailed) <- contr.sum(6)
data$cancer.detailed <- cancer.detailed


#################### Compare Dermatomyositis with Myositis (other)
#################### Cancer type
tab <- with(data, t(table(cancer, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:4], ' (C)'),
                   paste0(colnames(tab)[1:4], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 5:8] <- tab[3, 1:4] / sum(tab[3, 1:4])
print(tab)

# likelihood ratio test
subset <- which(is.na(data$derm) == FALSE)
mod1 <- multinom(cancer ~ derm + Sex + Age, data = data, subset = subset)
mod2 <- multinom(cancer ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)

#################### Compare Dermatomyositis with Myositis (other)
#################### Age
Age <- data$Age
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))

tab <- c(round(median(Age[which(data$derm == 'Myositis (other)')], na.rm = TRUE), 2), 
         round(median(Age[which(data$derm == 'Dermatomyositis')], na.rm = TRUE), 2), 
         round(median(Age, na.rm = TRUE), 2),
         iqr(Age[which(data$derm == 'Myositis (other)')]), 
         iqr(Age[which(data$derm == 'Dermatomyositis')]),
         iqr(Age))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Myositis (other)', 'Dermatomyositis', 'All')
print(tab)

# quantile regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
mod <- rq(Age ~ derm + Sex, tau = 0.50, data = data, subset = subset)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
print('coefficient estimate for Dermatomyositis and Myositis (other)')
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2) 
print(coef.ci)

#################### Compare Dermatomyositis with Myositis (other)
#################### Sex
tab <- with(data, t(table(Sex, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
sex.ind <- ifelse(data$Sex == 'F', 1, 0)
mod <- glm(sex.ind ~ derm + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50),
           subset = subset)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Compare Dermatomyositis with Myositis (other)
#################### Treatment
tab <- with(data, t(table(treatment, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:3], ' (C)'),
                   paste0(colnames(tab)[1:3], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 4:6] <- tab[3, 1:3] / sum(tab[3, 1:3])
print(tab)

# likelihood ratio test
subset <- which(is.na(data$derm) == FALSE)
mod1 <- multinom(treatment ~ derm + Sex + Age, subset = subset, data = data)
mod2 <- multinom(treatment ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)


#################### Compare Dermatomyositis with Myositis (other)
#################### Outcome
tab <- with(data, t(table(outcome, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:5], ' (C)'),
                   paste0(colnames(tab)[1:5], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, (1:5) + 5] <- tab[3, 1:5] / sum(tab[3, 1:5])
tab

# likelihood ratio test
mod1 <- multinom(outcome ~ Sex + Age, data = data, subset = which(is.na(data$derm) == FALSE))
mod2 <- multinom(outcome ~ derm + Sex + Age, data = data)
anova(mod1, mod2)

# logistic regression
death.ind <- ifelse(data$outcome == 'Death', 1, 0)
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
mod <- glm(death.ind ~ derm + Sex + Age, data = data, family = 'binomial', subset = subset)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Compare Dermatomyositis with Myositis (other)
#################### MSA
# prepare table
tab <- with(data, t(table(msa, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
msa.ind <- ifelse(data$msa == '+', 1, 0)
mod <- glm(msa.ind ~ derm + Sex + Age, family = 'binomial', data = data, subset = subset)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Compare Dermatomyositis with Myositis (other)
#################### AChR
tab <- with(data, t(table(AChR, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])

# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
achr.ind <- ifelse(data$AChR == '+', 1, 0)
mod <- glm(achr.ind ~ derm + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


#################### Compare Dermatomyositis with Myositis (other)
#################### Striational antibody
tab <- with(data, t(table(Striated, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)


# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
Striated.ind <- ifelse(data$Striated == '+', 1, 0)
mod <- glm(Striated.ind ~ derm + Sex + Age, data = data, subset = subset, family = 'binomial')
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Compare Dermatomyositis with Myositis (other)
#################### O+
tab <- with(data, t(table(Op, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
Op.ind <- ifelse(data$Op == '+', 1, 0)
mod <- glm(Op.ind ~ derm + Sex + Age, family = 'binomial', subset = subset, data = data)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Compare Dermatomyositis with Myositis (other)
#################### Peak Creatine kinase level
pck <- data$pck
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
tab <- c(round(median(pck[which(data$derm == 'Myositis (other)')], na.rm = TRUE), 2),
         round(median(pck[which(data$derm == 'Dermatomyositis')], na.rm = TRUE), 2), 
         round(median(pck, na.rm = TRUE), 2),
         iqr(pck[which(data$derm == 'Myositis (other)')]),
         iqr(pck[which(data$derm == 'Dermatomyositis')]), 
         iqr(pck))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Myositis (other)', 'Dermatomyositis', 'All')
print(tab)

# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
mod <- rq(pck ~ derm + Age + Sex, tau = 0.50, subset = subset, data = data)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2) 
print(coef.ci)


#################### Compare Dermatomyositis with Myositis (other)
#################### Interval
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
interval <- data$interval
tab <- c(round(median(interval[which(data$derm == 'Myositis (other)')], na.rm = TRUE), 2),
         round(median(interval[which(data$derm == 'Dermatomyositis')], na.rm = TRUE), 2), 
         round(median(interval, na.rm = TRUE), 2),
         round(mean(interval[which(data$derm == 'Myositis (other)')], na.rm = TRUE), 2),
         round(mean(interval[which(data$derm == 'Dermatomyositis')], na.rm = TRUE), 2), 
         round(mean(interval, na.rm = TRUE), 2),
         round(sd(interval[which(data$derm == 'Myositis (other)')], na.rm = TRUE), 2),
         round(sd(interval[which(data$derm == 'Dermatomyositis')], na.rm = TRUE), 2), 
         round(sd(interval, na.rm = TRUE), 2),
         iqr(interval[which(data$derm == 'Myositis (other)')]),
         iqr(interval[which(data$derm == 'Dermatomyositis')]), 
         iqr(interval))
tab <- matrix(tab, 4, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'Mean', 'Standard deviation', 'IQR')
colnames(tab) <- c('Myositis (other)', 'Dermatomyositis', 'All')
print(tab)

# median regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
mod <- rq(interval ~ derm + Age + Sex, tau = 0.50, data = data, subset = subset)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# regression coefficient 
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2) 
print(coef.ci)

# ordinary least squares regression
mod <- lm(interval ~ derm + Age + Sex, data = data, subset = subset)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# regression coefficients
coefs <- summary(mod)$coefficients[2, 1]
ses <- summary(mod)$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2) 
print(coef.ci)


#################### Compare Dermatomyositis with Myositis (other)
#################### MG-like syndrome
tab <- with(data, t(table(MGS, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])

# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
MGS.ind <- ifelse(data$MGS == '+', 1, 0)
mod <- glm(MGS.ind ~ derm + Sex + Age, family = 'binomial', subset = subset, data = data)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Compare Dermatomyositis with Myositis (other)
#################### Myocarditis
tab <- with(data, t(table(myocarditis, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
myocarditis.ind <- ifelse(data$myocarditis == '+', 1, 0)
mod <- glm(myocarditis.ind ~ derm + Sex + Age, family = 'binomial', data = data, subset = subset)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Compare Dermatomyositis with Myositis (other)
#################### pre-existing and de novo
tab <- with(data, t(table(pd, derm)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
subset <- which(data$derm %in% c('Dermatomyositis', 'Myositis (other)'))
pd.ind <- ifelse(data$pd == 'P', 1, 0)
mod <- glm(pd.ind ~ derm + Sex + Age, family = 'binomial', data = data, subset = subset)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)




#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Cancer type
tab <- with(data, table(cancer.detailed, myositis.subtypes.other2))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)',
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$cancer.detailed) == FALSE)
mod1 <- multinom(myositis.subtypes.other2 ~ cancer.detailed + Sex + Age,
                 data = data, subset = subset)
mod2 <- multinom(myositis.subtypes.other2 ~ Sex + Age, 
                 data = data, subset = subset)
anova(mod1, mod2)


# paired comparison
# Significance level of 0.05 / 15 should be used
levels <- levels(data$cancer.detailed)
for (i in 1:5) {
  for (j in (i + 1):6) {
    print(paste0('Comparison between ', levels[i], ' and ', levels[j]))
    subset2 <- with(data, which(is.na(cancer.detailed) == FALSE &
                                 cancer.detailed %in% c(levels[i], levels[j])))
    mod1 <- multinom(myositis.subtypes.other2 ~ cancer.detailed + Sex + Age, 
                     subset = subset2, data = data)
    mod2 <- multinom(myositis.subtypes.other2 ~ Sex + Age, 
                     subset = subset2, data = data)
    mmmm <- anova(mod1, mod2)
    p <- mmmm$`Pr(Chi)`[2]
    print(p)
  }
}


# prepare table
tab <- with(data, t(table(cancer.detailed, myositis.subtypes.other2)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:6], ' (C)'),
                   paste0(colnames(tab)[1:6], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, (1:6) + 6] <- tab[4, 1:6] / sum(tab[4, 1:6])
print(tab)

# likelihood ratio test
subset3 <- which(is.na(data$myositis.subtypes.other2) == FALSE)
mod1 <- multinom(cancer.detailed ~ Sex + Age, 
                 subset = subset3, data = data)
mod2 <- multinom(cancer.detailed ~ myositis.subtypes.other2 + Sex + Age,
                 subset = subset3, data = data)
anova(mod1, mod2)

# paired comparison
# Significance level of 0.05 / 3 should be used
levels <- levels(as.factor(data$myositis.subtypes.other2))
for (i in 1:2) {
  for (j in (i + 1):3) {
    
    print(paste0('Comparison between ', levels[i], ' and ', levels[j]))
    subset4 <- with(data, which(is.na(myositis.subtypes.other2) == FALSE &
                                  myositis.subtypes.other2 %in% c(levels[i], levels[j])))
    mod1 <- multinom(cancer.detailed ~ myositis.subtypes.other2 + Sex + Age, 
                     subset = subset4, data = data)
    mod2 <- multinom(cancer.detailed ~ Sex + Age, 
                     subset = subset4, data = data)
    mmmm <- anova(mod1, mod2)
    p <- mmmm$`Pr(Chi)`[2]
    print(p)
  }
}

#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Age
Age <- data$Age
m3 <- data$myositis.subtypes.other2
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
tab <- c(round(median(Age[which(m3 == 'Myocarditis (only C)')], na.rm = TRUE), 2), 
         round(median(Age[which(m3 == 'Myositis (no D and no C)')], na.rm = TRUE), 2),
         round(median(Age[which(m3 == 'Myositis with myocarditis')], na.rm = TRUE), 2),
         round(median(Age[which(is.na(m3) == FALSE)], na.rm = TRUE), 2),
         iqr(Age[which(m3 == 'Myocarditis (only C)')]), 
         iqr(Age[which(m3 == 'Myositis (no D and no C)')]),
         iqr(Age[which(m3 == 'Myositis with myocarditis')]),
         iqr(Age[which(is.na(m3) == FALSE)]))
tab <- matrix(tab, 2, 4, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Myocarditis (only C)', 'Myositis (no D and no C)', 'Myositis with myocarditis', 'All')
print(tab)

# median regression
subset <- which(m3 %in% c('Myocarditis (only C)', 'Myositis (no D and no C)'))
mod <- rq(Age ~ m3 + Sex, tau = 0.50, subset = subset, data = data)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimates
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# compare Age between Myocarditis (only C) and Myositis with myocarditis
subset2 <- which(m3 %in% c('Myocarditis (only C)', 'Myositis with myocarditis'))
mod <- rq(Age ~ m3 + Sex, tau = 0.50, subset = subset2, data = data)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimates
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)


# Compare Age between Myositis (no D and no C) and Myositis with myocarditis
subset3 <- which(m3 %in% c('Myositis (no D and no C)', 'Myositis with myocarditis'))
mod <- rq(Age ~ m3 + Sex, tau = 0.50, subset = subset3, data = data)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)


#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Sex
tab <- with(data, table(Sex, myositis.subtypes.other2))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)',
                                                 ' (P)', ' (P)', ' (P)'))

# likelihood ratio test
subset <- which(is.na(data$Sex) == FALSE)
mod1 <- multinom(myositis.subtypes.other2 ~ Age, subset = subset, data = data)
mod2 <- multinom(myositis.subtypes.other2 ~ Sex + Age, data = data)
anova(mod1, mod2)

# prepare table
tab <- with(data, t(table(Sex, myositis.subtypes.other2)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:4])
print(tab)

# multiple comparison
m3 <- as.factor(data$myositis.subtypes.other2)
sex.ind <- ifelse(data$Sex == 'F', 1, 0)
mod <- glm(sex.ind ~ m3 + Age, family = 'binomial', data = data)
mod.glht <- glht(mod, linfct = mcp("m3" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)

# p-value
print(round(pvs, 4))


#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Treatment
tab <- with(data, table(treatment, myositis.subtypes.other2))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)',
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$treatment) == FALSE)
mod1 <- multinom(myositis.subtypes.other2 ~ treatment + Sex + Age, data = data, subset = subset)
mod2 <- multinom(myositis.subtypes.other2 ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)

# prepare table
tab <- with(data, t(table(treatment, myositis.subtypes.other2)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:3], ' (C)'),
                   paste0(colnames(tab)[1:3], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 4:6] <- tab[4, 1:3] / sum(tab[4, 1:3])
print(tab)

# likelihood ratio test
subset2 <- which(is.na(data$myositis.subtypes.other2) == FALSE)
mod1 <- multinom(treatment ~ Sex + Age, subset = subset2, data = data)
mod2 <- multinom(treatment ~ myositis.subtypes.other2 + Sex + Age, subset = subset2, data = data)
anova(mod1, mod2)

#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Outcome
death <- ifelse(data$outcome == 'Death', '+', '-')
death.ind <- ifelse(data$outcome == 'Death', 1, 0)
tab <- with(data, table(outcome, myositis.subtypes.other2))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)',
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$outcome) == FALSE)
mod1 <- multinom(myositis.subtypes.other2 ~ outcome + Sex + Age, subset = subset, data = data)
mod2 <- multinom(myositis.subtypes.other2 ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)

# prepare table
tab <- with(data, t(table(outcome, myositis.subtypes.other2)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:5], ' (C)'),
                   paste0(colnames(tab)[1:5], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 6:10] <- tab[4, 1:5] / sum(tab[4, 1:5])
print(tab)


# likelihood ratio test
subset2 <- which(is.na(data$myositis.subtypes.other2) == FALSE)
mod1 <- multinom(outcome ~ Sex + Age, subset = subset2, data = data)
mod2 <- multinom(outcome ~ myositis.subtypes.other2 + Sex + Age, data = data, subset = subset2)
anova(mod1, mod2)


# paired comparison of outcome
# Significance level of 0.05 / 3 should be used instead
levels <- levels(as.factor(data$myositis.subtypes.other2))
for (i in 1:2) {
  for (j in (i + 1):3) {
    subset3 <- with(data, which(myositis.subtypes.other2 %in% c(levels[i], levels[j]) & 
                                  is.na(myositis.subtypes.other2) == FALSE))
    index <- which(m3 %in% c(levels[i], levels[j]))
    mod1 <- multinom(outcome ~ myositis.subtypes.other2 + Sex + Age, subset = subset3, data = data)
    mod2 <- multinom(outcome ~ Sex + Age, subset = subset3, data = data)
    mmmm <- anova(mod1, mod2)
    p <- mmmm$`Pr(Chi)`[2]
    print(p)
  }
}

# paired comparison
m3 <- as.factor(data$myositis.subtypes.other2)
mod <- glm(death.ind ~ m3 + Sex + Age, family = 'binomial', data = data)
mod.glht <- glht(mod, linfct = mcp("m3" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)

# p-value
print(round(pvs, 4))

# odds ratio for paired comparison
cis <- confint(mod.glht)$confint
cis[1, ] <- -cis[1, c(1, 3, 2)]
rownames(cis)[1] <- 'Myocarditis (only C) - Myositis (no D and no C)'
odds <- exp(cis)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Death and AChR within myocarditis
death <- ifelse(data$outcome == 'Death', '+', '-')
death <- ifelse(death == '+', 1, 0)
subset <- with(data, which(is.na(myositis.subtypes.other2) == FALSE & myocarditis == '+'))
tab <- with(data, table(AChR[subset], death[subset]))
tab <- cbind(tab, tab / rowSums(tab)) 
rownames(tab) <- c('AChR (-)', 'AChR (+)')
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))

# logistic regression
mod <- glm(death.ind ~ AChR + Sex + Age, family = 'binomial', 
           subset = subset, data = data, control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(AChR[subset], death[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
rownames(tab) <- c('Alive', 'Dead')
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])

# logistic regression
achr.ind <- ifelse(data$AChR == '+', 1, 0)
mod <- glm(achr.ind ~ death + Sex + Age, family = 'binomial', 
           subset = subset, data = data, control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Death and AChR within myositis (+) and myocarditis (-) and myositis (+) and myocarditis (+)
subset <- with(data, which(is.na(myositis.subtypes.other2) == FALSE & inflam.myopathy == '+'))
death <- ifelse(data$outcome == 'Death', '+', '-')
tab <- with(data, table(AChR[subset], death[subset]))
tab <- cbind(tab, tab / rowSums(tab)) 
rownames(tab) <- c('AChR (-)', 'AChR (+)')
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))

# logistic regression
death.ind <- ifelse(data$outcome == 'Death', 1, 0)
mod <- glm(death.ind ~ AChR + Sex + Age, family = 'binomial', 
           subset = subset, data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# logistic regression
# interaction effects
mod <- glm(death.ind ~ AChR * myocarditis + Sex + Age, family = 'binomial', 
           subset = subset, data = data, 
           control = glm.control(maxit = 50))

# p-value: interaction effects
print(coef(summary(mod))[6, 4])


# p-value: main effects of AChR
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


# p-value: main effects of myocarditis
print(coef(summary(mod))[3, 4])
# odds ratio
odds <- coef(summary(mod))[3, 1]
ses <- coef(summary(mod))[3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


# prepare table
tab <- with(data, t(table(AChR[subset], death[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
rownames(tab) <- c('Alive', 'Dead')
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
achr.ind <- ifelse(data$AChR == '+', 1, 0)
mod <- glm(achr.ind ~ death + Sex + Age, family = 'binomial', 
           subset = subset, data = data,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### MSA
tab <- with(data, table(msa, myositis.subtypes.other2))
tab <- cbind(tab, tab / sum(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$msa) == FALSE)
mod1 <- multinom(myositis.subtypes.other2 ~ msa + Sex + Age, subset = subset, data = data)
mod2 <- multinom(myositis.subtypes.other2 ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)

# prepare table
tab <- with(data, t(table(msa, myositis.subtypes.other2)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])

# logistic regression
m3 <- as.factor(data$myositis.subtypes.other2)
msa.ind <- ifelse(data$msa == '+', 1, 0)
mod <- glm(msa.ind ~ m3 + Age, data = data, family = 'binomial')

# multiple comparison
mod.glht <- glht(mod, linfct = mcp("m3" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)

# p-value
print(round(pvs, 4))

#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### AChR
tab <- with(data, t(table(AChR, myositis.subtypes.other2)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# Distribution of O+ within AChR tested and myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
test.AChR <- ifelse(is.na(data$AChR) == TRUE, 0, 1)
tab <- with(data, table(test.AChR, Op, myositis.subtypes.other2))
tab.pre <- c()
for (j in 1:3) {
  tab.pre[[j]] <- round(cbind(tab[, , j], tab[, , j] / rowSums(tab[, , j])), 2)
  colnames(tab.pre[[j]]) <- c(paste0(c('O (-)', 'O (+)'), ' (C)'), paste0(c('O (-)', 'O (+)'), ' (P)'))
  rownames(tab.pre[[j]]) <- c('Untested', 'Tested')
}
names(tab.pre) <- names(table(data$myositis.subtypes.other2))
print(tab.pre)

# Compare O+ between tested and untested for AChR
# combining odds ratio and p-values
pp <- c()
oddss <- c()

# logistic regression
Op.ind <- ifelse(data$Op == '+', 1, 0)
subset <- which(data$myositis.subtypes.other2 == 'Myocarditis (only C)')
mod <- glm(Op.ind ~ test.AChR + Sex + Age, data = data, 
           subset = subset, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
oddss <- rbind(oddss, odds)
pp <- c(pp, p)

# Compare O+ between tested and untested for Myositis (no D and no C)
# logistic regression
subset2 <- which(data$myositis.subtypes.other2 == 'Myositis (no D and no C)')
mod <- glm(Op.ind ~ test.AChR + Sex + Age, data = data, 
           subset = subset2, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
oddss <- rbind(oddss, odds)
pp <- c(pp, p)

# Compare O+ between tested and untested for Myositis with myocarditis
subset3 <- which(data$myositis.subtypes.other2 == 'Myositis with myocarditis')
mod <- glm(Op.ind ~ test.AChR + Sex + Age, data = data, 
           subset = subset3, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
oddss <- rbind(oddss, odds)
pp <- c(pp, p)


# Compare the distribution of AChR between O+ and O-
subset4 <- which(is.na(data$myositis.subtypes.other2) == FALSE)
tab <- with(data, t(table(AChR[subset4], Op[subset4])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
mod <- glm(Op.ind ~ AChR + Sex + Age, data = data, subset = subset4, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
oddss <- rbind(oddss, odds)
pp <- c(pp, p)
print(oddss)


#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### AChR level
AChR.level <- data$AChR.level
m3 <- data$myositis.subtypes.other2
tab <- c(round(median(AChR.level[which(m3 == 'Myocarditis (only C)')], na.rm = TRUE), 2), 
         round(median(AChR.level[which(m3 == 'Myositis (no D and no C)')], na.rm = TRUE), 2),
         round(median(AChR.level[which(m3 == 'Myositis with myocarditis')], na.rm = TRUE), 2),
         round(median(AChR.level[which(is.na(m3) == FALSE)], na.rm = TRUE), 2),
         iqr(AChR.level[which(m3 == 'Myocarditis (only C)')]), 
         iqr(AChR.level[which(m3 == 'Myositis (no D and no C)')]),
         iqr(AChR.level[which(m3 == 'Myositis with myocarditis')]),
         iqr(AChR.level[which(is.na(m3) == FALSE)]))
tab <- matrix(tab, 2, 4, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Myocarditis (only C)', 'Myositis (no D and no C)', 'Myositis with myocarditis', 'All')
print(tab)

# compare Myocarditis (only C) with Myositis (no D and no C)
subset <- which(m3 %in% c('Myocarditis (only C)', 'Myositis (no D and no C)'))
mod <- rq(AChR.level ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate for Myocarditis (only C) and Myositis (no D and no C)
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# Compare AChR level between Myocarditis (only C) and Myositis with myocarditis
subset <- which(m3 %in% c('Myocarditis (only C)', 'Myositis with myocarditis'))
mod <- rq(AChR.level ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

print('Compare AChR level between Myositis (no D and no C) and Myositis with myocarditis')

# Compare AChR level betweenMyositis (no D and no C) and Myositis with myocarditis
subset <- which(m3 %in% c('Myositis (no D and no C)', 'Myositis with myocarditis'))
mod <- rq(AChR.level ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# Coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)


#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Striational antibody
Striated.ind <- ifelse(data$Striated == '+', 1, 0)
m3 <- as.factor(data$myositis.subtypes.other2)
mod <- glm(Striated.ind ~ m3 + Age, family = 'binomial', data = data)

# paired comparison
mod.glht <- glht(mod, linfct = mcp("m3" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

# test for strational antibody
test.Striated <- ifelse(is.na(data$Striated) == TRUE, 0, 1)
# Distribution of O+ within Striational antibody tested and myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
tab <- with(data, table(test.Striated, Op, m3))
tab.pre <- c()
for (j in 1:3) {
  tab.pre[[j]] <- round(cbind(tab[, , j], tab[, , j] / rowSums(tab[, , j])), 2)
  colnames(tab.pre[[j]]) <- c(paste0(c('O (-)', 'O (+)'), ' (C)'), paste0(c('O (-)', 'O (+)'), ' (P)'))
  rownames(tab.pre[[j]]) <- c('Untested', 'Tested')
}
names(tab.pre) <- names(table(m3))
print(tab.pre)

# Compare the distribution of striational antibody between O+ and O-
subset <- which(is.na(m3) == FALSE)
tab <- with(data, t(table(Striated[subset], Op[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
mod <- glm(Op.ind ~ Striated + Sex + Age, data = data, subset = subset, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Peak creatine kinase level
pck <- data$pck
m3 <- data$myositis.subtypes.other2
subset <- which(is.na(m3) == FALSE)
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
tab <- c(round(median(pck[which(m3 == 'Myocarditis (only C)')], na.rm = TRUE), 2), 
         round(median(pck[which(m3 == 'Myositis (no D and no C)')], na.rm = TRUE), 2),
         round(median(pck[which(m3 == 'Myositis with myocarditis')], na.rm = TRUE), 2),
         round(median(pck[subset], na.rm = TRUE), 2),
         iqr(pck[which(m3 == 'Myocarditis (only C)')]), 
         iqr(pck[which(m3 == 'Myositis (no D and no C)')]),
         iqr(pck[which(m3 == 'Myositis with myocarditis')]),
         iqr(pck[subset]))
tab <- matrix(tab, 2, 4, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Myocarditis (only C)', 'Myositis (no D and no C)', 'Myositis with myocarditis', 'All')
print(tab)

# Compare Peak creatine kinase level between Myocarditis (only C) and Myositis (no D and no C)
subset <- which(m3 %in% c('Myocarditis (only C)', 'Myositis (no D and no C)'))
mod <- rq(pck ~ m3 + Sex + Age, data = data, tau = 0.50, subset = index)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4] * 3
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(1 - 0.05 / 6) * ses, coefs + qnorm(1 - 0.05 / 6) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# Compare Peak creatine kinase level between Myocarditis (only C) and Myositis with myocarditis
subset <- which(m3 %in% c('Myocarditis (only C)', 'Myositis with myocarditis'))
mod <- rq(pck ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4] * 3
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(1 - 0.05 / 6) * ses, coefs + qnorm(1 - 0.05 / 6) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# Compare Peak creatine kinase level between Myositis (no D and no C) and Myositis with myocarditis
subset <- which(m3 %in% c('Myositis (no D and no C)', 'Myositis with myocarditis'))
mod <- rq(pck ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4] * 3
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(1 - 0.05 / 6) * ses, coefs + qnorm(1 - 0.05 / 6) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)


#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Peak troponin
pt <- data$pt
m3 <- data$myositis.subtypes.other2
subset <- which(is.na(m3) == FALSE)
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))

tab <- c(round(median(pt[which(m3 == 'Myocarditis (only C)')], na.rm = TRUE), 2), 
         round(median(pt[which(m3 == 'Myositis (no D and no C)')], na.rm = TRUE), 2),
         round(median(pt[which(m3 == 'Myositis with myocarditis')], na.rm = TRUE), 2),
         round(median(pt[subset], na.rm = TRUE), 2),
         iqr(pt[which(m3 == 'Myocarditis (only C)')]), 
         iqr(pt[which(m3 == 'Myositis (no D and no C)')]),
         iqr(pt[which(m3 == 'Myositis with myocarditis')]),
         iqr(pt[subset]))
tab <- matrix(tab, 2, 4, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Myocarditis (only C)', 'Myositis (no D and no C)', 'Myositis with myocarditis', 'All')
print(tab)

# Compare Peak troponin between Myocarditis (only C) and Myositis (no D and no C)
subset <- which(m3 %in% c('Myocarditis (only C)', 'Myositis (no D and no C)'))
mod <- rq(pt ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4] * 3
print(p)

# Coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(1 - 0.05 / 6) * ses, coefs + qnorm(1 - 0.05 / 6) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# Compare Peak troponin between Myocarditis (only C) and Myositis with myocarditis
subset <- which(m3 %in% c('Myocarditis (only C)', 'Myositis with myocarditis'))
mod <- rq(pt ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(1 - 0.05 / 6) * ses, coefs + qnorm(1 - 0.05 / 6) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# Compare Peak troponin between Myositis (no D and no C) and Myositis with myocarditis
subset <- which(m3 %in% c('Myositis (no D and no C)', 'Myositis with myocarditis'))
mod <- rq(pt ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(1 - 0.05 / 6) * ses, coefs + qnorm(1 - 0.05 / 6) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)


#################### Compare myositis (no D and no C), myocarditis (only C), and myositis with myocarditis
#################### Interval
m3 <- data$myositis.subtypes.other2
interval <- data$interval
subset <- which(is.na(m3) == FALSE)
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
tab <- c(round(median(interval[which(m3 == 'Myocarditis (only C)')], na.rm = TRUE), 2), 
         round(median(interval[which(m3 == 'Myositis (no D and no C)')], na.rm = TRUE), 2),
         round(median(interval[which(m3 == 'Myositis with myocarditis')], na.rm = TRUE), 2),
         round(median(interval[subset], na.rm = TRUE), 2),
         round(mean(interval[which(m3 == 'Myocarditis (only C)')], na.rm = TRUE), 2), 
         round(mean(interval[which(m3 == 'Myositis (no D and no C)')], na.rm = TRUE), 2),
         round(mean(interval[which(m3 == 'Myositis with myocarditis')], na.rm = TRUE), 2),
         round(mean(interval[subset], na.rm = TRUE), 2),
         round(sd(interval[which(m3 == 'Myocarditis (only C)')], na.rm = TRUE), 2), 
         round(sd(interval[which(m3 == 'Myositis (no D and no C)')], na.rm = TRUE), 2),
         round(sd(interval[which(m3 == 'Myositis with myocarditis')], na.rm = TRUE), 2),
         round(sd(interval[subset], na.rm = TRUE), 2),
         iqr(interval[which(m3 == 'Myocarditis (only C)')]), 
         iqr(interval[which(m3 == 'Myositis (no D and no C)')]),
         iqr(interval[which(m3 == 'Myositis with myocarditis')]),
         iqr(interval[subset]))
tab <- matrix(tab, 4, 4, byrow = TRUE)
rownames(tab) <- c('Median', 'Mean', 'STD', 'IQR')
colnames(tab) <- c('Myocarditis (only C)', 'Myositis (no D and no C)', 'Myositis with myocarditis', 'All')
print(tab)

# Compare interval between Myocarditis (only C) and Myositis (no D and no C) (median)
index <- which(m3 %in% c('Myocarditis (only C)', 'Myositis (no D and no C)'))
mod <- rq(interval ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# Compare interval between Myocarditis (only C) and Myositis with myocarditis (median)
subset <- which(m3 %in% c('Myocarditis (only C)', 'Myositis with myocarditis'))
mod <- rq(interval ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# Compare interval between Myositis (no D and no C) and Myositis with myocarditis (median)
subset <- which(m3 %in% c('Myositis (no D and no C)', 'Myositis with myocarditis'))
mod <- rq(interval ~ m3 + Sex + Age, tau = 0.50, subset = subset, data = data)
# p-value (Significance level of 0.05 / 3 should be used instead) 
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

# analysis of variance
# Compare interval between Myositis (no D and no C) and Myositis with myocarditis (mean)
mod <- aov(interval ~ m3 + Sex + Age, data = data)

# p-value
# tukey paired comparison
mmmm <- TukeyHSD(mod, which = "m3")
round(mmmm$m3[, 4], 4)



#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Giant Cells and Granuloma
tab <- with(data, table(GG, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

tab <- with(data, t(table(GG, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Inflammatory cells
tab <- with(data, table(IC, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$IC) == FALSE)
mod1 <- multinom(m3 ~ IC + Sex + Age, subset = subset, data = data)
mod2 <- multinom(m3 ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)

# prepare table
tab <- with(data, t(table(IC, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
mmm <- as.factor(data$muscle.biopsy)
contrasts(mmm) <- contr.sum(3)
IC.ind <- ifelse(data$IC == '+', 1, 0)
mod <- glm(IC.ind ~ mmm + Sex + Age, data = data, family = 'binomial')
# paired comparison
mod.glht <- glht(mod, linfct = mcp("mmm" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Em
subset <- which(data$muscle.biopsy != 'Group 2')
# prepare table
tab <- with(data, table(Emp[subset], muscle.biopsy[subset]))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', 
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mmm.ind <- ifelse(data$muscle.biopsy == 'Group 1', 1, 0)
mod <- glm(mmm.ind ~ Emp + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound') 

# prepare table
tab <- with(data, t(table(Emp[subset], muscle.biopsy[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
Em.ind <- ifelse(data$Emp == '+', 1, 0)
mod <- glm(Em.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4] 
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 3 vs Group 1')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Pm+
subset <- which(data$muscle.biopsy != 'Group 2')
# prepare table
tab <- with(data, table(Pmp[subset], muscle.biopsy[subset]))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', 
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mmm.ind <- ifelse(data$muscle.biopsy == 'Group 1', 1, 0)
mod <- glm(mmm.ind ~ Pmp + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound') 

# prepare table
tab <- with(data, t(table(Pmp[subset], muscle.biopsy[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
Pm.ind <- ifelse(data$Pmp == '+', 1, 0)
mod <- glm(Pm.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4] 
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 3 vs Group 1')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Is+
subset <- which(data$muscle.biopsy != 'Group 2')
# prepare table
tab <- with(data, table(Isp[subset], muscle.biopsy[subset]))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', 
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mmm.ind <- ifelse(data$muscle.biopsy == 'Group 1', 1, 0)
mod <- glm(mmm.ind ~ Pmp + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound') 

# prepare table
tab <- with(data, t(table(Isp[subset], muscle.biopsy[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
Is.ind <- ifelse(data$Isp == '+', 1, 0)
mod <- glm(Is.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4] 
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 3 vs Group 1')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Pv+
subset <- which(data$muscle.biopsy != 'Group 2')
# prepare table
tab <- with(data, table(Pvp[subset], muscle.biopsy[subset]))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', 
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mmm.ind <- ifelse(data$muscle.biopsy == 'Group 1', 1, 0)
mod <- glm(mmm.ind ~ Pvp + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound') 

# prepare table
tab <- with(data, t(table(Pvp[subset], muscle.biopsy[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
Pv.ind <- ifelse(data$Pvp == '+', 1, 0)
mod <- glm(Pv.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4] 
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 3 vs Group 1')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### CD4
tab <- with(data, table(CD4, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))  
print(tab)

# likelihood ratio test
subset <- which(is.na(data$CD4) == FALSE) 
mod1 <- multinom(muscle.biopsy ~ CD4 + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)


# prepare table
tab <- with(data, t(table(CD4, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2]) 
print(tab)

# logistic regression
CD4.ind <- ifelse(data$CD4 == '+', 1, 0)
data$muscle.biopsy <- as.factor(data$muscle.biopsy)
mod <- glm(CD4.ind ~ muscle.biopsy + Sex + Age, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')

# pairwise comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### CD8
tab <- with(data, table(CD8, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))  
print(tab)

# likelihood ratio test
subset <- which(is.na(data$CD8) == FALSE) 
mod1 <- multinom(muscle.biopsy ~ CD8 + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)


# prepare table
tab <- with(data, t(table(CD8, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2]) 
print(tab)

# logistic regression
CD8.ind <- ifelse(data$CD8 == '+', 1, 0)
data$muscle.biopsy <- as.factor(data$muscle.biopsy)
mod <- glm(CD8.ind ~ muscle.biopsy + Sex + Age, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')

# pairwise comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### CD20B
tab <- with(data, table(CDB20, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))  
print(tab)

# likelihood ratio test
subset <- which(is.na(data$CDB20) == FALSE) 
mod1 <- multinom(muscle.biopsy ~ CDB20 + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)


# prepare table
tab <- with(data, t(table(CDB20, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2]) 
print(tab)

# logistic regression
CDB20.ind <- ifelse(data$CDB20 == '+', 1, 0)
data$muscle.biopsy <- as.factor(data$muscle.biopsy)
mod <- glm(CDB20.ind ~ muscle.biopsy + Sex + Age, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')

# pairwise comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### CD8>=CD4
tab <- with(data, table(CD.comp, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))  
print(tab)

# likelihood ratio test
subset <- which(is.na(data$CD.comp) == FALSE) 
mod1 <- multinom(muscle.biopsy ~ CD.comp + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)


# prepare table
tab <- with(data, t(table(CD.comp, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2]) 
print(tab)

# logistic regression
CD.comp.ind <- ifelse(data$CD.comp == '+', 1, 0)
data$muscle.biopsy <- as.factor(data$muscle.biopsy)
mod <- glm(CD.comp.ind ~ muscle.biopsy + Sex + Age, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')

# pairwise comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Macrophages
tab <- with(data, table(Macrophages, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$Macrophages) == FALSE)
mod1 <- multinom(muscle.biopsy ~ Macrophages + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)



# prepare table
tab <- with(data, t(table(Macrophages, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
Mp.ind <- ifelse(data$Macrophages == '+', 1, 0)
mod <- glm(Mp.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Plasma
tab <- with(data, table(plasma, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)
tab <- with(data, t(table(plasma, muscle.biopsy)))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
print(tab)

#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### MHC-I
tab <- with(data, table(mhc1, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)
tab <- with(data, t(table(mhc1, muscle.biopsy)))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
print(tab)

#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### MHC-II
tab <- with(data, table(mhc2, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)
tab <- with(data, t(table(mhc2, muscle.biopsy)))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
print(tab)





#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### IgG
tab <- with(data, table(igg, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

tab <- with(data, t(table(igg, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Complement
tab <- with(data, table(complement, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$complement) == FALSE)
mod1 <- multinom(muscle.biopsy ~ complement + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)



# prepare table
tab <- with(data, t(table(complement, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
complement.ind <- ifelse(data$complement == '+', 1, 0)
mod <- glm(complement.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(odds)



#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Focal Degeneration & Necrosis
tab <- with(data, table(FN, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))  
print(tab)

# likelihood ratio test
subset <- which(is.na(data$FN) == FALSE) 
mod1 <- multinom(muscle.biopsy ~ FN + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)


# prepare table
tab <- with(data, t(table(FN, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2]) 
print(tab)

# logistic regression
FN.ind <- ifelse(data$FN == '+', 1, 0)
data$muscle.biopsy <- as.factor(data$muscle.biopsy)
mod <- glm(FN.ind ~ muscle.biopsy + Sex + Age, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')

# pairwise comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Regeneration
tab <- with(data, table(Regeneration, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)')) 
print(tab)

# likelihood ratio test
subset <- which(is.na(data$Regeneration) == FALSE) 
mod1 <- multinom(muscle.biopsy ~ Regeneration + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)


# prepare table
tab <- with(data, t(table(Regeneration, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2]) 
print(tab)

# logistic regression
Regeneration.ind <- ifelse(data$Regeneration == '+', 1, 0)
data$muscle.biopsy <- as.factor(data$muscle.biopsy)
mod <- glm(Regeneration.ind ~ muscle.biopsy + Sex + Age, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')

# pairwise comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Increased Fiber Size Variability
tab <- with(data, table(IFSV, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))  
print(tab)

# likelihood ratio test
subset <- which(is.na(data$IFSV) == FALSE) 
mod1 <- multinom(muscle.biopsy ~ IFSV + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)


# prepare table
tab <- with(data, t(table(IFSV, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2]) 
print(tab)

# logistic regression
IFSV.ind <- ifelse(data$IFSV == '+', 1, 0)
data$muscle.biopsy <- as.factor(data$muscle.biopsy)
mod <- glm(IFSV.ind ~ muscle.biopsy + Sex + Age, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')

# pairwise comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))



#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Increased Connective Tissue
tab <- with(data, table(ICT, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))  
print(tab)

# likelihood ratio test
subset <- which(is.na(data$ICT) == FALSE) 
mod1 <- multinom(muscle.biopsy ~ ICT + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, data = data, subset = subset)
anova(mod1, mod2)


# prepare table
tab <- with(data, t(table(ICT, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2]) 
print(tab)

# logistic regression
ICT.ind <- ifelse(data$ICT == '+', 1, 0)
data$muscle.biopsy <- as.factor(data$muscle.biopsy)
mod <- glm(ICT.ind ~ muscle.biopsy + Sex + Age, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')

# pairwise comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### PD-L1
tab <- with(data, table(pdl1, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$pdl1) == FALSE)
mod1 <- multinom(muscle.biopsy ~ pdl1 + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)



# prepare table
tab <- with(data, t(table(pdl1, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
pdl1.ind <- ifelse(data$pdl1 == '+', 1, 0)
mod <- glm(pdl1.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(odds)



#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### O+
tab <- with(data, table(Op, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$Op) == FALSE)
mod1 <- multinom(muscle.biopsy ~ Op + Sex + Age, data = data, subset = subset)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)



# prepare table
tab <- with(data, t(table(Op, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
Op.ind <- ifelse(data$Op == '+', 1, 0)
mod <- glm(Op.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(odds)



#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### MG
tab <- with(data, table(MG, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$MG) == FALSE)
mod1 <- multinom(muscle.biopsy ~ MG + Sex + Age, data = data)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)

# prepare table
tab <- with(data, t(table(MG, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])


# logistic regression
mg.ind <- ifelse(data$MG == '+', 1, 0)
mod <- glm(mg.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(odds)

# paired comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

# odds ratio
cis <- confint(mod.glht)$confint
cis[3, ] <- - cis[3, ]
cis[3, ] <- cis[3, c(1, 3, 2)]
rownames(cis)[3] <- 'Group 2 - Group 3'
odds <- exp(cis)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### MG-like syndrome
tab <- with(data, table(MGS, muscle.biopsy))
tab <- cbind(tab, tab / rowSums(tab))
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset <- which(is.na(data$MGS) == FALSE)
mod1 <- multinom(muscle.biopsy ~ MGS + Sex + Age, data = data)
mod2 <- multinom(muscle.biopsy ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)

# prepare table
tab <- with(data, t(table(MGS, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
MGS.ind <- ifelse(data$MGS == '+', 1, 0)
mod <- glm(MGS.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2:3, 4]
names(p) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(p)

# odds ratio
odds <- coef(summary(mod))[2:3, 1]
ses <- coef(summary(mod))[2:3, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('Group 1 vs overall', 'Group 2 vs overall')
print(odds)

# paired comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

# odds ratio
cis <- confint(mod.glht)$confint
odds <- exp(cis)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Corticosteriods
tab <- with(data, t(table(cort, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
cort.ind <- ifelse(data$cort == '+', 1, 0)
mod <- glm(cort.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)
# paired comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

# odds ratio
cis <- confint(mod.glht)$confint
odds <- exp(cis)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Plasmapheresis
tab <- with(data, t(table(plasmapheresis, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
plasmapheresis.ind <- ifelse(data$plasmapheresis == '+', 1, 0)
mod <- glm(plasmapheresis.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)
# paired comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

# odds ratio
cis <- confint(mod.glht)$confint
odds <- exp(cis)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Respiratory support
tab <- with(data, t(table(rs, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
rs.ind <- ifelse(data$rs == '+', 1, 0)
mod <- glm(rs.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)
# paired comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

# odds ratio
cis <- confint(mod.glht)$confint
odds <- exp(cis)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Cardiac treatment
tab <- with(data, t(table(ct, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
ct.ind <- ifelse(data$ct == '+', 1, 0)
mod <- glm(ct.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)
# paired comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

# odds ratio
cis <- confint(mod.glht)$confint
odds <- exp(cis)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Muscle biopsy (myositis alone vs myocarditis alone vs myositis & myocarditis vs all)
#################### Non-steroid immunomodulators
tab <- with(data, t(table(imm, muscle.biopsy)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
imm.ind <- ifelse(data$imm == '+', 1, 0)
mod <- glm(imm.ind ~ muscle.biopsy + Sex + Age, family = 'binomial', data = data)
# paired comparison
mod.glht <- glht(mod, linfct = mcp("muscle.biopsy" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)
# p-value
print(round(pvs, 4))

# odds ratio
cis <- confint(mod.glht)$confint
odds <- exp(cis)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Cancer types
tab <- with(data, table(cancer.detailed, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ cancer.detailed + Sex + Age, family = 'binomial', 
           data = data, control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2:6, 4])

# odds ratio
odds <- coef(summary(mod))[2:6, 1]
ses <- coef(summary(mod))[2:6, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
rownames(odds) <- c('thymic vs overall', 'lung vs overall', 'melanoma vs overall',                    
                    'gastrointestinal vs overall', 'genitourinary vs overall')
print(odds)

# prepare table
tab <- with(data, t(table(cancer.detailed, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:6], ' (C)'),
                   paste0(colnames(tab)[1:6], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 7:12] <- tab[3, 1:6] / sum(tab[3, 1:6])
print(tab)

# likelihood ratio test
subset <- which(is.na(data$mg.type) == FALSE)
mod1 <- multinom(cancer.detailed ~ mg.type + Sex + Age, subset = subset, data = data)
mod2 <- multinom(cancer.detailed ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)

#################### MG alone vs MG with myositis/myocarditis
#################### O+
tab <- with(data, table(Op, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ Op + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(Op, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
Op.ind <- ifelse(data$Op == '+', 1, 0)
mod <- glm(Op.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### RNS
tab <- with(data, table(rns, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ rns + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(rns, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
rns.ind <- ifelse(data$rns == '+', 1, 0)
mod <- glm(rns.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)



#################### MG alone vs MG with myositis/myocarditis
#################### SFEMG
tab <- with(data, table(sfemg, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ sfemg + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(sfemg, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
sfemg.ind <- ifelse(data$sfemg == '+', 1, 0)
mod <- glm(sfemg.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### mEMG
tab <- with(data, table(mEMG, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ mEMG + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(mEMG, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
mEMG.ind <- ifelse(data$mEMG == '+', 1, 0)
mod <- glm(mEMG.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### achl
tab <- with(data, table(achl, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ achl + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(achl, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
achl.ind <- ifelse(data$achl == '+', 1, 0)
mod <- glm(achl.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### IP
tab <- with(data, table(ip, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ ip + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(ip, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
ip.ind <- ifelse(data$ip == '+', 1, 0)
mod <- glm(ip.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### SME+/EMG+
tab <- with(data, table(sem, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ sem + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(sem, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
sem.ind <- ifelse(data$sem == '+', 1, 0)
mod <- glm(sem.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Anti-AChR
tab <- with(data, table(AChR, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ AChR + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(AChR, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
AChR.ind <- ifelse(data$AChR == '+', 1, 0)
mod <- glm(AChR.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Anti-AChR for de novo MG
subset <- which(data$denova == 'Denovo')
tab <- with(data, table(AChR[subset], mg.type[subset]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ AChR + Sex + Age, family = 'binomial', data = data, 
           subset = subset,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(AChR[subset], mg.type[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
AChR.ind <- ifelse(data$AChR == '+', 1, 0)
mod <- glm(AChR.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           subset = subset, control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)



#################### MG alone vs MG with myositis/myocarditis
#################### Anti-AChR Level
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
mg.type <- data$mg.type
subset <- which(is.na(data$mg.type) == FALSE)
AChR.level <- data$AChR.level
tab <- c(round(median(AChR.level[which(mg.type == 'MG alone')], na.rm = TRUE), 2),
         round(median(AChR.level[which(mg.type == 'MG with myositis/myocarditis')], na.rm = TRUE), 2),
         round(median(AChR.level[subset], na.rm = TRUE), 2),
         iqr(AChR.level[which(mg.type == 'MG alone')]),
         iqr(AChR.level[which(mg.type == 'MG with myositis/myocarditis')]),
         iqr(AChR.level[subset]))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('MG alone', 'MG with myositis/myocarditis', 'All')
print(tab)

# median regression
mod <- rq(AChR.level ~ mg.type + Age + Sex, tau = 0.50, data = data)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2) 
print(coef.ci)



#################### MG alone vs MG with myositis/myocarditis
#################### Anti-striational
tab <- with(data, table(Striated, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ Striated + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(Striated, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
Striated.ind <- ifelse(data$Striated == '+', 1, 0)
mod <- glm(Striated.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- -coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)




#################### MG alone vs MG with myositis/myocarditis
#################### Anti PD1/PDL1 monotherapy
tab <- with(data, table(PD1, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ PD1 + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(PD1, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
PD1.ind <- ifelse(data$PD1 == '+', 1, 0)
mod <- glm(PD1.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)



#################### MG alone vs MG with myositis/myocarditis
#################### Anti CTLA-4 monotherapy
tab <- with(data, table(CTLA4, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ CTLA4 + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(CTLA4, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
CTLA4.ind <- ifelse(data$CTLA4 == '+', 1, 0)
mod <- glm(CTLA4.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### MG alone vs MG with myositis/myocarditis
#################### Anti CTLA-4 and PD1/PDL1 combination therapy
tab <- with(data, table(CTLA4_PD1, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ CTLA4_PD1 + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


tab <- with(data, t(table(CTLA4_PD1, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
CTLA4_PD1.ind <- ifelse(data$CTLA4_PD1 == '+', 1, 0)
mod <- glm(CTLA4_PD1.ind ~ mg.type + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Distribution of treatment categories
tab <- with(data, t(table(treatment, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:3], ' (C)'),
                   paste0(colnames(tab)[1:3], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 4:6] <- tab[3, 1:3] / sum(tab[3, 1:3])
print(tab)

# likelihood ratio test
subset <- which(is.na(data$mg.type) == FALSE)
mod1 <- multinom(treatment ~ mg.type + Sex + Age, subset = subset, data = data)
mod2 <- multinom(treatment ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)


#################### MG alone vs MG with myositis/myocarditis
#################### Hepatitis
tab <- with(data, table(HP, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ HP + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(HP, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
HP.ind <- ifelse(data$HP == '+', 1, 0)
mod <- glm(HP.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)



#################### MG alone vs MG with myositis/myocarditis
#################### Corticosteriods
tab <- with(data, table(cort, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ cort + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(cort, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
cort.ind <- ifelse(data$cort == '+', 1, 0)
mod <- glm(cort.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Plasmapheresis
tab <- with(data, table(plasmapheresis, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ plasmapheresis + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(plasmapheresis, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
plasmapheresis.ind <- ifelse(data$plasmapheresis == '+', 1, 0)
mod <- glm(plasmapheresis.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Cardiac treatment
tab <- with(data, table(ct, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ ct + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(ct, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
ct.ind <- ifelse(data$ct == '+', 1, 0)
mod <- glm(ct.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Non-steroid immunomodulators 
tab <- with(data, table(imm, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ imm + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(imm, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
imm.ind <- ifelse(data$imm == '+', 1, 0)
mod <- glm(imm.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Outcome
tab <- with(data, t(table(outcome, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:5], ' (C)'),
                   paste0(colnames(tab)[1:5], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, (1:5) + 5] <- tab[3, 1:5] / sum(tab[3, 1:5])
print(tab)

# likelihood ratio test
outcome <- as.factor(data$outcome)
subset <- which(is.na(data$mg.type) == FALSE)
mod1 <- multinom(outcome ~ mg.type + Sex + Age, subset = subset, data = data)
mod2 <- multinom(outcome ~ Sex + Age, subset = subset, data = data)
anova(mod1, mod2)

# logistic regression
death.ind <- ifelse(outcome == 'Death', 1, 0)
mod <- glm(death.ind ~ mg.type + Sex + Age, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Age
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
Age <- data$Age
subset <- which(is.na(data$mg.type) == FALSE)
mg.type <- data$mg.type
tab <- c(round(median(Age[which(mg.type == 'MG alone')], na.rm = TRUE), 2),
         round(median(Age[which(mg.type == 'MG with myositis/myocarditis')], na.rm = TRUE), 2), 
         round(median(Age[subset], na.rm = TRUE), 2),
         iqr(Age[which(mg.type == 'MG alone')]),
         iqr(Age[which(mg.type == 'MG with myositis/myocarditis')]),
         iqr(Age[subset]))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('MG alone', 'MG with myositis/myocarditis', 'All')
print(tab)

# median regression
mod <- rq(Age ~ mg.type + Sex, tau = 0.50)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)

#################### MG alone vs MG with myositis/myocarditis
#################### Sex
tab <- with(data, table(Sex, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(Sex, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
sex.ind <- ifelse(data$Sex == 'F', 1, 0)
mod <- glm(sex.ind ~ mg.type + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### MG alone vs MG with myositis/myocarditis
#################### MSA
tab <- with(data, table(msa, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ msa + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(msa, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
msa.ind <- ifelse(data$msa == '+', 1, 0)
mod <- glm(msa.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)



#################### MG alone vs MG with myositis/myocarditis
#################### Anti-MuSK
tab <- with(data, table(musk, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ musk + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(musk, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
musk.ind <- ifelse(data$musk == '+', 1, 0)
mod <- glm(msa.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### Anti-LRP4
tab <- with(data, table(LRP4, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ LRP4 + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(LRP4, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
LRP4.ind <- ifelse(data$LRP4 == '+', 1, 0)
mod <- glm(LRP4.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### pre-existing and de novo MG
tab <- with(data, table(denova, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ denova + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(denova, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
denova.ind <- ifelse(data$denova == 'Denovo', 1, 0)
mod <- glm(denova.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### cardiac involvement
tab <- with(data, table(Cardiac, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ Cardiac + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(Cardiac, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
Cardiac.ind <- ifelse(data$Cardiac == '+', 1, 0)
mod <- glm(Cardiac.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### MG alone vs MG with myositis/myocarditis
#################### confirmed diagnosis of MG
tab <- with(data, table(confirm_diag, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ confirm_diag + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(confirm_diag, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
confirm_diag.ind <- ifelse(data$confirm_diag == '+', 1, 0)
mod <- glm(confirm_diag.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)



#################### MG alone vs MG with myositis/myocarditis
#################### interval
mg.type <- data$mg.type
subset <- which(is.na(mg.type) == FALSE)
interval <- data$interval
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
tab <- c(round(median(interval[which(mg.type == 'MG alone')], na.rm = TRUE), 2), 
         round(median(interval[which(mg.type == 'MG with myositis/myocarditis')], na.rm = TRUE), 2),
         round(median(interval[subset], na.rm = TRUE), 2),
         round(mean(interval[which(mg.type == 'MG alone')], na.rm = TRUE), 2), 
         round(mean(interval[which(mg.type == 'MG with myositis/myocarditis')], na.rm = TRUE), 2),
         round(mean(interval[subset], na.rm = TRUE), 2),
         round(sd(interval[which(mg.type == 'MG alone')], na.rm = TRUE), 2), 
         round(sd(interval[which(mg.type == 'MG with myositis/myocarditis')], na.rm = TRUE), 2),
         round(sd(interval[subset], na.rm = TRUE), 2),
         iqr(interval[which(mg.type == 'MG alone')]), 
         iqr(interval[which(mg.type == 'MG with myositis/myocarditis')]),
         iqr(interval[subset]))
tab <- matrix(tab, 4, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'Mean', 'STD', 'IQR')
colnames(tab) <- c('MG alone', 'MG with myositis/myocarditis', 'All')
print(tab)
print('Median regression')

# p-value
mod <- rq(interval ~ mg.type + Age + Sex, tau = 0.50, data = data)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary.rq(mod, se = 'ker')$coefficients[2, 1]
ses <- summary.rq(mod, se = 'ker')$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)


# mean regression
mod <- lm(interval ~ mg.type + Age + Sex, tau = 0.50, data = data)
# p-value
p <- coef(summary(mod, se = 'ker'))[2, 4]
print(p)

# coefficient estimate
coefs <- summary(mod)$coefficients[2, 1]
ses <- summary(mod)$coefficients[2, 2]
coef.ci <- cbind(coefs, coefs - qnorm(.975) * ses, coefs + qnorm(.975) * ses)
colnames(coef.ci) <- c('Coefficient estimate', 'Lower bound', 'Upper bound')
coef.ci <- round(coef.ci, 2)
print(coef.ci)


#################### MG alone vs MG with myositis/myocarditis
#################### fatigable or flunctuating
tab <- with(data, table(ff, mg.type))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mg.res.ind <- ifelse(data$mg.type == 'MG alone', 1, 0)
mod <- glm(mg.res.ind ~ ff + Sex + Age, family = 'binomial', data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# prepare table
tab <- with(data, t(table(ff, mg.type)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
ff.ind <- ifelse(data$ff == '+', 1, 0)
mod <- glm(ff.ind ~ mg.type + Sex + Age, family = 'binomial', data = data,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Combine myositis and MG
#################### Death and AChR
subset <- which(data$non.derm.IBM == '-')
death <- ifelse(data$outcome == 'Death', "+", '-')
death.ind <- ifelse(death == '+', 1, 0)
tab <- with(data, table(AChR[subset], death[subset]))
tab <- cbind(tab, tab / rowSums(tab)) 
rownames(tab) <- c('AChR (-)', 'AChR (+)')
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mod <- glm(death.ind ~ AChR + Sex + Age, data = data, 
           family = 'binomial',
           subset = subset, control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# logistic regression
# interaction effects between AChR and myocarditis
mod <- glm(death.ind ~ AChR * myocarditis + Sex + Age, data = data, family = 'binomial', 
           subset = subset,
           control = glm.control(maxit = 50))

# p-value (main effects of AChR)
print(coef(summary(mod))[2, 4])
# p-value (main effect of myocarditis)
print(coef(summary(mod))[3, 4])
# p-value (interactions effects between AChR and myocarditis)
print(coef(summary(mod))[6, 4])
# odds ratio (interactions effects between AChR and myocarditis)
print(exp(coef(summary(mod))[6, 1]))

# logistic regression
# interaction effects between AChR and myositis
mod <- glm(death.ind ~ AChR * my + Sex + Age, data = data, family = 'binomial', 
           subset = subset, control = glm.control(maxit = 50))

# p-value (main effect of AChR)
print(coef(summary(mod))[2, 4])
# p-value (main effect of myositis) 
print(coef(summary(mod))[3, 4])
# p-value (interactions effects between AChR and myositis) 
print(coef(summary(mod))[6, 4])

# logistic regression
# interaction effects between AChR and MG
mod <- glm(death.ind ~ AChR * MG + Sex + Age, data = data, family = 'binomial', 
           subset = subset, control = glm.control(maxit = 50))

# p-value (main effect of AChR)
print(coef(summary(mod))[2, 4])
# p-value (main effect of MG) 
print(coef(summary(mod))[3, 4])
# p-value (interactions effects between AChR and MG) 
print(coef(summary(mod))[6, 4])


print('Interaction effects between AChR and striational antibody on death')
mod <- glm(death.ind ~ AChR * Striated + Sex + Age, data = data, family = 'binomial', 
           subset = subset,
           control = glm.control(maxit = 50))

# p-value (main effect of AChR)
print(coef(summary(mod))[2, 4])
# p-value (main effect of striational antibody) 
print(coef(summary(mod))[3, 4])
# p-value (interactions effects between AChR and striational antibody) 
print(coef(summary(mod))[6, 4])



# prepare table
tab <- with(data, t(table(AChR[subset], death[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
rownames(tab) <- c('Alive', 'Dead')
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
AChR.ind <- ifelse(data$AChR == '+', 1, 0)
mod <- glm(AChR.ind ~ death + Sex + Age, data = data, family = 'binomial', 
           subset = subset,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Combine myositis and MG
#################### Outcome and AChR
subset <- which(data$non.derm.IBM == '-')
tab <- with(data, t(table(outcome[subset], AChR[subset])))
tab <- cbind(tab, tab / rowSums(tab))
rownames(tab) <- c('AChR (-)', 'AChR (+)')
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset2 <- with(data, which(is.na(AChR) == FALSE & non.derm.IBM == '-'))
mod1 <- multinom(outcome ~ Sex + Age, data = data, subset = subset2)
mod2 <- multinom(outcome ~ AChR + Sex + Age, subset = subset2, data = data)
anova(mod1, mod2)

# likelihood ratio test
# interaction effects
subset3 <- with(data, which(is.na(AChR) == FALSE & is.na(myocarditis) == FALSE & non.derm.IBM == '-'))
mod1 <- multinom(outcome ~ AChR + myocarditis + Sex + Age, data = data, 
                 subset = subset3)
mod2 <- multinom(outcome ~ AChR * myocarditis + Sex + Age, data = data, 
                 subset = subset3)
anova(mod1, mod2)


#################### Combine myositis and MG
#################### Death and Striational antibody
subset <- which(data$non.derm.IBM == '-')
death <- ifelse(data$outcome == 'Death', '+', '-')
death.ind <- ifelse(death == '+', 1, 0)
tab <- with(data, table(Striated[subset], death[subset]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
rownames(tab) <- c('striational antibody (-)', 'striational antibody (+)')

# logistic regression
mod <- glm(death.ind ~ Striated + Sex + Age, family = 'binomial', 
           subset = subset, data = data,
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# prepare table
tab <- with(data, t(table(Striated[subset], death[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
rownames(tab) <- c('Alive', 'Dead')
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
Striated.ind <- ifelse(data$Striated == '+', 1, 0)
mod <- glm(Striated.ind ~ death + Sex + Age, family = 'binomial', 
           subset = subset, data = data, 
           control = glm.control(maxit = 50))


# p-value
print(coef(summary(mod))[2, 4])
# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Combine myositis and MG
#################### Death and Cardiovascular treatment
subset <- which(data$non.derm.IBM == '-')
death <- ifelse(data$outcome == 'Death', '+', '-')
death.ind <- ifelse(death == '+', 1, 0)
tab <- with(data, table(ct[subset], death[subset]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
rownames(tab) <- c('Cardiovascular treatment (-)', 'Cardiovascular treatment (+)')
print(tab)

# logistic regression
mod <- glm(death.ind ~ ct + Sex + Age, family = 'binomial', 
           subset = subset, data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)



#################### Combine myositis and MG
#################### Outcome and Cardiovascular treatment
subset <- which(data$non.derm.IBM == '-')
tab <- with(data, t(table(outcome[subset], ct[subset])))
tab <- cbind(tab, tab / rowSums(tab))
rownames(tab) <- c('Cardiovascular treatment (-)', 'Cardiovascular treatment (+)')
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset2 <- with(data, which(non.derm.IBM == '-' & is.na(ct) == FALSE))
mod1 <- multinom(outcome ~ Sex + Age, subset = subset2, data = data)
mod2 <- multinom(outcome ~ ct + Sex + Age, data = data,
                 subset = subset2)
anova(mod1, mod2)

#################### Combine myositis and MG
#################### Death and respiratory support
subset <- which(data$non.derm.IBM == '-')
death <- ifelse(data$outcome == 'Death', '+', '-')
death.ind <- ifelse(death == '+', 1, 0)
tab <- with(data, table(rs[subset], death[subset]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)',
                                                 ' (P)', ' (P)'))
rownames(tab) <- c('respiratory support (-)', 'respiratory support (+)')
print(tab)

# logistic regression
mod <- glm(death.ind ~ rs + Sex + Age, family = 'binomial', 
           subset = subset, data = data, 
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Combine myositis and MG
#################### Outcome and respiratory support
subset <- which(data$non.derm.IBM == '-')
tab <- with(data, t(table(outcome[subset], rs[subset])))
tab <- cbind(tab, tab / rowSums(tab))
rownames(tab) <- c('respiratory support (-)', 'respiratory support (+)')
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)', ' (P)', ' (P)'))
print(tab)

# likelihood ratio test
subset2 <- with(data, which(non.derm.IBM == '-' & is.na(rs) == FALSE))
mod1 <- multinom(outcome ~ Sex + Age, subset = subset2, data = data)
mod2 <- multinom(outcome ~ rs + Sex + Age, data = data,
                 subset = subset2)
anova(mod1, mod2)

#################### Combine myositis and MG
#################### Death and other predictors
subset <- which(data$non.derm.IBM == '-')
death <- ifelse(data$outcome == 'Death', '+', '-')
death.ind <- ifelse(death == '+', 1, 0)
mod <- glm(death.ind ~ AChR + ct + Striated + rs + Sex + Age, family = 'binomial', 
           subset = subset, data = data,
           control = glm.control(maxit = 50))

# Odds ratio
OR <- exp(coef(mod))

# 95% confidence intervals for odds ratio
CI <- exp(confint(mod))

# p-values
pvals <- summary(mod)$coefficients[, 4]

# Combine into one table
results <- cbind(OR, CI, pvals)
results <- results[2:5, ]
colnames(results) <- c("Odds ratio", "2.5 %", "97.5 %", "p-value")
rownames(results) <- c('AChR', 'Cardiovascular treatment', 'StrAbs', 'Respiratory support')
round(results, 3)









#################### Comparison with Table 1 in Fenioux et al (2023)
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')

# prepare table
tab <- table(tet.my)
tab <- c(tab, tab / sum(tab))
names(tab) <- mapply(paste0, names(tab), c(' (C)', ' (C)', ' (P)', " (P)"))
tab <- rbind(tab, c(767, 28, 767 / (767 + 28), 28 / (767 + 28)))
rownames(tab) <- c('Table 2', 'Fenioux et. al (2023)')


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Female
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')

# prepare table
tab <- with(data, table(tet.my, Sex[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(252, 764 - 252, 252 / 764, (764 - 252) / (764),
                           14, 13, 14 / 27, 13 / 27), byrow = TRUE,
                         2, 4))  

require("DescTools")
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Sex = c("F", "M"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)

tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab, t(apply(tab[, 1:2], 1, function(x) x / sum(x))))
colnames(tab.pre) <- c('Cancer subtypes', 'Study', 'F (C)', 'M (C)', 'F (P)', 'M (P)')
rownames(tab.pre) <- NULL
print(tab.pre)
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Age
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
Age <- data$Age
myocarditis <- data$myocarditis
tet <- ifelse(data$cancer.detailed == 'thymic', '+', '-')
tab <- c(round(median(Age[which(myocarditis == '+' & tet == '-')], na.rm = TRUE), 2),
         round(median(Age[which(myocarditis == '+' & tet == '+')], na.rm = TRUE), 2),
         round(median(Age[which(myocarditis == '+')], na.rm = TRUE), 2),
         iqr(Age[which(myocarditis == '+' & tet == '-')]),
         iqr(Age[which(myocarditis == '+' & tet == '+')]),
         iqr(Age[which(myocarditis == '+')]))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Other cancers', 'TET', 'All')
print(tab)


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Prior history of MG
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')

tab <- with(data, table(tet.my, MG.PD[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(723 - 3, 3, (723 - 3) / (723), 3 / 723, 
                           23, 2, 23 / 25, 2 / 25), byrow = TRUE,
                         2, 4))  


tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Preexisting.MG = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)

tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
tab.pre

BreslowDayTest(x = tabs, OR = 1, correct = TRUE)


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Anti-PD-1/PD-L1
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, PD1[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(759 - 587, 587, (759 - 587) / (759), 587 / 759, 
                           1, 26, 1 / 27, 26 / 27), byrow = TRUE,
                         2, 4))  

tabs <- xtabs(freq ~ .,
              cbind(expand.grid(PD1 = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)

tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# comparison
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Anti-CTLA-4
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, CTLA4[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(759 - 13, 13, (759 - 13) / (759), 13 / 759, 
                           27, 0, 27 / 27, 0 / 27), byrow = TRUE,
                         2, 4))  

tabs <- xtabs(freq ~ .,
              cbind(expand.grid(CTLA4 = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)


tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)



#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Anti CTLA-4 and PD1/PDL1 combination therapy
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, CTLA4_PD1[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(759 - 155, 155, (759 - 155) / (759), 155 / 759, 
                           26, 1, 26 / 27, 1 / 27), byrow = TRUE,
                         2, 4))  

tabs <- xtabs(freq ~ .,
              cbind(expand.grid(CTLA4_PD1 = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)


tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)


# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### MG-like syndrome
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, MGS[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(767 - 157, 157, (767 - 157) / (767), 157 / 767, 
                           10, 18, 10 / 28, 18 / 28), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(MG.syndrome = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)



#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Myositis/rhabdomyolysis
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, MR[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(767 - 280, 280, (767 - 280) / (767), 280 / 767, 
                           7, 21, 7 / 28, 21 / 28), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Myositis_rhabdomyolysis = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Definite myocarditis
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, dm[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(696 - 221, 221, (696 - 221) / (696), 221 / 696, 
                           22 - 9, 9, (22 - 9) / 22, (9) / 22), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Definite_myocarditis = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Probable myocarditis
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, pm[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(696 - 262, 262, (696 - 262) / (696), 262 / 696, 
                           22 - 4, 4, (22 - 4) / 22, (4) / 22), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Probable_myocarditis = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)



#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Possible myocarditis
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, pm2[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(696 - 204, 204, (696 - 204) / (696), 204 / 696, 
                           22 - 7, 7, (22 - 7) / 22, (7) / 22), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Possible_myocarditis = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)



#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Peak troponin
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
pt_norm <- data$pt_norm
myocarditis <- data$myocarditis
tet <- ifelse(data$cancer.detailed == 'thymic', '+', '-')
tab <- c(round(median(pt_norm[which(myocarditis == '+' & tet == '-')], na.rm = TRUE), 2),
         round(median(pt_norm[which(myocarditis == '+' & tet == '+')], na.rm = TRUE), 2),
         round(median(pt_norm[which(myocarditis == '+')], na.rm = TRUE), 2),
         iqr(pt_norm[which(myocarditis == '+' & tet == '-')]),
         iqr(pt_norm[which(myocarditis == '+' & tet == '+')]),
         iqr(pt_norm[which(myocarditis == '+')]))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Other cancers', 'TET', 'All')
print(tab)


#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Peak CK
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
pck_norm <- data$pck_norm
myocarditis <- data$myocarditis
tet <- ifelse(data$cancer.detailed == 'thymic', '+', '-')
tab <- c(round(median(pck_norm[which(myocarditis == '+' & tet == '-')], na.rm = TRUE), 2),
         round(median(pck_norm[which(myocarditis == '+' & tet == '+')], na.rm = TRUE), 2),
         round(median(pck_norm[which(myocarditis == '+')], na.rm = TRUE), 2),
         iqr(pck_norm[which(myocarditis == '+' & tet == '-')]),
         iqr(pck_norm[which(myocarditis == '+' & tet == '+')]),
         iqr(pck_norm[which(myocarditis == '+')]))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('Other cancers', 'TET', 'All')
print(tab)



#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Corticosteriods
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, cort[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(767 - 625, 625, (767 - 625) / (767), 625 / 767, 
                           28 - 27, 27, (28 - 27) / 28, 27 / 28), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Corticosteriods = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)




#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Abatacept
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, abatacept[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(767 - 82, 82, (767 - 82) / (767), 82 / 767, 
                           28 - 5, 5, (28 - 5) / 28, 5 / 28), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Abatacept = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)



#################### Comparison with Table 1 in Fenioux et al (2023)
#################### Intravenous immune globulin
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, iig[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(767 - 71, 71, (767 - 71) / (767), 71 / 767, 
                           28 - 13, 13, (28 - 13) / 28, 13 / 28), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(Intravenous_immune_globulin = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)



#################### Comparison with Table 1 in Fenioux et al (2023)
#################### plasmapheresis
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
tab <- with(data, table(tet.my, plasmapheresis[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(767 - 82, 82, (767 - 82) / (767), 82 / 767, 
                           28 - 7, 7, (28 - 7) / 28, 7 / 28), byrow = TRUE,
                         2, 4))  
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(plasmapheresis = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)



#################### Comparison with Table 1 in Fenioux et al (2023)
#################### death
cancer.type <- as.character(data$cancer.detailed)
cancer.type[which(is.na(cancer.type) == TRUE)] <- 'N/A'
tet <- ifelse(cancer.type == 'thymic', '+', '-')
tet.my <- tet[which(data$myocarditis == '+')]
tet.my <- ifelse(tet.my == '+', 'TET', 'Other cancers')
death <- ifelse(data$outcome == 'Death', '+', '-')
tab <- with(data, table(tet.my, death2[which(myocarditis == '+')]))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, matrix(c(738 - 122, 122, (738 - 122) / (738), 122 / 738, 
                           28 - 10, 10, (28 - 10) / 28, 10 / 28), byrow = TRUE,
                         2, 4)) 
tabs <- xtabs(freq ~ .,
              cbind(expand.grid(death2 = c("-", "+"),
                                Cancer = c("Other cancers", "TET"),
                                study = c("Table 2", "Fenioux et. al (2023)")),
                    freq = as.numeric(t(tab[, 1:2])))
)
tab.pre <- data.frame(c('Other cancers', 'TET', 'Other cancers', 'TET'),
                      c('Table 2', 'Table 2',
                        "Fenioux et. al (2023)", "Fenioux et. al (2023)"), 
                      tab)
colnames(tab.pre) <- c('Cancer subtypes', 'Study', '- (C)', '+ (C)', '- (P)', '+ (P)')
rownames(tab.pre) <- NULL
print(tab.pre)

# Breslow-Day Test
BreslowDayTest(x = tabs, OR = 1, correct = TRUE)

