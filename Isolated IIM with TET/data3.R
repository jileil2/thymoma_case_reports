##################### set working directory
setwd("/Users/jileilin/Desktop/Research/thymomal case reports/")

##################### load packages
library('quantreg')
library('nnet')

data <- age
data <- as.data.frame(data)
colnames(data) <- 'Age'
data$Sex <- Sex
data$Diagnosis.Thymoma <- dt
write.csv(data, 'data3.csv')

##################### read data
data1 <- read.csv('data1.csv')

##################### Age within thymoma and thymic carcinoma
age <- data$Age
dt <- data$Diagnosis.Thymoma
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))

tab <- c(round(median(age[which(dt == 'thymoma')], na.rm = TRUE), 2),
         round(median(age[which(dt == 'thymic carcinoma')], na.rm = TRUE), 2), 
         round(median(age, na.rm = TRUE), 2),
         iqr(age[which(dt == 'thymoma')]),
         iqr(age[which(dt == 'thymic carcinoma')]),
         iqr(age))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('thymoma', 'thymic carcinoma', 'All')
knitr::kable(tab, caption = "Age within thymoma and thymic carcinoma")

# median regression
mod <- rq(age ~ dt + Sex, tau = 0.50, subset = which(dt != 'no thymic tumor'))
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

##################### Sex within thymoma and thymic carcinoma
Sex <- data$Sex
dt <- data$Diagnosis.Thymoma
tab <- t(table(Sex, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'), paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
sex.ind <- ifelse(Sex == 'F', 1, 0)
mod <- glm(sex.ind ~ dt + age, family = 'binomial', subset = which(dt != 'no thymic tumor'),
           control = glm.control(maxit = 50))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

##################### Myositis distribution
my <- data$myositis
dt <- data$Diagnosis.Thymoma
age <- data$Age
subset <- which(dt != 'no thymic tumor')
tab <- t(table(my, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:3], ' (C)'),
                   paste0(colnames(tab)[1:3], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])


# likelihood ratio test
mod1 <- multinom(my ~ dt + Sex + age, subset = subset)
mod2 <- multinom(my ~ Sex + age, subset = subset)
anova(mod1, mod2)

##################### Cardiac involvement
Cardiac <- data$Cardiac
dt <- data$Diagnosis.Thymoma
age <- data$Age
Sex <- data$Sex
subset <- which(dt != 'no thymic tumor')

tab <- t(table(Cardiac, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])

# logistic regression
cardiac.ind <- ifelse(Cardiac == 'Y', 1, 0)
mod <- glm(cardiac.ind ~ dt + age + Sex, family = 'binomial', 
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

##################### O+
Op <- data$Op
dt <- data$Diagnosis.Thymoma
subset <- which(dt != 'no thymic tumor')
tab <- t(table(Op, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])

# logistic regression
Op.ind <- ifelse(Op == '+', 1, 0)
mod <- glm(Op.ind ~ dt + age + Sex, family = 'binomial', subset = subset, 
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


##################### Striational antibodies
tab <- with(data, t(table(Striated, Diagnosis.Thymoma)))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
print(tab)


##################### MSA
tab <- with(data, t(table(MSA, Diagnosis.Thymoma)))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

##################### CK+
subset <- which(data$Diagnosis.Thymoma != 'no thymic tumor')
tab <- with(data, t(table(CK, Diagnosis.Thymoma)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])

# logistic regression
CKp.ind <- ifelse(data$CK == '+', 1, 0)
mod <- glm(CKp.ind ~ Diagnosis.Thymoma + Age + Sex, family = 'binomial', 
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


##################### CK level
pck <- data$CK.level
dt <- data$Diagnosis.Thymoma
subset <- which(dt != 'no thymic tumor')
age <- data$Age
Sex <- data$Sex
tab <- c(round(median(pck[which(dt == 'thymoma')], na.rm = TRUE), 2),
         round(median(pck[which(dt == 'thymic carcinoma')], na.rm = TRUE), 2), 
         round(median(pck, na.rm = TRUE), 2),
         iqr(pck[which(dt == 'thymoma')]),
         iqr(pck[which(dt == 'thymic carcinoma')]),
         iqr(pck))
tab <- matrix(tab, 2, 3, byrow = TRUE)
rownames(tab) <- c('Median', 'IQR')
colnames(tab) <- c('thymoma', 'thymic carcinoma', 'All')

# median regression
mod <- rq(pck ~ dt + age + Sex, tau = 0.50, subset = subset)
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


##################### Corticosteriods
cort <- data$cort
dt <- data$Diagnosis.Thymoma
subset <- which(dt != 'no thymic tumor')
age <- data$Age
Sex <- data$Sex
tab <- t(table(cort, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
cort.ind <- ifelse(cort == '+', 1, 0)
mod <- glm(cort.ind ~ dt + age + Sex, family = 'binomial', 
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


##################### plasmapheresis
plasmapheresis <- data$plasmapheresis
dt <- data$Diagnosis.Thymoma
subset <- which(dt != 'no thymic tumor')
age <- data$Age
Sex <- data$Sex
tab <- t(table(plasmapheresis, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
plasmapheresis.ind <- ifelse(plasmapheresis == '+', 1, 0)
mod <- glm(plasmapheresis.ind ~ dt + age + Sex, family = 'binomial', 
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


##################### Respiratory support
rs <- data$rs
dt <- data$Diagnosis.Thymoma
subset <- which(dt != 'no thymic tumor')
age <- data$Age
Sex <- data$Sex
tab <- t(table(rs, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
rs.ind <- ifelse(rs == '+', 1, 0)
mod <- glm(rs.ind ~ dt + age + Sex, family = 'binomial', 
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


##################### Cardiovascular treatments
ct <- data$ct
dt <- data$Diagnosis.Thymoma
subset <- which(dt != 'no thymic tumor')
age <- data$Age
Sex <- data$Sex
tab <- t(table(ct, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
ct.ind <- ifelse(ct == '+', 1, 0)
mod <- glm(ct.ind ~ dt + age + Sex, family = 'binomial', 
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


##################### Non-steroid immunomodulators
imm <- data$imm
dt <- data$Diagnosis.Thymoma
subset <- which(dt != 'no thymic tumor')
age <- data$Age
Sex <- data$Sex
tab <- t(table(ct, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[4] <- 'All'
tab[4, 3:4] <- tab[4, 1:2] / sum(tab[4, 1:2])
print(tab)

# logistic regression
imm.ind <- ifelse(imm == '+', 1, 0)
mod <- glm(imm.ind ~ dt + age + Sex, family = 'binomial', 
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


##################### outcome
outcome <- data$outcome
dt <- data$Diagnosis.Thymoma
subset <- which(dt != 'no thymic tumor')
age <- data$Age
Sex <- data$Sex
tab <- t(table(outcome, dt))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', " (C)",
                                                 ' (P)', ' (P)', ' (P)', " (P)"))

# likelihood ratio test
mod1 <- multinom(outcome ~ dt + Sex + age, subset = subset)
mod2 <- multinom(outcome ~ Sex + age, subset = subset)
anova(mod1, mod2)

# logistic regression
death.ind <- ifelse(outcome == 'Death', 1, 0)
mod <- glm(death.ind ~ dt + age + Sex, family = 'binomial', 
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

##################### compare table 1 and table 3 in terms of thymoma and thymic carcinoma



##################################################### sex (table 3)
Sex3 <- dat3$Sex

##################################################### modeling
dt <- data$Diagnosis.Thymoma
dt[which(dt == 'no thymic tumor')] <- NA
dt.table1 <- data1$thymoma.type
dat.ana <- data.frame(age = c(data1$Age.TP, data$Age),
                      sex = c(data1$Sex, data$Sex),
                      TP.type = ifelse(c(dt.table1, dt) == 'thymoma', 1, 0),
                      group = c(rep('table 1', nrow(data1)),
                                rep('table 3', nrow(data))))

tab <- t(table(dat.ana$TP.type, dat.ana$group))
colnames(tab) <- c('thymic carcinoma', 'thymoma')
tab <- cbind(tab, tab / rowSums(tab))
rownames(tab) <- c('table 1', 'table 2') 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', 
                                                 ' (P)', ' (P)'))
print(tab)

# logistic regression
mod <- glm(TP.type ~ group + age + sex, data = dat.ana, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- -coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)