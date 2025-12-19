##################### set working directory
setwd("/Users/jileilin/Desktop/Research/thymomal case reports/")

##################### load packages
library("strucchange")
library('nnet')
library('quantreg')
library('multcomp')

#################### read data
data <- read.csv('data1.csv')[, -1]


#################### Gender distribution

## logistic regression
sex.ind <- ifelse(data$Sex == 'F', 1, 0)
mod <- glm(sex.ind ~ EOMG, family = 'binomial')

## p-value
p <- coef(summary(mod))[2, 4]
print(p)

## odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Age distribution
out <- c(summary(na.omit(data$Age)), sd(na.omit(data$Age)))
names(out)[length(out)] <- 'sd'
print(out)

# male
out <- c(summary(na.omit(data$Age[data$Sex == 'M'])), sd(na.omit(data$Age[data$Sex == 'M'])))
names(out)[length(out)] <- 'sd'
print(out)

# female
out <- c(summary(na.omit(data$Age[data$Sex == 'F'])), sd(na.omit(data$Age[data$Sex == 'F'])))
names(out)[length(out)] <- 'sd'
print(out)

#################### Age of myositis onset
out <- c(summary(na.omit(data$Age.MY)), sd(na.omit(data$Age.MY)))
names(out)[length(out)] <- 'sd'
print(out)

# male
out <- c(summary(na.omit(data$Age.MY[data$Sex == 'M'])), sd(na.omit(data$Age.MY[data$Sex == 'M'])))
names(out)[length(out)] <- 'sd'
print(out)

# female
out <- c(summary(na.omit(data$Age.MY[data$Sex == 'F'])), sd(na.omit(data$Age.MY[data$Sex == 'F'])))
names(out)[length(out)] <- 'sd'
print(out)


#################### Proportion of EOMG and LOMG
# EOMG: earlier than 45
Count <- table(data$EOMG.45)
Prop <- Count / length(na.omit(data$EOMG.45))
Count <- as.character(Count)
Summary <- rbind(round(Prop, 4), Count)
rownames(Summary) <- c("Proportion", "Count")
colnames(Summary) <- c("EOMG", "LOMG")
print(Summary)

# EOMG: earlier than 50
Count <- table(data$EOMG.50)
Prop <- Count / length(na.omit(data$EOMG.50))
Count <- as.character(Count)
Summary <- rbind(round(Prop, 4), Count)
rownames(Summary) <- c("Proportion", "Count")
colnames(Summary) <- c("EOMG", "LOMG")
print(Summary)


#################### Prevalence of TET
#### by EOMG
tt <- table(data$TET, data$EOMG)
Count <- cbind(table(data$TET), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
TET.ind <- ifelse(data$TET == '1', 1, 0)
mod <- glm(TET.ind ~ data$EOMG + data$Sex, family = 'binomial')

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

#### by gender
tt <- table(data$TET, data$Sex)
Count <- cbind(table(data$TET), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)

colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
TET.ind <- ifelse(data$TET == '1', 1, 0)
mod <- glm(TET.ind ~ data$Sex + data$Age.TP, family = 'binomial')

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



#################### Prevalence of thymoma within TET
### between EOMG and LOMG
index.thymoma <- which(data$TET == 1)
tt <- table(data$thymoma[index.thymoma], EOMG[index.thymoma])
Count <- cbind(table(data$thymoma[index.thymoma]), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
TP.ind <- ifelse(data$thymoma == '+', 1, 0)
mod <- glm(TP.ind ~ data$EOMG + data$Sex, subset = index.thymoma, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


### between gender
tt <- table(data$thymoma[index.thymoma], data$Sex[index.thymoma])
Count <- cbind(table(TP.fp[index.thymoma]), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
Sex <- as.factor(data$Sex)
Sex <- relevel(Sex, ref = 'M')
mod <- glm(TP.ind ~ Sex + data$Age.TP, subset = index.thymoma, family = 'binomial')

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


#################### Subtypes of thymoma within TET
# thymoma subtype
Table <- table(data$thymoma.subtype)
print(Table)

# thymoma subtype by EOMG
tt <- table(data$thymoma.subtype, data$EOMG)
Count <- cbind(table(subtype.thymoma), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- round(cbind(t(Count), Prop),3)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# compare A/AB and B between EOMG and LOMG
subtype.TP.ind <- ifelse(data$thymoma.subtype == 'A/AB', 1, 0)
mod <- glm(subtype.TP.ind ~ data$EOMG + data$Sex, family = 'binomial')

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

#################### change-point detection (progression of mortality rate over time)
death <- 1 * (data$outcome == "Death")
year <- dat$Year
idx <- order(year)
year <- year[idx]
death <- death[idx]
unique.year <- as.vector(unique(year))
x <- death
X_ts <- ts(x, start = 1, end = length(x), frequency = 1)
bp <- breakpoints(X_ts ~ 1)
change_points <- breakpoints(bp)
breakpoint <- change_points$breakpoints
print(year[breakpoint])

#################### Compare outcome/mortality between thymoma and no thymoma
index <- which(data$year >= 1993 & is.na(data$thymoma) == FALSE)
tt <- table(data$outcome[index], data$thymoma[index])
Count <- cbind(table(data$outcome[index]), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table) <- c("All", 'Absent', 'Present')
print(Table)

# likelihood ratio test
mod1 <- nnet::multinom(outcome ~ Sex + Age.TP, subset = index, data = data)
mod2 <- nnet::multinom(outcome ~ thymoma + Sex + Age.TP, subset = index, , data = data)
anova(mod1, mod2)

### compare mortality between thymoma presence
death.ind <- ifelse(data$outcome == 'Death', 1, 0)
mod <- glm(death.ind ~ thymoma + Age.TP + Sex, subset = index, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

#################### The distribution of the thymoma subtypes in the 5 myositis subtypes
tt <- table(data$thymoma.subtype, data$myositis.subtype)
Count <- cbind(table(data$thymoma.subtype), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

### likelihood ratio test (adjust for multiple comparison)


# compare thymoma subtype between polymyositis and dermatomyositis
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'Dermatomyositis'))
mod1 <- nnet::multinom(thymoma.subtype ~ Sex + Age.TP, subset = subset, data = data)
mod2 <- nnet::multinom(thymoma.subtype ~ myositis.subtype + Sex + Age.TP, subset = subset, data = data)
anova(mod1, mod2)

# compare thymoma subtype between GCM/GrM and dermatomyositis
subset <- which(data$myositis.subtype %in% c('GCM/GrM', 'Dermatomyositis'))
mod1 <- nnet::multinom(thymoma.subtype ~ Sex + Age.TP, subset = subset, data = data)
mod2 <- nnet::multinom(thymoma.subtype ~ myositis.subtype + Sex + Age.TP, subset = subset, data = data)
anova(mod1, mod2)

# compare thymoma subtype between polymyositis and GCM/GrM
subset <- which(data$myositis.subtype %in% c('GCM/GrM', 'Polymyositis'))
mod1 <- nnet::multinom(thymoma.subtype ~ Sex + Age.TP, subset = subset, data = data)
mod2 <- nnet::multinom(thymoma.subtype ~ myositis.subtype + Sex + Age.TP, subset = subset, data = data)
anova(mod1, mod2)


#################### MG subtypes
## Separate by EOMG and LOMG
tt <- table(data$MG.subtype, data$EOMG)
Count <- cbind(table(data$MG.subtype), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- round(cbind(t(Count), Prop),3)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

## likelihood ratio test
mod1 <- nnet::multinom(MG.subtype ~ Sex, data = data)
mod2 <- nnet::multinom(MG.subtype ~ Sex + EOMG, data = data)
anova(mod1, mod2)


#################### Compare the temporal relationship between EOMG and LOMG
subset <- which(is.na(data$EOMG) == FALSE)
mod1 <- nnet::multinom(Earlier.onset ~ Sex, subset = subset, data = data)
mod2 <- nnet::multinom(Earlier.onset ~ Sex + EOMG, subset = subset, data = data)
anova(mod1, mod2)



#################### Lag times between primary and subsequent disorders
### run those two lines together
mod <- glm(lag ~ EOMG + Sex, family = poisson(link = "log"), data = data)
summary(mod)



#################### Prevalence of AChR between E/LOMG
tt <- table(data$AChR, data$EOMG)
Count <- cbind(table(Y), tt)
Prop <- t(Count) / apply(Count,2,sum)
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)", "(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
AChR.ind <- ifelse(data$AChR == 'Pos', 1, 0)
mod <- glm(AChR.ind ~ EOMG + Sex, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### AChR level by EOMG and LOMG
Level <- data$AChR.Level
idx.po <- which(data$AChR == 'Pos')
X <- as.numeric(Level[idx.po])
X <- X[-which(is.na(X))]
Table <- c(sum(!is.na(X)), summary(X), sd(X, na.rm = T))
X.EOMG <- as.numeric(Level[which(data$AChR =="Pos" & EOMG == "EOMG")])
X.LOMG <- as.numeric(Level[which(data$AChR=="Pos" & EOMG == "LOMG")])
Table = rbind(Table,
              c(sum(!is.na(X.EOMG)), summary(X.EOMG[-which(is.na(X.EOMG))]), sd(X.EOMG, na.rm = T)),
              c(sum(!is.na(X.LOMG)), summary(X.LOMG[-which(is.na(X.LOMG))]), sd(X.LOMG, na.rm = T)))
rownames(Table) <- c("All", "EOMG", "LOMG")
colnames(Table)[1] <- "Count"
print(Table)


#################### median of Titers between EOMG and LOMG
indx <- which(data$AChR == 'Pos' & is.na(data$EOMG) == FALSE)

### median regression adjusted for sex
mod <- rq(AChR.Level ~ EOMG + Sex, tau = 0.50, data = data, 
          subset = indx)
# p-value
coef(summary(mod, se = 'ker'))[2, 4]

# Median difference between EOMG and LOMG (regression coefficient, adjusted for gender) 
# with 95% confidence intervals
coef <- coef(summary(mod, se = 'ker'))[2, 1]
se <- coef(summary(mod, se = 'ker'))[2, 2]
pp <- c(coef, coef - se * qnorm(.975), coef + se * qnorm(.975))
names(pp) <- c('Estimate', 'Lower bound', 'Upper bound')
print(round(pp, 2))

### median regression 
mod <- rq(AChR.Level ~ EOMG, tau = 0.50, data = data)
# p-value
print(coef(summary(mod, se = 'ker'))[2, 4])

# Median difference between EOMG and LOMG
# with 95% confidence intervals
coef <- coef(summary(mod, se = 'ker'))[2, 1]
se <- coef(summary(mod, se = 'ker'))[2, 2]
pp <- c(coef, coef - se * qnorm(.975), coef + se * qnorm(.975))
names(pp) <- c('Estimate', 'Lower bound', 'Upper bound')
print(round(pp, 2))

#################### Striational Antibody with thymoma
tt <- table(data$Striated, data$thymoma)
Count <- cbind(table(Y2), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
st.ind <- ifelse(data$Striated == '+', 1, 0)
mod <- glm(st.ind ~ thymoma + Sex + Age.TP, family = 'binomial', data = data)

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


#################### MSA and EOMG
tt <- table(data$MSA, data$EOMG)
Count <- cbind(table(data$MSA), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
MSA.ind <- ifelse(data$MSA == '+', 1, 0)
mod <- glm(MSA.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### MAA prevalence between EOMG and LOMG
tt <- table(data$MAA, data$EOMG)
Count <- cbind(table(data$MAA), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), Prop)
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
MAA.ind <- ifelse(data$MAA == '+', 1, 0)
mod <- glm(MAA.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

#################### Cardiac involvement by EOMG/LOMG
Count <- table(data$Cardiac)
Prop <- Count / length(data$Cardiac)
Count <- as.character(Count)
Summary <- rbind(round(Prop, 4), Count)
rownames(Summary) <- c("Proportion", "Count")
colnames(Summary) <- c("No Cardiac", "Yes")
print(Summary)

tt <- table(data$Cardiac, data$EOMG)
Count <- cbind(table(data$Cardiac), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
Cardiac.ind <- ifelse(data$Cardiac == 'Yes', 1, 0)
mod <- glm(Cardiac.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


#################### Cardiac involvement and thymoma
tt <- table(data$Cardiac, data$thymoma)
Count <- cbind(table(data$Cardiac), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
Cardiac.ind <- ifelse(data$Cardiac == 'Yes', 1, 0)
subset <- which(data$diagnose.thymoma == 1)
mod <- glm(Cardiac.ind ~ thymoma + Age.TP + Sex, data = data, subset = subset, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


#################### Transaminitis between EOMG and LOMG
tt <- table(data$transaminitis, data$EOMG)
Count <- cbind(table(transaminitis), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
transaminitis.ind <- ifelse(data$transaminitis == 'elevated', 1, 0)
mod <- glm(transaminitis.ind ~ EOMG + Age + Sex, family = 'binomial', data = data,
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


#################### compare outcome distribution between EOMG and LOMG
idx <- which(data$year >= 1993)
tt <- table(data$outcome[idx], data$EOMG[idx])
Count <- cbind(table(data$outcome[idx]), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# likelihood ratio test
mod1 <- nnet::multinom(outcome ~ Sex, data = data, subset = intersect(idx, which(is.na(EOMG) == FALSE)))
mod2 <- nnet::multinom(outcome ~ EOMG + Sex, data = data, subset = idx)
anova(mod1, mod2)

# logistic regression (mortality rate)
death.ind <- ifelse(data$outcome == 'Death', 1, 0)
mod <- glm(death.ind ~ EOMG + Sex, subset = idx, data = data, family = 'binomial')

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


#################### Corticosteriods between EOMG and LOMG
index <- which(data$year >= 1993)
tab <- t(table(data$cort[index], data$EOMG[index]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)


# logistic regression
cort.ind <- ifelse(data$cort == '+', 1, 0)
mod <- glm(cort.ind ~ EOMG + Sex, family = 'binomial', subset = index, data = data, control = glm.control(maxit = 50))

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


#################### Corticosteriods between EOMG and LOMG
index <- which(data$year >= 1993)
tab <- t(table(data$plasmapheresis[index], data$EOMG[index]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
plasmapheresis.ind <- ifelse(data$plasmapheresis == '+', 1, 0)
mod <- glm(plasmapheresis.ind ~ EOMG + Sex, family = 'binomial', subset = index, data = data, control = glm.control(maxit = 50))

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


#################### Resipiratory support between EOMG and LOMG
index <- which(data$year >= 1993)
tab <- t(table(data$rs[index], data$EOMG[index]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
rs.ind <- ifelse(data$rs == '+', 1, 0)
mod <- glm(rs.ind ~ EOMG + Sex, family = 'binomial', subset = index, data = data,
           control = glm.control(maxit = 50))

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


#################### Cardiovascular treatments between EOMG and LOMG
index <- which(data$year >= 1993)
tab <- t(table(data$ct[index], data$Cardiac[index]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])

# logistic regression
ct.ind <- ifelse(data$ct == '+', 1, 0)
mod <- glm(ct.ind ~ EOMG + Sex, family = 'binomial', data = data,
           subset = index, control = glm.control(maxit = 50))

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


#################### Non-steroid immunomodulators between EOMG and LOMG
index <- which(data$year >= 1993)
tab <- t(table(data$imm[index], data$EOMG[index]))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

# logistic regression
imm.ind <- ifelse(data$imm == '+', 1, 0)
mod <- glm(imm.ind ~ EOMG + Sex, family = 'binomial', data = data,
           subset = index, control = glm.control(maxit = 50))

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

#################### GG positive rate between EOMG and LOMG
subset <- which(data$GG != 'nobiopsy')
tt <- table(data$GG[subset], data$EOMG[subset])
Count <- cbind(table(GG3), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
GG.ind <- ifelse(data$GG == "+", 1, 0)
mod <- glm(GG.ind ~ EOMG + Sex, subset = subset, data = data, family = 'binomial')

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



#################### mortality rate by thymoma after breakpoint
index <- which(data$year >= 1993)
tt <- table(data$outcome[index], data$thymoma[index])
Count <- cbind(rowSums(tt), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)", "(P)"), each = nrow(tt)))
rownames(Table) <- c("All", 'Absent', 'Present')
print(Table)

# compare outcome distribution by thymoma
mod1 <- nnet::multinom(outcome ~ Sex + Age.TP, data = data, subset = index)
mod2 <- nnet::multinom(outcome ~ TP.fp + Sex + Age.TP, data = data, subset = index)
anova(mod1, mod2)


#################### compare mortality rate between thymoma presence
TP.ind <- ifelse(data$thymoma == '+', 1, 0)
mod <- glm(death.ind ~ thymoma + Sex + Age.TP, subset = index, data = data, family = 'binomial')

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

# logistic regression
death.ind <- ifelse(data$outcome == 'Death', 1, 0)
mod <- glm(death.ind ~ thymoma, subset = index, data = data, family = 'binomial')

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


#################### distribution of giant cells and granuloma
Count <- with(data, rbind(table(Giant2[ind.biopsy == 1]), 
                          table(Granuloma[ind.biopsy == 1]), 
                          table(Eosinophils[ind.biopsy == 1])))
Prop <- Count / rowSums(Count)
Summary <- cbind(Count[,2], Prop[,2])
rownames(Summary) <- c("Giant Cells", "Granuloma", "Eosinophils")
colnames(Summary) <- c("Positive Count", "Positive proportion")
print(Summary)

Summary <- with(data, table(Giant.cells[ind.biopsy == 1], 
                            Granuloma[ind.biopsy == 1]))
print(Summary)

# logistic regression
Giant.cells.ind <- ifelse(data$Giant.cells == 'Pos', 1, 0)
mod <- glm(Giant.cells.ind ~ Granuloma + Age + Sex, family = 'binomial', 
           subset = which(ind.biopsy == 1), data = data,
           control = glm.control(maxit = 50))

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

# count of GG
table(data$GG)

#################### inflammatory cells between EOMG and LOMG
tt <- with(data, table(Inflam[ind.biopsy == 1], EOMG[ind.biopsy == 1]))
Count <- cbind(table(data$Inflam), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
Inflam.ind <- ifelse(data$Inflam == 'Pos', 1, 0)
mod <- glm(Inflam.ind ~ EOMG + Sex, subset = which(ind.biopsy == 1), data = data, family = 'binomial')

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

#################### inflammatory locations between EOMG and LOMG
vec <- with(data, apply(apply(cbind(Em, Pm, Pv, Is), 2, function(x) ifelse(x == '+', 1, 0)), 2, sum))
Base <- with(data, which(Em == '+' | Pm == '+' | Pv == '+' | Is == '+'))
vec <- c(vec, length(Base))
names(vec)[5] <- 'nBase'
print(vec)


# logistic regression (Em)
Em.ind <- ifelse(data$Em == '+', 1, 0)
mod <- glm(Em.ind ~ EOMG + Sex, family = 'binomial', subset = Base, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# logistic regression (Pm)
Pm.ind <- ifelse(data$Pm == '+', 1, 0)
mod <- glm(Pm.ind ~ EOMG + Sex, family = 'binomial', subset = Base, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# logistic regression (Pv)
Pv.ind <- ifelse(data$Pv == '+', 1, 0)
mod <- glm(Pv.ind ~ EOMG + Sex, family = 'binomial', subset = Base, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# logistic regression (Is)
Is.ind <- ifelse(data$Is == '+', 1, 0)
mod <- glm(Is.ind ~ EOMG + Sex, family = 'binomial', subset = Base, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')



#################### inflammatory locations by myositis categories

inflam.location <- with(data, apply(cbind(Em, Pm, Pv, Is), 2, function(x) ifelse(x == '+', 1, 0)))
Base <- with(data, which(Em == '+' | Pm == '+' | Pv == '+' | Is == '+'))
types <- na.omit(unique(data$myositis.subtype))
counts <- c()
nBases <- c()
for (j in 1:5) {
  idx.j <- which(data$myositis.subtype == types[j])
  vec.j <- apply(inflam.location[idx.j, ], 2, sum)
  nBases.j <- sum(apply(inflam.location[idx.j, ], 1, function(x) ifelse(max(x) == 1, 1, 0)))
  counts <- rbind(counts, vec.j)
  nBases <- c(nBases, nBases.j)
}
Prop <- c()
for (j in 1:4) {
  Prop <- rbind(Prop, counts[j, ] / nBases[j])
}
Prop <- rbind(Prop, 0)
colnames(counts) <- colnames(Prop) <- c("Em+", "Pm+", "Pv+", "Is+")
Table <- cbind(counts, Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = ncol(counts)))
print(Table)

# logistic regression
subset <- intersect(Base, which(subtype.MY %in% c('Polymyositis', 'GCM/GrM')))

# logistic regression (Em)
Em.ind <- ifelse(data$Em == '+', 1, 0)
mod <- glm(Em.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# sensitivity analysis
mod <- glm(Em.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', 
           subset = setdiff(subset, c(55, 118, 181)), data = data)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# logistic regression (Pm)
Pm.ind <- ifelse(data$Pm == '+', 1, 0)
mod <- glm(Pm.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', subset = subset, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# sensitivity analysis
mod <- glm(Pm.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', 
           subset = setdiff(subset, c(55, 118, 181)), data = data)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# logistic regression (Pv)
Pv.ind <- ifelse(data$Pv == '+', 1, 0)
mod <- glm(Pv.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', subset = Base, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# sensitivity analysis
mod <- glm(Pv.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', 
           subset = setdiff(subset, c(55, 118, 181)), data = data)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)


# logistic regression (Is)
Is.ind <- ifelse(data$Is == '+', 1, 0)
mod <- glm(Is.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', subset = subset, data = data)

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

# sensitivity analysis
mod <- glm(Is.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', 
           subset = setdiff(subset, c(55, 118, 181)), data = data)
# p-value
p <- coef(summary(mod))[2, 4]
print(p)

### Compare Group 1 (Dermatomyositis) and Group 2 (Polymyositis and GCM/GrM)

# logistic regression (Em)
mod <- glm(Em.ind ~ group1 + Sex + Age.MY, family = 'binomial', subset = Base, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# logistic regression (Pm)
Pm.ind <- ifelse(data$Pm == '+', 1, 0)
mod <- glm(Pm.ind ~ group1 + Sex + Age.MY, family = 'binomial', subset = Base, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# logistic regression (Pv)
Pv.ind <- ifelse(data$Pv == '+', 1, 0)
mod <- glm(Pv.ind ~ group1 + Sex + Age.MY, family = 'binomial', subset = Base, data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# logistic regression (Is)
Is.ind <- ifelse(data$Is == '+', 1, 0)
mod <- glm(Is.ind ~ group1 + Sex + Age.MY, family = 'binomial', subset = Base, data = data)

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


#################### Compare median CK level between EOMG and LOMG
mod <- rq(CK.level ~ EOMG + Sex, data = data, tau = 0.50)
# p-value
print(coef(summary(mod, se = 'ker'))[2, 4])

# Median difference of CK levels between EOMG and LOMG
coef <- coef(summary(mod, se = 'ker'))[2, 1]
se <- coef(summary(mod, se = 'ker'))[2, 2]
pp <- c(coef, coef - se * qnorm(.975), coef + se * qnorm(.975))
names(pp) <- c('Estimate', 'Lower bound', 'Upper bound')
print(round(pp, 2))

#################### O+/O- rates
Count <- table(data$Oplus)
Prop <- Count / sum(Count)
Summary <- rbind(Count, Prop)
rownames(Summary) <- c("Count","Proportion")
print(Summary)

#################### O+ by EOMG
tt <- with(data, table(Oplus, EOMG))
Count <- cbind(table(data$Oplus), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
rownames(Table)[-1] <- paste(colnames(tt))
print(Table)

# logistic regression
O.ind <- ifelse(data$Oplus == '+', 1, 0)
mod <- glm(O.ind ~ EOMG + Sex, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- -coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')



#################### O+ by thymoma
tt <- with(data, table(Oplus, TET))
Count <- cbind(rowSums(tt), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), round(Prop, 3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
rownames(Table)[-1] <- paste(colnames(tt))
print(Table)

# logistic regression
O.ind <- ifelse(data$Oplus == '+', 1, 0)
mod <- glm(O.ind ~ TET + Sex + Age.MG, family = 'binomial', data = data)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- -coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


#################### RNS
Count <- table(data$RNS)
Prop <- Count / sum(Count)
Summary <- rbind(Count, Prop)
rownames(Summary) <- c("Count","Proportion")
print(Summary)

# RNS by EOMG
tt <- with(data, table(RNS, EOMG))
Count <- cbind(table(MG3), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
rownames(Table)[-1] <- paste(colnames(tt))
print(Table)

# logistic regression
RNS.ind <- ifelse(data$RNS == '+', 1, 0)
mod <- glm(RNS.ind ~ EOMG + Sex, data = data, family = 'binomial')

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


#################### AchI
Count <- table(data$AchI)
Prop <- Count/sum(Count)
Summary <- rbind(Count, Prop)
rownames(Summary) <- c("Count","Proportion")
print(Summary)

# AchI by EOMG
tt <- with(data, table(AchI, EOMG))
Count <- cbind(table(MG3), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
rownames(Table)[-1] <- paste(colnames(tt))
print(Table)

# logistic regression
AchI.ind <- ifelse(data$AchI == '+', 1, 0)
mod <- glm(AchI.ind ~ EOMG + Sex, data = data, family = 'binomial')

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

#################### SFEMG
Count <- table(data$SFEMG)
Prop <- Count / sum(Count)
Summary <- rbind(Count, Prop)
rownames(Summary) <- c("Count", "Proportion")
print(Summary)

# SFEMG by EOMG
tt <- with(data, table(SFEMG, EOMG))
Count <- cbind(table(data$SFEMG), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
rownames(Table)[-1] <- paste(colnames(tt))
print(Table)

# logistic regression
SFEMG.ind <- ifelse(data$SFEMG == '+', 1, 0)
mod <- glm(SFEMG.ind ~ EOMG + Sex, data = data, family = 'binomial')

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


#################### mEMG
Count <- table(data$mEMG)
Prop <- Count / sum(Count)
Summary <- rbind(Count, Prop)
rownames(Summary) <- c("Count","Proportion")
print(Summary)

# mEMG by EOMG
tt <- with(data, table(mEMG, EOMG))
Count <- cbind(table(data$mEMG), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), round(Prop,3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
rownames(Table)[-1] <- paste(colnames(tt))
print(Table)


# logistic regression
mEMG.ind <- ifelse(data$mEMG == '+', 1, 0)
mod <- glm(mEMG.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# logistic regression
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


#################### EMG
Count <- table(data$EMG)
Prop <- Count / sum(Count)
Summary <- rbind(Count, Prop)
rownames(Summary) <- c("Count","Proportion")
print(Summary)

# EMG by EOMG
tt <- with(data, table(EMG, EOMG))
Count <- cbind(table(data$EMG), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), round(Prop, 3))
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
rownames(Table)[-1] <- paste(colnames(tt))
print(Table)

# logistic regression
EMG.ind <- ifelse(data$EMG == '+', 1, 0)
mod <- glm(EMG.ind ~ EOMG + Sex, data = data, family = 'binomial')

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



#################### myositis subtypes
Count <- table(data$myositis.subtype)
Prop <- Count / sum(Count)
Table <- cbind(Count, Prop)
tt <- with(data, table(myositis.subtype, Cardiac))
Count.C <- tt[, 2]
Prop.C <- Count.C / (tt[, 2] + tt[, 1])
Table2 <- cbind(Table, Count.C, Prop.C)
print(Table2)


#################### thymoma and myositis subtypes
tt <- with(data, table(myositis.subtype, thymoma))
Count.T <- tt[,2]
Prop.T <- Count.T / (tt[, 2] + tt[, 1])
Table3 <- cbind(Table, tt, Prop.T)
print(Table3)

#################### thymoma between group 1 (Dermatomyositis) and group 2 (Polymyositis and GCM/GrM)
tt <- with(data, table(group1, thymoma))
table <- cbind(tt, tt / rowSums(tt))
colnames(table) <- c('- (C)', '+ (C)', '- (P)', "+ (P)")
print(table)

# logistic regression
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM'))
thymoma.ind <- ifelse(data$thymoma == '+', 1, 0)
mod <- glm(thymoma.ind ~ myositis.subtype + Sex + Age.MY, subset = subset, 
           data = data, family = 'binomial')

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

# sensitivity analysis
mod <- glm(thymoma.ind ~ myositis.subtype + Sex + Age.MY, 
           subset = setdiff(subset, c(55, 118, 181)), 
           data = data, family = 'binomial')
p <- coef(summary(mod))[2, 4]
print(p)

# logistic regression
# Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
thymoma.ind <- ifelse(data$thymoma == '+', 1, 0)
mod <- glm(thymoma.ind ~ group2 + Sex + Age.MY, data = data, family = 'binomial')

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

# sensitivity analysis
mod <- glm(thymoma.ind ~ group2 + Sex + Age.MY, 
           subset = (1:206)[-c(55, 118, 181)], data = data, family = 'binomial')
p <- coef(summary(mod))[2, 4]
print(p)

#################### myocarditis and thymoma
tt <- with(data, table(thymoma, myocarditis))
Count <- cbind(table(data$EOMG), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
thymoma.ind <- ifelse(data$thymoma == '+', 1, 0)
mod <- glm(thymoma.ind ~ myocarditis + Sex + Age, data = data, family = 'binomial')

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

#################### myocarditis and cardiac involvement
tt <- with(data, table(thymoma, Cardiac))
Count <- cbind(table(data$Cardiac), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] <- "All"
print(Table)

# logistic regression
thymoma.ind <- ifelse(data$thymoma == '+', 1, 0)
mod <- glm(thymoma.ind ~ Cardiac + Sex + Age, data = data, family = 'binomial')

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


#################### O+ and myositis categories
tt <- with(data, table(myositis.subtype, Oplus))
Prop <- tt[, 2] / rowSums(tt)
Table <- cbind(tt, Prop)
print(Table)

# logistic regression
# compare O+ between Polymyositis and GCM/GrM
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM'))
O.ind <- ifelse(data$Oplus == '+', 1, 0)
mod <- glm(O.ind ~ subtype.MY + Sex + Age.MG, subset = subset, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')


# compare between Group 1 (Dermatomyositis) and Group 2 (Polymyositis and GCM/GrM)
mod <- glm(O.ind ~ group1 + Sex + Age.MG, data = data, family = 'binomial')

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




#################### EOMG and myositis categories
tt <- with(data, table(myositis.subtype, EOMG))
Count <- cbind(table(data$myositis.subtype), tt)
Total <- colSums(Count)
Prop <- t(Count) / Total
Count <- matrix(as.character(Count), nrow = 3, byrow = T)
Prop <- matrix(as.character(round(Prop, 3)), nrow = 3, byrow = F)
Table <- rbind(Count, Prop)
colnames(Table) <- rownames(tt)
rownames(Table) <- rep(c("All", "EOMG", "LOMG"), 2)
print(Table)

# likelihood ratio test 
mod1 <- nnet::multinom(subtype.MY ~ Sex)
mod2 <- nnet::multinom(subtype.MY ~ EOMG + Sex)
anova(mod1, mod2)

#################### EOMG within each myositis category
tt <- with(data, table(myositis.subtype, EOMG))
Prop <- tt / rowSums(tt)
Table <- cbind(tt, Prop)
Table

# logistic regression
# compare EOMG between Group 1 (Dermatomyositis) and Group 2 (Polymyositis and GCM/GrM)
EOMG.ind <- ifelse(EOMG == 'EOMG', 1, 0)
mod <- glm(EOMG.ind ~ group1 + Sex + Age.MG, data = data, family = 'binomial')

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

# logistic regression
# compare EOMG between between polymyositis and GCM/GrM
subset <- which(data$myositis.subtype %in% c('GCM/GrM', 'Polymyositis')) 
mod <- glm(EOMG.ind ~ myositis.subtype + Sex + Age.MG, subset = subset, data = data, family = 'binomial')

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


# logistic regression
# compare EOMG between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
mod <- glm(EOMG.ind ~ group2 + Sex + Age.MG, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

#################### Cardiac involvement within each myositis category
tt <- with(data, table(myositis.subtype, Cardiac))
Prop <- tt / rowSums(tt)
Table <- cbind(tt, Prop)
print(Table)


# logistic regression
# compare EOMG between between polymyositis and GCM/GrM
subset <- which(data$myositis.subtype %in% c('GCM/GrM', 'Polymyositis'))
Cardiac.ind <- ifelse(data$Cardiac == 'Yes', 1, 0)
mod <- glm(Cardiac.ind ~ myositis.subtype + Sex + Age.MY, subset = subset, data = data, family = 'binomial')

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

# sensitivity analysis
mod <- glm(Cardiac.ind ~ myositis.subtype + Sex + Age.MY, subset = setdiff(subset, c(55, 118, 181)), 
           data = data, family = 'binomial')
p <- coef(summary(mod))[2, 4]
print(p)

# logistic regression
# compare EOMG between between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
mod <- glm(Cardiac ~ group2 + Sex + Age.MY, data = data, family = 'binomial')

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

# logistic regression
# compare EOMG between between Group 1 (Dermatomyositis) and Group 2 (Polymyositis and GCM/GrM)
mod <- glm(Cardiac ~ group1 + Sex + Age.MY, data = data, family = 'binomial')

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


#################### AChR within myositis subtypes
tt <- with(data, table(myositis.subtype, AChR))
Prop <- tt / rowSums(tt)
Table <- cbind(tt, Prop)
print(Table)


#################### median of AChR within myositis subtypes
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
Count <- with(data, tapply(AChR.Level, myositis.subtype, function(x) sum(!is.na(x))))
median_by_subtype <- tapply(AChR.Level, myositis.subtype, median, na.rm = TRUE)
iqr_by_subtype <- tapply(AChR.Level, subtype.MY, iqr)
sd_by_subtype <- tapply(AChR.Level, subtype.MY, sd, na.rm=TRUE)
Table <- cbind(Count, round(median_by_subtype, 2), round(sd_by_subtype, 2), iqr_by_subtype)
rownames(Table) <- names(Count)
colnames(Table) <- c("Count", "Median_AchR", 'SD', 'IQR')

# median regression
# Compare the AChR titer between Group 1 (Dermatomyositis) with Group 2 (Polymyositis and GCM/GrM)
mod <- rq(AChR.Level ~ group2 + Age.MY + Sex, data = data)
# p-value
coef(summary(mod, se = 'ker'))[2, 4]

# median in Group 1 and Group 2
with(data, tapply(AChR.Level, group2, function(x) median(x, na.rm = TRUE)))


#################### striational antibody within myositis subtype
tt <- with(data, table(myositis.subtype, Striated))
Prop <- tt / rowSums(tt)
Table <- cbind(tt, Prop)
print(Table)

# logistic regression
# between polymyositis and GCM
as.ind <- ifelse(data$Striated == '+', 1, 0)
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM'))
mod <- glm(as.ind ~ subtype.MY + Sex + Age.MG, subset = subset, family = 'binomial')
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

# logistic regression
# between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
as.ind <- ifelse(data$Striated == '+', 1, 0)
mod <- glm(as.ind ~ group2 + Sex + Age.MG, data = data, family = 'binomial')

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


#################### MSA within myositis subtype
tt <- with(data, table(myositis.subtype, MSA))
Prop <- tt / rowSums(tt)
Table <- cbind(tt, Prop)
print(Table)

# logistic regression
# between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
msa.ind <- ifelse(data$MSA == '+', 1, 0)
mod <- glm(msa.ind ~ group2 + Sex + Age.MY, data = data, family = 'binomial')

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
knitr::kable(round(odds, 2), caption="Odds ratio with 95% confidence intervals (adjusted for gender and age)")


# logistic regression
# between Group 1 (Dermatomyositis) and Group 2 (Polymyositis and GCM/GrM)
mod <- glm(msa.ind ~ group1 + Sex + Age.MY, data = data, family = 'binomial')

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


#################### outcome within each myositis subtype (after breakpoint)
subset <- which(data$year >= 1993)
# Likelihood ratio tests
subset2 <- intersect(which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM')), subset)
mod1 <- nnet::multinom(outcome ~ Sex + Age.max, subset = subset2, data = data)
mod2 <- nnet::multinom(outcome ~ myositis.subtype + Sex + Age.max, subset = subset2, data = data)
anova(mod1, mod2)


# death
tt <- table(death[subset], subtype.MY[subset])
Count <- cbind(table(death[subset]), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- rbind(t(Count), round(Prop,3))
rownames(Table)[c(1, 5 + 2)] = "All"
print(Table)

# likelihood ratio test
death.ind <- ifelse(outcome == 'Death', 1, 0)
subset3 <- intersect(subset, which(is.na(data$myositis.subtype) == FALSE))
mod1 <- glm(death.ind ~ Sex + Age.max, family = 'binomial', data = data, subset = subset3)
mod2 <- glm(death.ind ~ myositis.subtype + Sex + Age.max, family = 'binomial', data = data, subset = subset3)
anova(mod1, mod2)

# logistic regression
# between Polymyositis and GCM/GrM 
mod <- glm(death.ind ~ myositis.subtype + Sex + Age.max, family = 'binomial', data = data, subset = subset2)

# p-value
p <- coef(summary(mod))[2, 4]
print(p)

# odds ratio
odds <- -coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

# logistic regression
# between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
subset4 <- intersect(subset, which(is.na(data$group1) == FALSE))
mod1 <- nnet::multinom(outcome ~ Sex + Age.max, subset = subset4, data = data)
mod2 <- nnet::multinom(outcome ~ group2 + Sex + Age.max, subset = subset4, data = data)
anova(mod1, mod2)


# logistic regression
# mortality rate between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
mod <- glm(death.ind ~ group2 + Sex + Age.max, subset = subset, data = data, family = 'binomial')

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

# logistic regression
# mortality rate between Group 1 (polymyositis and GCM/GrM) and Group 2 (Dermatomyositis)
mod <- glm(death.ind ~ group1 + Sex + Age.max, subset = subset, data = data, family = 'binomial')

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


#################### Corticosteriods within each myositis subtype (after breakpoint)
subset <- which(data$year >= 1993)
tab <- with(data, t(table(cort[subset], myositis.subtype[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[6] <- 'All'
tab[6, 3:4] <- tab[6, 1:2] / sum(tab[6, 1:2])
print(tab)

# logistic regression
cort.ind <- ifelse(data$cort == '+', 1, 0)
data$myositis.subtype <- as.factor(data$myositis.subtype)
mod <- glm(cort.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', subset = subset, data = data, control = glm.control(maxit = 50))

# pairwise comparison of myositis
mod.glht <- glht(mod, linfct = mcp("myositis.subtype" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)

# p-values
print(round(pvs, 4))


# logistic regression
# between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
cort.ind <- ifelse(data$cort == '+', 1, 0)
mod <- glm(cort.ind ~ group2 + Sex + Age.MY, family = 'binomial', subset = subset, data = data, control = glm.control(maxit = 50))

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



#################### Plasmapheresis within each myositis subtype (after breakpoint)
subset <- which(data$year >= 1993)
tab <- with(data, t(table(plasmapheresis[subset], myositis.subtype[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[6] <- 'All'
tab[6, 3:4] <- tab[6, 1:2] / sum(tab[6, 1:2])
print(tab)

# logistic regression
plasmapheresis.ind <- ifelse(data$plasmapheresis == '+', 1, 0)
data$myositis.subtype <- as.factor(data$myositis.subtype)
mod <- glm(plasmapheresis.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', subset = subset, data = data, control = glm.control(maxit = 50))

# pairwise comparison of myositis
mod.glht <- glht(mod, linfct = mcp("myositis.subtype" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)

# p-values
print(round(pvs, 4))


# logistic regression
# between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
plasmapheresis.ind <- ifelse(data$plasmapheresis == '+', 1, 0)
mod <- glm(plasmapheresis.ind ~ group2 + Sex + Age.MY, family = 'binomial', subset = subset, data = data, control = glm.control(maxit = 50))

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






#################### Respiratory support within each myositis subtype (after breakpoint)
subset <- which(data$year >= 1993)
tab <- with(data, t(table(rs[subset], myositis.subtype[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[6] <- 'All'
tab[6, 3:4] <- tab[6, 1:2] / sum(tab[6, 1:2])
print(tab)

# logistic regression
rs.ind <- ifelse(data$rs == '+', 1, 0)
data$myositis.subtype <- as.factor(data$myositis.subtype)
mod <- glm(rs.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', subset = subset, data = data, control = glm.control(maxit = 50))

# pairwise comparison of myositis
mod.glht <- glht(mod, linfct = mcp("myositis.subtype" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)

# p-values
print(round(pvs, 4))


# logistic regression
# between Group 1 (polymyositis and GCM/GrM) and Group 2 (OM, XOther, and Dermatomyositis)
rs.ind <- ifelse(data$rs == '+', 1, 0)
mod <- glm(rs.ind ~ group2 + Sex + Age.MY, family = 'binomial', subset = subset, data = data, control = glm.control(maxit = 50))

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


#################### Cardiovascular treatments within each myositis subtype (after breakpoint)
subset <- which(data$year >= 1993)
tab <- with(data, t(table(ct[subset], myositis.subtype[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[6] <- 'All'
tab[6, 3:4] <- tab[6, 1:2] / sum(tab[6, 1:2])
print(tab)

# logistic regression
# compare 'GCM/GrM' with 'Polymyositis'
ct.ind <- ifelse(data$ct == '+', 1, 0)
subset2 <- which(data$myositis.subtype %in% c('GCM/GrM', 'Polymyositis') & data$year >= 1993)
mod <- glm(ct.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', 
           subset = subset2, data = data, control = glm.control(maxit = 50))

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

#################### Non-steroid immunomodulators within each myositis subtype (after breakpoint)
subset <- which(data$year >= 1993)
tab <- with(data, t(table(imm[subset], myositis.subtype[subset])))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[6] <- 'All'
tab[6, 3:4] <- tab[6, 1:2] / sum(tab[6, 1:2])

# logistic regression
imm.ind <- ifelse(data$imm == '+', 1, 0)
data$myositis.subtype <- as.factor(data$myositis.subtype)
mod <- glm(imm.ind ~ myositis.subtype + Sex + Age.MY, family = 'binomial', 
           subset = subset, data = data, control = glm.control(maxit = 50))

# multiple comparison
mod.glht <- glht(mod, linfct = mcp("myositis.subtype" = "Tukey"))   
summ <- summary(mod.glht)
pvs <- summ$test$pvalues
names(pvs) <- names(summ$test$coefficients)

# p-value
print(round(pvs, 4))


#################### MG subtypes within myositis subtype
tt <- with(data, table(myositis.subtype, MG.subtype))
Table <- tt / rowSums(tt)
print(Table)


# likelihood ratio test
# compare MG Subtype between polymyositis and Dermatomyositis
# the significance level should be divided by 3 to adjust for multiple comparison
subset1 <- which(data$myositis.subtype %in% c('Polymyositis', 'Dermatomyositis') & data$year >= 1993)
mod1 <- nnet::multinom(MG.subtype ~ Sex + Age.MG, subset = subset1, data = data)
mod2 <- nnet::multinom(MG.subtype ~ myositis.subtype + Sex + Age.MG, subset = subset1, data = data)
anova(mod1, mod2)

# likelihood ratio test
# compare MG Subtype between GCM/GrM and Dermatomyositis
# the significance level should be divided by 3 to adjust for multiple comparison
subset2 <- which(data$myositis.subtype %in% c('GCM/GrM', 'Dermatomyositis') & data$year >= 1993)
mod1 <- nnet::multinom(MG.subtype ~ Sex + Age.MG,  subset = subset2, data = data)
mod2 <- nnet::multinom(MG.subtype ~ myositis.subtype + Sex + Age.MG, subset = subset2, data = data)
anova(mod1, mod2)

# likelihood ratio test
# compare MG Subtype between GCM/GrM and Polymyositis
# the significance level should be divided by 3 to adjust for multiple comparison
subset3 <- which(data$myositis.subtype %in% c('GCM/GrM', 'Polymyositis') & data$year >= 1993)
mod1 <- nnet::multinom(MG.subtype ~ Sex + Age.MG, subset = subset3, data = data)
mod2 <- nnet::multinom(MG.subtype ~ myositis.subtype + Sex + Age.MG,  subset = subset3, data = data)
anova(mod1, mod2)


#################### temporal sequence within each myositis subtype
tt <- with(data, table(Earlier.onset, myositis.subtype))
Count <- cbind(table(Y), tt)
Total <- colSums(Count)
Prop <- t(Count)/Total
Table <- cbind(t(Count), Prop)
colnames(Table) <- paste(colnames(Table), rep(c("(C)","(P)"), each = nrow(tt)))
rownames(Table)[1] = "All"
print(Table)

# likelihood ratio test
# Compare between Group 1 (Dermatomyositis) with Group 2 (polymyositis and GCM)
Age.mmyg <- with(data, apply(cbind(Age.MG, Age.MY), 1, max))
subset <- which(is.na(data$group1) == FALSE & data$year >= 1993)
mod1 <- nnet::multinom(Y ~ Sex + Age.mmyg, subset = subset, data = data)
mod2 <- nnet::multinom(Y ~ group1 + Sex + Age.mmyg, subset = subset, data = data)
anova(mod1, mod2)


#################### CD-4 between EOMG and LOMG
tt <- with(data, table(EOMG, CD4))
Table <- cbind(tt, tt / c(rowSums(tt)[1], rowSums(tt)[2],
                          rowSums(tt)[1], rowSums(tt)[2]))
print(Table)

# logistic regression
CD4.ind <- ifelse(data$CD4 == '+', 1, 0)
mod <- glm(CD4.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')

#################### CD-4 among myositis subtypes
tt <- with(data, table(myositis.subtype, CD4))
Table <- cbind(tt, apply(tt, 2, function(x) x / rowSums(tt)))
print(Table)


#################### CD-8 between EOMG and LOMG
tt <- with(data, table(EOMG, CD8))
Table <- cbind(tt, apply(tt, 2, function(x) x / rowSums(tt)))
print(Table)

# logistic regression
CD8.ind <- ifelse(data$CD8 == '+', 1, 0)
mod <- glm(CD8.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### CD-8 among myositis subtypes
tt <- with(data, table(myositis.subtype, CD8))
Table <- cbind(tt, apply(tt, 2, function(x) x / rowSums(tt)))
print(Table)

#################### B-plus between EOMG and LOMG
tt <- table(EOMG, Bp)
Table <- cbind(tt, tt / c(rowSums(tt)[1], rowSums(tt)[2],
                          rowSums(tt)[1], rowSums(tt)[2]))
print(Table)

# logistic regression
Bp.ind <- ifelse(data$Bp == '+', 1, 0)
mod <- glm(Bp.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### B-plus among myositis subtypes
tt <- with(data, table(myositis.subtype, Bplus))
Table <- cbind(tt, apply(tt, 2, function(x) x / rowSums(tt)))
print(Table)


#################### Compare CD8 and CD4
Count <- table(data$CD.comp)
Prop <- Count / sum(Count)
Summary <- cbind(Count, Prop)
rownames(Summary) <- c("CD8>=CD4", "CD4>CD8")
print(Summary)

#################### Compare CD8 and CD4 between EOMG and LOMG
tt <- with(data, table(EOMG, CD.comp))
Table <- cbind(tt, tt / c(rowSums(tt)[1], rowSums(tt)[2],
                          rowSums(tt)[1], rowSums(tt)[2]))

# logistic regression
CD.comp.ind <- ifelse(data$CD.comp == '+', 1, 0)
mod <- glm(CD.comp.ind ~ EOMG + Sex, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Compare CD8 and CD4 among myositis subtypes
tt <- with(data, table(myositis.subtype, CD.comp))
bbb <- apply(tt, 1, function(x) if (max(x) == 0) 0 else x / sum(x))
bbb <- do.call('rbind', bbb)
Table <- cbind(tt, bbb)
print(Table)


#################### Macrophages between EOMG and LOMG
tt <- with(data, table(EOMG, macrophages))
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))

# logistic regression
macrophages.ind <- ifelse(data$macrophages == '+', 1, 0)
mod <- glm(macrophages.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Macrophages among myositis subtypes
tt <- with(data, table(myositis.subtype, macrophages))
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
print(Table)

#################### Plasma between EOMG and LOMG
tt <- table(EOMG, Plasma)  

# logistic regression
Plasma.ind <- ifelse(data$Plasma == '+', 1, 0)
mod <- glm(Plasma.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Plasma among myositis subtypes
tt <- with(data, table(myositis.subtype, Plasma))
print(tt)


#################### MHC-I between EOMG and LOMG
tt <- with(data, table(EOMG, MHC1))
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
print(Table)

# logistic regression
MHC1.ind <- ifelse(data$MHC1 == '+', 1, 0)
mod <- glm(MHC1.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### MHC-I among myositis subtypes
tt <- with(data, table(myositis.subtype, MHC1))
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
print(Table)


#################### MHC-II between EOMG and LOMG
tt <- with(data, table(EOMG, MHC2))
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
print(Table)

#################### MHC-II among myositis subtypes
tt <- with(data, table(myositis.subtype, MHC2))
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
print(Table)


#################### IgG between EOMG and LOMG
tt <- with(data, table(EOMG, IgG))
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
print(Table)

# logistic regression
complement.ind <- ifelse(data$IgG == '+', 1, 0)
mod <- glm(complement.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### IgG among myositis subtypes
tt <- with(data, table(myositis.subtype, IgG))
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
print(Table)

#################### Lag time within myositis subtypes
XX <- with(data, data.frame(Lag = lag, Subtype = group1))

summary_table <- XX %>%
  group_by(Subtype) %>%
  summarise(
    Count = n(),
    Mean = mean(Lag, na.rm = TRUE),
    Median = median(Lag, na.rm = TRUE),
    SD = sd(Lag, na.rm = TRUE),
    Min = min(Lag, na.rm = TRUE),
    Max = max(Lag, na.rm = TRUE),
    Q.25 = quantile(Lag, probs=0.25, na.rm=TRUE),
    Q.75 = quantile(Lag, probs=0.75, na.rm=TRUE)
  )
print(summary_table)

# poisson regression
Age.mmyg <- with(data, apply(cbind(Age.MG, Age.MY), 1, max))
subset <- (1:206)[-19]
mod <- glm(lag ~ group1 + Sex + Age.mmyg, subset = subset, data = data,
           family = poisson(link = "log"))

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
coef <- -coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
coef <- cbind(coef, coef - qnorm(.975) * ses, coef + qnorm(.975) * ses) 
colnames(coef) <- c('Estimate', 'Lower bound', 'Upper bound')
print(coef)

#################### Giant cells and myocarditis
subset <- which(data$GG != 'nobiopsy')
Summary <- with(data, table(GG[subset], myocarditis[subset]))
print(Summary)

# logistic regression
GG.ind <- ifelse(data$GG == '+', 1, 0)
mod <- glm(GG.ind ~ myocarditis + Sex + Age, subset = subset, family = 'binomial', data = data)

# p-value
coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Giant cells and myocarditis
#################### interaction effects between Giant cells and thymoma
subset <- which(data$GG != 'nobiopsy' & data$myocarditis_subset == 1)
myocarditis.ind <- ifelse(data$myocarditis == 'Myocarditis', 1, 0)
mod <- glm(myocarditis.ind ~ thymoma * GG + Sex + Age, family = 'binomial',
           subset = subset,
           data = data, control = glm.control(maxit = 50))

# p-value
coef(summary(mod))[c(2, 3, 6), 4]

# odds ratio
odds <- coef(summary(mod))[c(2, 3, 6), 1]
ses <- coef(summary(mod))[c(2, 3, 6), 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Eosinophils and EOMG
tt <- with(data, table(EOMG, Eosinophils))
tt <- rbind(colSums(tt), tt)
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
colnames(Table) <- c('Eosinophils - (C)', 'Eosinophils + (C)', 
                     'Eosinophils - (P)', 'Eosinophils + (P)')
rownames(Table)[1] <- c('All')
print(Table)

# logistic regression
Eosinophils.ind <- ifelse(data$Eosinophils == '+', 1, 0)
mod <- glm(Eosinophils.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
print(coef(summary(mod))[2, 4])

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Eosinophils and EOMG
tt <- with(data, table(myositis.subtype, Eosinophils))
tt <- rbind(colSums(tt), tt)
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
colnames(Table) <- c('Eosinophils - (C)', 'Eosinophils + (C)', 
                     'Eosinophils - (P)', 'Eosinophils + (P)')
rownames(Table)[1] <- c('All')
print(Table)

# comparison among myositis subcategories
subset <- which(is.na(data$myositis.subtype) == FALSE)
# likelihood ratio test
# logistic regression
mod1 <- glm(Eosinophils.ind ~ Sex + Age.MY, family = 'binomial', 
            control = glm.control(maxit = 50), subset = subset)
mod2 <- glm(Eosinophils.ind ~ subtype.MY + Sex + Age.MY, family = 'binomial', 
            control = glm.control(maxit = 50))
anova(mod1, mod2)

#################### Eosinophils between polymyositis and GCM
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM'))

# logistic regression
Eosinophils.ind <- ifelse(data$Eosinophils == "+", 1, 0)
mod <- glm(Eosinophils.ind ~ myositis.subtype + Sex + Age.MY, 
           subset = subset, data = data, family = 'binomial', 
           control = glm.control(maxit = 50))

# p-value
coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

#################### Focal Necrosis between EOMG and LOMG
tt <- with(data, table(EOMG, FN))
tt <- rbind(colSums(tt), tt)
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
colnames(Table) <- c('focal necrosis - (C)', 'focal necrosis + (C)', 
                     'focal necrosis - (P)', 'focal necrosis + (P)')
rownames(Table)[1] <- c('All')
print(Table)

# logistic regression
FN.ind <- ifelse(data$FN == "+", 1, 0)
mod <- glm(FN.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Focal Necrosis among myositis subtypes
tt <- with(data, table(myositis.subtype, FN))
tt <- rbind(colSums(tt), tt)
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
colnames(Table) <- c('focal necrosis - (C)', 'focal necrosis + (C)', 
                     'focal necrosis - (P)', 'focal necrosis + (P)')
rownames(Table)[1] <- c('All')
print(Table)

# compare between Polymyositis and GCM/GrM
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM'))

# logistic regression
FN.ind <- ifelse(data$FN == "+", 1, 0)
mod <- glm(FN.ind ~ myositis.subtype + Sex + Age.MY, data = data, 
           subset = subset, family = 'binomial')

# p-value
coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

# compare between Group 2 (determatomyositis) with Group 1 (polymyositis and GCM)
mod <- glm(FN.ind ~ group1 + Sex + Age, data = data, family = 'binomial')

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

#################### Regeneration and EOMG
tt <- with(data, table(EOMG, Regeneration))
tt <- rbind(colSums(tt), tt)
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
colnames(Table) <- c('regeneration - (C)', 'regeneration + (C)', 
                     'regeneration - (P)', 'regeneration + (P)')
rownames(Table)[1] <- c('All')

# logistic regression
Regeneration.ind <- ifelse(data$Regeneration == "+", 1, 0)
mod <- glm(Regeneration.ind ~ EOMG + Sex, data = data, family = 'binomial')

# p-value
coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds) 
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)


#################### Regeneration among myositis subtypes
tt <- with(data, table(myositis.subtype, Regeneration))
tt <- rbind(colSums(tt), tt)
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
colnames(Table) <- c('regeneration - (C)', 'regeneration + (C)', 
                     'regeneration - (P)', 'regeneration + (P)')
rownames(Table)[1] <- c('All')
print(Table)


# compare between Polymyositis and GCM/GrM
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM'))

# logistic regression
Regeneration.ind <- ifelse(data$Regeneration == "+", 1, 0)
mod <- glm(Regeneration.ind ~ myositis.subtype + Sex + Age.MY, data = data, 
           subset = subset, family = 'binomial')

# p-value
coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

# compare between Group 2 (determatomyositis) with Group 1 (polymyositis and GCM)
mod <- glm(Regeneration.ind ~ group1 + Sex + Age, data = data, family = 'binomial')

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


#################### Increased Fiber Size Variability among myositis subtypes
tt <- with(data, table(myositis.subtype, IFSV))
tt <- rbind(colSums(tt), tt)
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
colnames(Table) <- c('regeneration - (C)', 'regeneration + (C)', 
                     'regeneration - (P)', 'regeneration + (P)')
rownames(Table)[1] <- c('All')
print(Table)


# compare between Polymyositis and GCM/GrM
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM'))

# logistic regression
IFSV.ind <- ifelse(data$IFSV == "+", 1, 0)
mod <- glm(IFSV.ind ~ myositis.subtype + Sex + Age.MY, data = data, 
           subset = subset, family = 'binomial')

# p-value
coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

# compare between Group 2 (determatomyositis) with Group 1 (polymyositis and GCM)
mod <- glm(IFSV.ind ~ group1 + Sex + Age, data = data, family = 'binomial')

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


#################### Increased Connective Tissue among myositis subtypes
tt <- with(data, table(myositis.subtype, ICT))
tt <- rbind(colSums(tt), tt)
Table <- cbind(tt, apply(tt[, 1:2], 2, function(x) x / rowSums(tt[, 1:2])))
colnames(Table) <- c('regeneration - (C)', 'regeneration + (C)', 
                     'regeneration - (P)', 'regeneration + (P)')
rownames(Table)[1] <- c('All')
print(Table)


# compare between Polymyositis and GCM/GrM
subset <- which(data$myositis.subtype %in% c('Polymyositis', 'GCM/GrM'))

# logistic regression
ICT.ind <- ifelse(data$ICT == "+", 1, 0)
mod <- glm(ICT.ind ~ myositis.subtype + Sex + Age.MY, data = data, 
           subset = subset, family = 'binomial')

# p-value
coef(summary(mod))[2, 4]

# odds ratio
odds <- coef(summary(mod))[2, 1]
ses <- coef(summary(mod))[2, 2] 
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)
colnames(odds) <- c('Odds ratio', 'Lower bound', 'Upper bound')
print(odds)

# compare between Group 2 (determatomyositis) with Group 1 (polymyositis and GCM)
mod <- glm(ICT.ind ~ group1 + Sex + Age, data = data, family = 'binomial')

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
