##################### set working directory
setwd("/Users/jileilin/Desktop/Research/thymomal case reports/")

##################### load packages
library('quantreg')

#################### read data
data <- read.csv('data4.csv')[, -1]

#################### Age at presentation
Age <- data$Age
iqr <- function(x) paste0(round(quantile(na.omit(x), .25), 2), '-', 
                          round(quantile(na.omit(x), .75), 2))
tab <- rbind(c(round(median(Age, na.rm = TRUE), 2),
               iqr(Age)),
             c(round(median(Age[which(Sex == 'F')], na.rm = TRUE), 2),
               iqr(Age[which(Sex == 'F')])),
             c(round(median(Age[which(Sex == 'M')], na.rm = TRUE), 2),
               iqr(Age[which(Sex == 'M')])))
colnames(tab) <- c('Median', 'IQR')
rownames(tab) <- c('All', 'Female', 'Male')
tab <- as.data.frame(tab)
tab$Median <- as.numeric(tab$Median)
print(tab)

#################### Sex
Sex <- data$Sex
tab <- table(Sex)
names(tab) <- c('Female', 'Male')
print(tab)

#################### Cause leading to discovery of thymoma
tab <- with(data, table(Sex, Cause))
tab <- cbind(tab, tab / rowSums(tab)) 
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', 
                                                 ' (P)', ' (P)'))
rownames(tab) <- c('Female', 'Male', 'All')
print(tab)

#################### WHO Subtype
tab <- with(data, table(Sex, who.subtype))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- mapply(paste0, colnames(tab), c(' (C)', ' (C)', ' (C)', ' (C)', ' (C)', 
                                                 ' (P)', ' (P)', ' (P)', ' (P)', ' (P)'))
rownames(tab) <- c('Female', 'Male')
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 6:10] <- tab[3, 1:5] / sum(tab[3, 1:5])
print(tab)

#################### Clinical MG
tab <- with(data, t(table(Clinical.MG, Sex)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)


#################### Presence of AChR
tab <- with(data, t(table(AChR, Sex)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:2], ' (C)'),
                   paste0(colnames(tab)[1:2], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 3:4] <- tab[3, 1:2] / sum(tab[3, 1:2])
print(tab)

#################### Maximum dimension of Size (before) and Local Symptoms
bsize <- data$bsize
ls <- data$Local.Symptoms
tab <- data.frame(matrix(0, 3, 4))
colnames(tab) <- c('Mean', 'SD', 'Median', 'IQR')
tab$Mean <- c(round(mean(bsize, na.rm = TRUE), 2),
              round(mean(bsize[ls == 'Existent'], na.rm = TRUE), 2),
              round(mean(bsize[ls == 'None'], na.rm = TRUE), 2))
tab$SD <- c(round(sd(bsize, na.rm = TRUE), 2),
            round(sd(bsize[ls == 'Existent'], na.rm = TRUE), 2),
            round(sd(bsize[ls == 'None'], na.rm = TRUE), 2))
tab$Median <- c(round(median(bsize, na.rm = TRUE), 2),
                round(median(bsize[ls == 'Existent'], na.rm = TRUE), 2),
                round(median(bsize[ls == 'None'], na.rm = TRUE), 2))
tab$IQR <- c(iqr(bsize),
             iqr(bsize[ls == 'Existent']),
             iqr(bsize[ls == 'None']))
rownames(tab) <- c('All', 'Existent', 'None')
colnames(tab) <- c('Mean', 'SD', 'Median', 'IQR')
print(tab)

# ordinary least squares
mod <- lm(bsize ~ ls)
# p-value
print(coef(summary(mod))[2, 4])


# median regression
mod <- rq(bsize ~ ls, tau = 0.50)
# p-value
print(coef(summary(mod, se = 'ker'))[2, 4])



#################### Pathological Findings
index <- which(colnames(data) %in% c('Hemorrhage', 'Necrosis', 'Cyst',
                                     'Fibrosis', 'Pleural.Effusion'))
pf <- c()
for (i in index) {
  yi <- data[, i]
  yi <- ifelse(yi == '+', 1, 0)
  pf <- rbind(pf, c(sum(yi), mean(yi)))
}
rownames(pf) <- colnames(data)[index]
colnames(pf) <- c('Positive cases', 'Positive rate')
print(pf)


#################### Masaoka Stage
tab <- with(data, t(table(Masaoka.Stage, Sex)))
tab <- cbind(tab, tab / rowSums(tab)) 
colnames(tab) <- c(paste0(colnames(tab)[1:3], ' (C)'),
                   paste0(colnames(tab)[1:3], ' (P)'))
tab <- rbind(tab, apply(tab, 2, sum))
rownames(tab)[3] <- 'All'
tab[3, 4:6] <- tab[3, 1:3] / sum(tab[3, 1:3])
print(tab)


#################### Duration
duration <- data$duration
tab <- c(round(mean(duration, na.rm = TRUE), 2),
         round(sd(duration, na.rm = TRUE), 2),
         round(median(duration, na.rm = TRUE), 2),
         iqr(duration))
names(tab) <- c('Mean', 'SD', 'Median', 'IQR')
print(tab)

#################### Recurrence
recurrence <- data$recurrence
tab <- c(table(recurrence), table(recurrence) / sum(table(recurrence)))
names(tab) <- c('N (C)', 'Y (C)', 'N (P)', 'Y (P)')
print(tab)