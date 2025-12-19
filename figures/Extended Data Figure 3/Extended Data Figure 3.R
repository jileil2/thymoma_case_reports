##################### set working directory
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
subset <- which(data$non.derm.IBM == '-')
cancer.detailed <- data$cancer.detailed
ll <- c()
for (i in 1:n) {
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


#################### forest plot for sex
sex.ind <- ifelse(data$Sex == 'F', 1, 0)
mod <- glm(sex.ind ~ cancer.detailed + Age, family = 'binomial', 
           data = data, subset = subset)

# p-value
p <- coef(summary(mod))[2:6, 4]

# odds ratio
odds <- coef(summary(mod))[2:6, 1]
ses <- coef(summary(mod))[2:6, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)

# groupwise proportion
pp <- tapply(sex.ind[subset], cancer.detailed[subset], function(x) mean(x, na.rm = TRUE))[1:5]
nn <- tapply(sex.ind[subset], cancer.detailed[subset], function(x) sum(x, na.rm = TRUE))[1:5]
NN <- tapply(sex.ind[subset], cancer.detailed[subset], function(x) length(na.omit(x)))[1:5]
prop <- paste0(nn, '/', NN, ', ', round(pp * 100))

# forest plot
mean <- odds[, 1]
lower <- odds[, 2]
upper <- odds[, 3]
pvs <- p
pvs <- round(pvs, 4)
base_data <- tibble::tibble(mean  = mean, lower = lower, upper = upper,
                            OR = round(mean, 2), pvs = pvs,
                            prop = prop, 
                            variable = c('TET', 'Lung cancer', 
                                         'Melanoma', 'Gastrointestinal cancer',
                                         'Genitourinary cancer'))


base_data |>
  forestplot(labeltext = c(variable, prop, OR, pvs),
             align = c("l", "c", "c", "c"),     # <-- Align columns
             xlim = c(0, 12),
             title = expression(bold("Sex (overall female proportion: 176/489, 36%)")),
             xlab = expression(bold("Association with female versus male")),
             xlog = TRUE,
             boxsize = 0.35,
             lwd.ci = 5,                       # Thicker CI lines
             xticks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
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
  fp_add_header(variable = c("", "Cancer type\n"),
                OR = c("", "Adj-OR\n"),
                pvs = c("", "p-value\n"),
                prop = c("", "Incidence\n(n/N, %)")) |>
  fp_set_zebra_style("#EFEFEF")

#################### violin plot for age
Age <- data$Age
non.derm.IBM <- data$non.derm.IBM
subset <- which(non.derm.IBM == '-')

# quantile regression
mod <- rq(Age ~ cancer.detailed + Sex, tau = 0.50, subset = subset, data = data)
# p-value
p <- coef(summary.rq(mod, se = 'ker'))[2:6, 4]
names(p) <- c('thymic vs overall', 'lung vs overall', 'melanoma vs overall',                     
              'gastrointestinal vs overall', 'genitourinary vs overall')
print(p)

# Build data
p <- c(p, NA)
df <- with(data, data.frame(cancer = cancer.detailed[subset],
                            age = Age[subset]))
df <- na.omit(df)

# Main 5 cancer types
main_groups <- c("thymic", "lung", "melanoma", "gastrointestinal", "genitourinary")
df_main <- df[df$cancer %in% main_groups, ]
df_main$group <- factor(df_main$cancer, levels = main_groups)
df_all <- df
df_all$group <- "All"
df_plot <- rbind(df_main, df_all)
df_plot$group <- factor(df_plot$group, levels = c(main_groups, "All"))

# Colors (keep your palette)
cols <- c('magenta', 'deepskyblue', 'indianred1',
          'royalblue', 'slateblue', 'darkslategray3')
names(cols) <- levels(df_plot$group)

top_y <- max(df$age, na.rm = TRUE) + 15
pval_df <- data.frame(group = levels(df_plot$group),
                      label = paste0("p = ", format(p, digits = 2, scientific = FALSE)),
                      y = top_y)
x_labs <- c(thymic = "TET", lung = "Lung cancer",
            melanoma = "Melanoma", 
            gastrointestinal = "Gastrointestinal cancer",
            genitourinary = "Genitourinary cancer",
            All = "All")

# Plot
p_age <- ggplot(df_plot, aes(x = group, y = age)) +
  geom_violin(
    aes(color = group),
    fill = "white",
    trim = FALSE,
    draw_quantiles = c(0.25, 0.5, 0.75),
    linewidth = 1.2,
    show.legend = FALSE
  ) +
  geom_jitter(
    aes(color = group),
    position = position_jitter(width = 0.15),
    size = 1.2,
    alpha = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cols) +
  geom_text(data = pval_df, aes(x = group, y = y, label = label), vjust = 0) +
  scale_x_discrete(labels = x_labs) +
  labs(x = "", y = "Age (year)", title = "Age") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, size = 1.2),
    axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1,
                                   margin = margin(t = 6)),
    axis.text.y     = element_text(size = 14),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    plot.title      = element_text(hjust = 0.5, size = 18),
    plot.margin     = margin(t = 6, r = 8, b = 16, l = 8))
print(p_age)

#################### forest plot for ICI therapy
cancer.type2 <- as.character(data$cancer.detailed)
cancer.type2[which(cancer.type2 == 'thymic')] <- NA
ll <- c()
for (i in 1:nrow(data)) {
  if (is.na(cancer.type2[i]) == TRUE) {
    ll[i] <- NA
  } else if (cancer.type2[i] == 'thymic') {
    ll[i] <- NA
  } else if (cancer.type2[i] == 'lung') {
    ll[i] <- 1
  } else if (cancer.type2[i] == 'melanoma') {
    ll[i] <- 2
  } else if (cancer.type2[i] == 'gastrointestinal') {
    ll[i] <- 3
  } else if (cancer.type2[i] == 'genitourinary') {
    ll[i] <- 4
  } else if (cancer.type2[i] == 'others') {
    ll[i] <- 5
  }
}
cancer.type2 <- reorder(cancer.type2, ll)
contrasts(cancer.type2) <- contr.sum(5)

# logistic regression
treat.ind <- ifelse(data$treatment == 'PD-1 & CTLA-4', 0, 1)
mod <- glm(treat.ind ~ cancer.type2 + Sex + Age, 
           family = 'binomial', data = data, subset = subset)

# p-value odds ratio
p <- coef(summary(mod))[2:5, 4]
odds <- coef(summary(mod))[2:5, 1]
ses <- coef(summary(mod))[2:5, 2]
odds <- cbind(odds, odds - qnorm(.975) * ses, odds + qnorm(.975) * ses)
odds <- exp(odds)

# forest plot
odds <- rbind(c(NA, NA, NA), odds)
mean <- odds[, 1]
lower <- odds[, 2]
upper <- odds[, 3]
pvs <- p
pvs <- round(pvs, 4)
pvs <- c(NA, pvs)

# groupwise proportion
pp <- tapply(treat.ind[subset], cancer.detailed[subset], function(x) mean(x, na.rm = TRUE))[1:5]
nn <- tapply(treat.ind[subset], cancer.detailed[subset], function(x) sum(x, na.rm = TRUE))[1:5]
NN <- tapply(treat.ind[subset], cancer.detailed[subset], function(x) length(na.omit(x)))[1:5]
prop <- paste0(nn, '/', NN, ', ', round(pp * 100))


base_data <- tibble::tibble(mean  = mean,
                            lower = lower,
                            upper = upper,
                            OR = round(mean, 2),
                            pvs = pvs,
                            prop = prop, 
                            variable = c('TET (100%) *', 
                                         'Lung cancer (89%) *', 
                                         'Melanoma (56%) *', 
                                         'Gastrointestinal cancer (91%) *',
                                         'Genitourinary cancer (76%) *'))

base_data |>
  forestplot(labeltext = c(variable, prop, OR, pvs),
             align = c("l", "c", "c", "c"),     # <-- Align columns
             xlim = c(0, 12),
             title = expression(bold("ICI therapy (overall prevalence of ICI monotherapy: 399/491, 81\\%)")),
             xlab = expression(bold("Association with ICI Monotherapy versus ICI combination therapy"))
             ,
             xlog = TRUE,
             boxsize = .35,
             lwd.ci = 5,                       # Thicker CI lines
             xticks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
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
  fp_add_header(variable = c("", "Cancer type\n"),
                OR = c("", "Adj-OR\n"),
                pvs = c("", "p-value\n"),
                prop = c("", "Monotherapy\n(n/N, %)")) |>
  fp_set_zebra_style("#EFEFEF")