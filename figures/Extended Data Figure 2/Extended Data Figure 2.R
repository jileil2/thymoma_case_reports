##################### set working directory
setwd("/Users/jileilin/Desktop/Research/thymomal case reports/")

##################### load packages
library("strucchange")
library('nnet')
library('quantreg')
library('multcomp')
library('forestplot')

#################### read data
data <- read.csv('data1.csv')[, -1]

#################### Age at MG Onset 
df <- with(data, data.frame(Age = na.omit(Age)))

ggplot(df, aes(x = Age)) +
  geom_histogram(
    binwidth = 10,
    fill = "white",
    color = "black",
    size = 2
  ) +
  labs(
    x = expression(bold("Age (All, N = 191)")),
    y = "Frequency",
    title = expression(bold("Age at MG Onset"))
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 2, color = "black"),
    axis.text = element_text(color = "black", size = 17),       # Bigger tick labels
    axis.title.x = element_text(color = "black", face = "bold", size = 18),
    axis.title.y = element_text(color = "black", face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 3),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )

#################### Age at MG Onset (Female)
df_female <- with(data, data.frame(Age = na.omit(Age[Sex == "F"])))

ggplot(df_female, aes(x = Age)) +
  geom_histogram(
    binwidth = 10,
    fill = "white",
    color = "magenta",   # Match base R border color
    size = 2             # Border thickness
  ) +
  labs(
    x = expression(bold("Age (Female, N = 115)")),
    y = "Frequency",
    title = expression(bold("Age at MG Onset"))
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 2, color = "black"),
    axis.text = element_text(color = "black", size = 17),       # Larger tick labels
    axis.title.x = element_text(color = "black", face = "bold", size = 18),
    axis.title.y = element_text(color = "black", face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 3),  # Thick box around plot
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )

dev.off()

#################### Age at MG Onset (Male)
df_male <- with(data, data.frame(Age = na.omit(Age[Sex == "M"])))

ggplot(df_male, aes(x = Age)) +
  geom_histogram(
    binwidth = 10,
    fill = "white",
    color = "deepskyblue",   # Match base R border color
    size = 2                 # Border thickness
  ) +
  labs(
    x = expression(bold("Age (Male, N = 76)")),
    y = "Frequency",
    title = expression(bold("Age at MG Onset"))
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 2, color = "black"),
    axis.text = element_text(color = "black", size = 17),
    axis.title.x = element_text(color = "black", face = "bold", size = 18),
    axis.title.y = element_text(color = "black", face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 3),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )

#################### Age at Myositis Onset
df <- with(data, data.frame(Age = na.omit(Age)))

ggplot(df, aes(x = Age)) +
  geom_histogram(
    binwidth = 10,
    fill = "white",
    color = "black",
    size = 2
  ) +
  labs(
    x = expression(bold("Age (All, N = 191)")),
    y = "Frequency",
    title = "Age at Myositis Onset"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 2, color = "black"),
    axis.text = element_text(color = "black", size = 17),
    axis.title.x = element_text(color = "black", face = "bold", size = 18),
    axis.title.y = element_text(color = "black", face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 3),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )

dev.off()

##################### Age at Myositis Onset (Female)
df_female <- with(data, data.frame(Age = na.omit(Age[Sex == "F"])))

ggplot(df_female, aes(x = Age)) +
  geom_histogram(
    binwidth = 10,
    fill = "white",
    color = "magenta",
    size = 2
  ) +
  labs(
    x = expression(bold("Age (Female, N = 115)")),
    y = "Frequency",
    title = "Age at Myositis Onset"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 2, color = "black"),
    axis.text = element_text(color = "black", size = 17),
    axis.title.x = element_text(color = "black", face = "bold", size = 18),
    axis.title.y = element_text(color = "black", face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 3),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )
dev.off()


##################### Age at Myositis Onset (Male)
df_male <- with(data, data.frame(Age = na.omit(Age[Sex == "M"])))

ggplot(df_male, aes(x = Age)) +
  geom_histogram(
    binwidth = 10,
    fill = "white",
    color = "deepskyblue",
    size = 2
  ) +
  labs(
    x = expression(bold("Age (Male, N = 76)")),
    y = "Frequency",
    title = "Age at Myositis Onset"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 2, color = "black"),
    axis.text = element_text(color = "black", size = 17),
    axis.title.x = element_text(color = "black", face = "bold", size = 18),
    axis.title.y = element_text(color = "black", face = "bold", size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 3),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )
