library(ggplot2)
library(dplyr)
library(reshape2)
library(gplots)
library(RColorBrewer)
source("functions.R")

setwd("~/projects/phd/diff_expr/")

# SETTINGS ----------------------------------------------------------------
# ggplot2 settings
theme_set(theme_bw())
theme_update(text = element_text(size=48),
             plot.margin = unit(c(1,1,1,1), "cm"),
             panel.border = element_rect(size = 2),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

# Color palette
custom_colour <- colorRampPalette(brewer.pal(11, "RdYlGn"))(50)

# Ovarian data set 1
# Not log2
ds1 <- read.table('data/ovarian_cancer/GSE18521/processed/processed_original.tsv', header = T, row.names = 1)
control_ds1 <- ds1[54:63]
tumour_ds1 <- ds1[1:53]

# Ovarian data set 2 was pre-processed using RMA
# Log2 and quantile normalised
ds2 <- read.table('data/ovarian_cancer/GSE26712/processed/processed_original.tsv', header = T, row.names = 1)
control_ds2 <- ds2[1:10]
tumour_ds2 <- ds2[11:195]

exp_control_ds2 <- 2^(control_ds2)
exp_tumour_ds2 <- 2^(tumour_ds2)

# Plot pdf to see whether it has been quantile normalised
lcon_ds1 <- log2(control_ds1)

colnames(lcon_ds1) <- paste0("N", 1:10)
mtlgcon_ds1 <- melt(lcon_ds1, variable.name = "ID")
colnames(control_ds1) <- paste0("N", 1:10)
mtcon_ds1 <- melt(control_ds1, variable.name = "ID")

ggplot(mtcon_ds1, aes(x=value, color=ID)) + 
  geom_density()

mean(exp_tumour_ds2[90,])
hist(data.matrix(exp_tumour_ds2[90,]), breaks=30)

# Calculate fold change
control_mean <- apply(control_ds1, 1, mean)
tumour_mean <- apply(tumour_ds1, 1, mean)
fc <- tumour_mean/control_mean

# Visualise rank stability within controls
raw_sd <- apply(control_ds1, 1, sd)
raw_mean <- apply(control_ds1, 1, mean)
raw_cv <- raw_sd/raw_mean

# Mean coeff var plot
plot(log2(raw_mean), log2(raw_cv))

# Rank values
control_rank <- apply(control_ds1, 2, rank)
rank_sd <- apply(control_rank, 1, sd)
rank_mean <- apply(control_rank, 1, mean)
rank_cv <- rank_sd/rank_mean

# Quantile norm
ctrl_ds1 <- data.matrix(control_ds1)
qnctrl_ds1 <- qnorm(ctrl_ds1)
qnraw_sd <- apply(qnctrl_ds1, 1, sd)
qnraw_mean <- apply(qnctrl_ds1, 1, mean)
qnraw_cv <- qnraw_sd/qnraw_mean

# Manipulate for ggplot
df <- cbind.data.frame(c(raw_cv, qnraw_cv, rank_cv),
                       rep(c("Raw values", "Quantile normalised", "Rank values"), each=11774))
colnames(df) <- c("value","ID")

reorder(df$ID, rep(c(1,2,3), each=11774), order = T)


# Density plot ------------------------------------------------------------
qnctrl1 <- melt(qnctrl_ds1)[,2:3]
lgqnctrl1 <- melt(log2(qnctrl_ds1))[,2:3]

ggplot(lgqnctrl1, aes(x=value, color=Var2)) + 
  geom_density()

# Boxplot -----------------------------------------------------------------
boxplot(cbind(raw_cv,rank_cv, qnraw_cv), outline=F)

# Violin plot -------------------------------------------------------------
ggplot(df, aes(reorder(df$ID, rep(c(1,2,3), each=11774), order = T), y=value)) + 
  geom_violin(aes(fill=reorder(df$ID, rep(c(1,2,3), each=11774), order = T))) +
  labs(y = "Coefficient of variation") +
  theme(legend.position="none",
        axis.title.x=element_blank())
# axis.text.x=element_blank(),
# axis.ticks.x=element_blank())

ggsave("stability.eps", width = 20, height = 20, units = "cm")

ggplot(df, aes(reorder(df$ID, rep(c(1,2,3), each=11774), order = T), y=value)) + 
  geom_violin(aes(fill=reorder(df$ID, rep(c(1,2,3), each=11774), order = T))) +
  labs(y = "Coefficient of variation") +
  theme(legend.position="none",
        axis.title.x=element_blank())
# axis.text.x=element_blank(),
# axis.ticks.x=element_blank())

ggplot(df, aes(ID, value)) + 
  geom_violin(aes(fill=ID))

ggplot(mtcars, aes(factor(cyl), mpg)) + 
  geom_violin(aes(fill=factor(cyl))) +
  labs(y = "Length") +
  theme(legend.position="none",
        axis.title.x=element_blank())

# 2D Scatterplot ----------------------------------------------------------
# R base plot
plot(sum_chips, col = rep(1:6, each = 5), pch = 19)
# text(sum_chips, labels,
#      cex = 0.6, srt = 0, adj = c(0,-10))

# ggplot2
head(mtcars)

ggplot(mtcars, aes(x=wt, y=mpg)) +
  geom_point()

ggsave(filename = "~/Desktop/plot.eps", width = 10, height = 10)
# 3D Scatterplot ----------------------------------------------------------
library(rgl)
rgl.open()
rgl.bg(color="white")
plot3d(x,y,z)

persp3d(volcano)
rgl.postscript("dump/volcano.eps")

f <- function(x, y) {
  z = ((x^2)+(3*y^2)) * exp(-(x^2)-(y^2))
}

# Plot a 3D function surface plot
plot3d(f,
       xlim = c(-3, 3), ylim = c(-3, 3),)

# Plot a 3D density plot
# Use MASS package to calculate densities
plot_3d <- MASS::kde2d(A,M)
persp3d(plot_3d)

# PCA ---------------------------------------------------------------------
pc_cars <- prcomp(mtcars[,3:6], center = T, scale. = T)
pc_cars_p <- prcomp(mtcars[,3:6], center = T, scale. = F)
pc_cars$x[,1]

par(mfrow = c(1,2))
plot(pc_cars$x[,1], pc_cars$x[,2])
plot(pc_cars_p$x[,1], pc_cars_p$x[,2])

# Clustering -------------------------------------------------
# Import ovarian cancer data set 1
# Not log2
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/GSE18521_entrez.tsv",
                            header = T, row.names = 1)
control_data1 <- ovarian_data1[1:10]
tumour_data1 <- ovarian_data1[11:63]

# Ovarian cancer data set 2
# Pre-processed using RMA - Log2 and quantile normalised
ovarian_data2 <- read.table('data/ovarian_cancer/GSE26712/processed/GSE27612_entrez.tsv',
                            header = T, row.names = 1)
# Reverse the log2
ovarian_data2a <- 2^ovarian_data2
control_data2 <- ovarian_data2a[1:10]
tumour_data2 <- ovarian_data2a[11:195]

# Significant genes for data set 1
pvalue_1 <- row_ttest(control_data1, tumour_data1)
# Assign row indices
names(pvalue_1) <- 1:length(pvalue_1)
row_index1 <- as.numeric(names(sort(pvalue_1)[1:20]))
head(sort(pvalue_1))
# Remove integer row names as they affect the indexing
fltr_data1 <- ovarian_data1[row_index1,]
# Save as EPS file
setEPS()
postscript("heatmap1.eps", width = 10, height = 6)
# Set margins around figure
# par(mar = rep(4, 4))
heatmap.2(data.matrix(fltr_data1),
          cellnote = round(data.matrix(fltr_data1), 1),
          notecol = "black",
          notecex = 0.2,
          trace = "none",
          key = F,
          col = custom_colour,
          lwid = c(0.5,4),
          lhei = c(0.5,4),
          # labCol = c(paste0("N",1:10), paste0("P",1:53)),
          scale = "row")
dev.off()

# Significant genes for data set 2
pvalue_2 <- row_ttest(control_data2, tumour_data2)
# Assign row indices
names(pvalue_2) <- 1:length(pvalue_2)
row_index2 <- as.numeric(names(sort(pvalue_2)[1:20]))
# Remove integer row names as they affect the indexing
fltr_data2 <- ovarian_data2[row_index2,]
# Save as EPS file
setEPS()
postscript("heatmap2.eps", width = 12, height = 6)
# Set margins around figure
# par(mar = rep(4, 4))
heatmap.2(data.matrix(fltr_data2),
          trace = "none",
          key = F,
          col = custom_colour,
          cexCol = 0.4,
          lwid = c(0.5,4),
          lhei = c(0.8,4),
          labCol = c(paste0("N",1:10), paste0("P",1:185)),
          scale = "row")
dev.off()

# AFTER REMOVING
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/GSE18521_entrez1.tsv",
                            header = T, row.names = 1)
control_data1 <- ovarian_data1[1:10]
tumour_data1 <- ovarian_data1[11:62]
# Significant genes for data set 1
pvalue_1 <- row_ttest(control_data1, tumour_data1)
# Assign row indices
names(pvalue_1) <- 1:length(pvalue_1)
row_index1 <- as.numeric(names(sort(pvalue_1)[1:20]))
head(sort(pvalue_1))
# Remove integer row names as they affect the indexing
fltr_data1 <- ovarian_data1[row_index1,]
# Save as EPS file
setEPS()
postscript("heatmap1.eps", width = 10, height = 6)
# Set margins around figure
# par(mar = rep(4, 4))
heatmap.2(data.matrix(fltr_data1),
          cellnote = round(data.matrix(fltr_data1), 1),
          notecol = "black",
          notecex = 0.2,
          trace = "none",
          key = F,
          col = custom_colour,
          lwid = c(0.5,4),
          lhei = c(0.5,4),
          labCol = c(paste0("N",1:10), paste0("P",1:52)),
          scale = "row")
dev.off()

# Multiple plots ----------------------------------------------------------
# Mix of base plots and ggplots
# Settings for base plot
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    mar = c(6,2,2,2))
# Total probe intensities for each chip
sum_chips <- apply(df, 2, sum)
plot(sum_chips,
     col = rep(1:6, each = 5), pch = 19)
# text(sum_chips, labels,
#      cex = 0.6, srt = 0, adj = c(0,-10))
sum_plot <- recordPlot()
# Boxplots of chips
boxplot(df, las = 2)
boxplot <- recordPlot()
# Plot density curves
pdf <- plot_pdf(df)
# Plot multiple plots
multiplot <- plot_grid(boxplot, sum_plot, pdf,
                       nrow = 3)
save_plot("dump/before.eps", before_norm,
          base_height = 10, base_width = 8)

# SAVE --------------------------------------------------------------------
save_eps(boxplot(df), "dump/test.eps")

# Plot PDF ----------------------------------------------------------------
density_nospike <- plot_density(qnorm_nospike) + xlim(0, 7)
  geom_vline(xintercept = max_classA) +
  geom_vline(xintercept = max_classB) +
  stat_function(fun = dgamma, n = 500,
                args = list(shape = gamma_shape, rate = gamma_rate),
                col = "blue") +
  xlim(0, 7)
