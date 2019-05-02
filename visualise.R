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

head(mtcars)

ggplot(mtcars, aes(x=wt, y=mpg)) +
  geom_point()

ggsave(filename = "~/Desktop/plot.eps", width = 10, height = 10)

### FUNCTIONS ###
# Used in GFS function. Bins score with range [0,1] into intervals
# E.g. 4 Intervals: Binned into 0.2, 0.4, 0.6, 0.8
bin <- function(score, num_intervals) {
  for (i in 1:num_intervals) {
    if (score <= i/num_intervals) {
      return (i/(num_intervals+1))
    }
  }
}

# Gene Fuzzy Scoring function transforms gene expression values
# Wilson Goh's paper
# Dense rank is used
GFS <- function(A, upper=0.05, lower=0.15, num_intervals=0) {
  print(sprintf("Top %.2f of expressed genes are assigned GFS scores of 1", upper))
  print(sprintf("Genes below the top %.2f of expressed genes are assigned GFS scores of 0", lower))
  # Rank function ranks largest value as 1 [-A is used]
  # Handle NaN?
  ranked_A <- apply(-A, 2, dense_rank)
  rownames(ranked_A) <- rownames(A)
  # Returns [1,] = upper, [2,] = lower
  qtile <- apply(ranked_A, 2, quantile, probs=c(upper, lower), names=F)
  
  if (num_intervals <= 0) {
    for (c in 1:ncol(ranked_A)) {
      # Calculate qtile range
      q_range <- qtile[2,c] - qtile[1,c]
      for (r in 1:nrow(ranked_A)) {
        if (ranked_A[r,c] <= qtile[1,c]) {
          # Assign 1s
          ranked_A[r,c] <- 1
        } else if (ranked_A[r,c] > qtile[2,c]){
          # Assign 0s
          ranked_A[r,c] <- 0
        } else {
          # Assign score
          score <- (qtile[2,c] - ranked_A[r,c]) / q_range
          ranked_A[r,c] <- score
        }
      }
    }
  } else {
    # Discrete intervals
    for (c in 1:ncol(ranked_A)) {
      # Calculate qtile range
      q_range <- qtile[2,c] - qtile[1,c]
      for (r in 1:nrow(ranked_A)) {
        if (ranked_A[r,c] <= qtile[1,c]) {
          # Assign 1s
          ranked_A[r,c] <- 1
        } else if (ranked_A[r,c] > qtile[2,c]){
          # Assign 0s
          ranked_A[r,c] <- 0
        } else {
          # Assign score
          score <- (qtile[2,c] - ranked_A[r,c]) / q_range
          # Round off score
          ranked_A[r,c] <- bin(score, num_intervals)
        }
      }
    }
  }
  return (ranked_A)
}

# Quantile normalisation
# Takes in array where columns are samples and rows are genes
qnorm <- function(arr) {
  sort_arr <- apply(arr,2,sort)
  ref_distr <- apply(sort_arr,1,mean)
  rank_arr <- apply(arr,2,dense_rank)
  qnorm_arr <- apply(rank_arr,c(1,2), function(x) ref_distr[x])
  return(qnorm_arr)
}

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

# Breast cancer data set
# control <- read.table('data/yeoh_2002/processed/processed_normal.tsv', header = T, row.names = 1)
# patient <- read.table('data/yeoh_2002/processed/processed_TEL-AML1.tsv', header = T, row.names = 1)

colnames(control_ds2) <- paste0("N", 1:10)
mcon_ds2 <- melt(control_ds2, variable.name = "ID")

# Plot pdf to see whether it has been quantile normalised
lcon_ds1 <- log2(control_ds1)
colnames(lcon_ds1) <- paste0("N", 1:10)
mtlgcon_ds1 <- melt(lcon_ds1, variable.name = "ID")

colnames(control_ds1) <- paste0("N", 1:10)
mtcon_ds1 <- melt(control_ds1, variable.name = "ID")
head(mtcon_ds1)

ggplot(mtcon_ds1, aes(x=value, color=ID)) + 
  geom_density()

mean(exp_tumour_ds2[90,])
hist(data.matrix(exp_tumour_ds2[90,]), breaks=30)

# Calculate fold change
control_mean <- apply(control_ds1, 1, mean)
tumour_mean <- apply(tumour_ds1, 1, mean)
fc <- tumour_mean/control_mean

# Visulise heatmap
arr1 <- data.matrix(control_ds1[fc > 2 | fc < 0.5,])
heatmap(arr1)

# Visualise rank stability within controls
raw_sd <- apply(control_ds1, 1, sd)
raw_mean <- apply(control_ds1, 1, mean)
raw_cv <- raw_sd/raw_mean

# Mean coeff var plot
plot(log2(raw_mean), log2(raw_cv))

# Rank vlaues
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
df <- cbind.data.frame(c(raw_cv, qnraw_cv, rank_cv), rep(c("Raw values", "Quantile normalised", "Rank values"), each=11774))
colnames(df) <- c("value","ID")

reorder(df$ID, rep(c(1,2,3), each=11774), order = T)

# Boxplot
boxplot(cbind(raw_cv,rank_cv, qnraw_cv), outline=F)

qnctrl1 <- melt(qnctrl_ds1)[,2:3]
lgqnctrl1 <- melt(log2(qnctrl_ds1))[,2:3]

ggplot(lgqnctrl1, aes(x=value, color=Var2)) + 
  geom_density()

# Violin plot
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

# Manipulate for base plot
par(mfrow=c(1,1))
boxplot(cpare)

control_var <- apply(control_ds1, 1, var)
control_var <- apply(control_ds1, 1, var)

# Visualise whether low expression genes have high coeff var
plot(raw_mean, raw_cv)

n <- rnorm(1000,0,1)
p <- rnorm(1000,1,1)

p1 <- sapply(n, function(x) x-p[1])
p2 <- sapply(n, function(x) x-p[2])
p3 <- sapply(n, function(x) x-p[3])
p4 <- sapply(n, function(x) x-p[4])
p5 <- sapply(n, function(x) x-p[5])
p6 <- sapply(n, function(x) x-p[6])
p7 <- sapply(n, function(x) x-p[7])
p8 <- sapply(n, function(x) x-p[8])
p9 <- sapply(n, function(x) x-p[9])
p10 <- sapply(n, function(x) x-p[10])

plot(density(c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)))
par(new=T)
curve(dnorm(x,-0.65,0.5), xlim=c(-3,3))


# PCA ---------------------------------------------------------------------
pc_cars <- prcomp(mtcars[,3:6], center = T, scale. = T)
pc_cars_p <- prcomp(mtcars[,3:6], center = T, scale. = F)
pc_cars$x[,1]

par(mfrow = c(1,2))
plot(pc_cars$x[,1], pc_cars$x[,2])
plot(pc_cars_p$x[,1], pc_cars_p$x[,2])

# 3D ----------------------------------------------------------------------
x <- 1:20
y <- 1:20

f <- function(x,y) {
  return (2*x + 3*y)
}

z <- outer(x,y,f)

# 3D surface plot using persp
plt <- persp(x,y,z, ticktype = "detailed", xlim=c(-5,5))

# Adding points and lines to persp
trans3d <- function(x,y,z, pmat) {
  tr <- cbind(x,y,z,1) %*% pmat
  list(x = tr[,1]/tr[,4], y= tr[,2]/tr[,4])
}

coord <- trans3d(1,1,f(1,1), plt)
points(trans3d(1,1,f(1,1), plt))

z2 <- sapply(1:length(x),function(n)f(x[n],y[n]))
lines(trans3d(x,y,z2,plt),col="red",lwd=2)
lines(trans3d(c(-10,10,10,-10,-10),c(-10,-10,10,10,-10),c(2,2,8,8,2), pmat), col="blue")

# Plot 3D scatter
library(rgl)
rgl.open()
rgl.bg(color="white")
plot3d(x,y,z)


# CLUSTERING -------------------------------------------------
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