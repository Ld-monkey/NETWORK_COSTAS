# Install all packages.
#BiocManager::install(c("minet", "Rgraphviz"))
#install.packages(c("corrplot", "ggplot2"), type="source", repos="https://cloud.r-project.org/")

# make sure you have the packages below before you try to import them.
options(jupyter.plot_mimetypes = c("text/plain", "image/png" ))
library("minet")
library("Rgraphviz")
library("ggplot2")
library("corrplot")
library("gplots")
library("RColorBrewer")
# Simulated gene expression data
data(syn.data)
dim(syn.data)

# Simulated gene expression data
data(syn.data)
dim(syn.data)

# The columns are the names of the genes.
head(syn.data)

# These are the names of the conditions
rownames(syn.data)

# Also source a file with some small functions that I developed to help the, → transition to cytoscape.
source("rFunctions_TPnote.r")

# --------------------- Part 1 ---------------------

options(repr.plot.width = 6, repr.plot.height = 7)
# Explore the whole dataset for Pearson correlation
heatmap.2(t(as.matrix(syn.data)), # We transpose the dataset to have genes as → rows and conditions as columns.
          trace="none",
          distfun=function(x) as.dist(1-abs(cor(t(x)))),
          col=rev(heat.colors(32)),
          main="Heatmap of GE data with 
          Pearson method", key.title="")

### Get the correlation matrix, cluster , viualise and infere the network
# Pearson corelation
corM <- cor(syn.data) # == corM<- cor(syn.data, method = "pearson")

# We will visualise the correlation matrix with the corrplot function.
options(repr.plot.width = 5, repr.plot.height = 6)
corrplot(corM, type="upper", diag=FALSE,
         order="hclust", hclust.method="ward.D",
         cl.pos="n", tl.cex=0.5,
         main="Pearson correlation plot of GE data.", mar=c(0,0,1,0), 
         cex.main=0.75)
# the function applies a hierarchical clustering of the matrix of the pearson → correlation coefficient.  

# Do the same thing for (Spearman and Kendall rank correlations are both available in R).
# With euclidian distance dist.
# and by and by trying to rise the expression matrix to a power (i.e. like the “soft
# thresholding” approach we show in the WGCRN TP)

# Spearman Headmap
heatmap.2(t(as.matrix(syn.data)), # We transpose the dataset to have genes as → rows and conditions as columns.
          trace="none",
          distfun=function(x) as.dist(cor(t(x), method = "spearman")),
          col=rev(heat.colors(32)),
          main="Heatmap of GE data 
          with Spearman method", key.title="")

# Spearman correlation 
corSpearman <- cor(syn.data, method = "spearman")

# Plot Spearman correantion of GE data.
corrplot(corSpearman, type="upper", diag=FALSE,
         order="hclust", hclust.method="ward.D",
         cl.pos="n", tl.cex=0.5,
         main="Spearman correlation plot of GE data.", mar=c(0,0,1,0), 
         cex.main=0.75)

# Kendall Headmap
heatmap.2(t(as.matrix(syn.data)), # We transpose the dataset to have genes as → rows and conditions as columns.
          trace="none",
          distfun=function(x) as.dist(cor(t(x), method = "kendall")),
          col=rev(heat.colors(32)),
          main="Heatmap of GE data 
          with Kendall method", key.title="")

# kendall correlation 
corKendal <- cor(syn.data, method = "kendall")

# Plot kendal correantion of GE data.
corrplot(corKendal, type="upper", diag=FALSE,
         order="hclust", hclust.method="ward.D",
         cl.pos="n", tl.cex=0.5,
         main="Kendall correlation plot of GE data.", mar=c(0,0,1,0), 
         cex.main=0.75)

# by rising the correlation matrix to an exponent
corExponant <- cor(syn.data)^2

# 1.2.1 Thresholding the correlation matrix – Infer the network.

# A function to apply a threshold value on a matrix m
matrix_threshold <- function(m, t, both = TRUE, ...) {
  diag(m) <- 0;
  m[m >= -t & m <= t] <-0;
  return(m)
}

# Apply a threshold t = 0.75.
corMT <- matrix_threshold(corM, t=0.75)

options(repr.plot.width = 5, repr.plot.height = 6)

# plot the network
corrplot(corMT, tl.cex=0.5, 
         main="Network 'adjacency' matrix", 
         mar=c(0,0,1,0),
         cex.main=0.75)

# 1.3 - Apply 3 differents thresholds [0.6, 0.8, 0.99]
corMT_one <- matrix_threshold(corM, t=0.6) # 0.6
corMT_two <- matrix_threshold(corM, t=0.8) # 0.8
corMT_third <- matrix_threshold(corM, t=0.99) # 0.99

# plot all networks.
# 0.6
corrplot(corMT_one, tl.cex=0.5, 
         main="Network 'adjacency' matrix", 
         mar=c(0,0,1,0),
         cex.main=0.75)
# 0.8
corrplot(corMT_two, tl.cex=0.5, 
         main="Network 'adjacency' matrix", 
         mar=c(0,0,1,0),
         cex.main=0.75)
# 0.99
corrplot(corMT_third, tl.cex=0.5, 
         main="Network 'adjacency' matrix", 
         mar=c(0,0,1,0),
         cex.main=0.75)
