library(edgeR)

# sample names and replicate count
tissue1 <- "ectodermalcell"
tissue2 <- "H1"
tissue1_count <- 6
tissue2_count <- 2

# read in raw count data
setwd("~/Differential-Gene-Expression-Analysis/limma/")
counts <- read.delim("counts", skip=1, header=TRUE, sep='\t', row.names = 1)

# exclude extra information
counts = subset(counts, select=-c( Chr, Start, End, Strand, Length))
names(counts)

# create DGEList object
dge <- DGEList((counts))

# Get normalization factors
dge <- calcNormFactors(dge)
head(dge)

# drop low-expressed genes (zero count in all samples)
cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
dge <- dge[-drop,] 
dim(d) # number of genes left

# specify model to be fitted
sample <- factor(rep(c(tissue1, tissue2), times = c(tissue1_count, tissue2_count)))
design.mat <- model.matrix(~0+sample)
colnames(design.mat) <- levels(sample)
dim(design.mat) # total num. replicates * num samples 

# Counts are transformed to logCPM based on the calculated normalization factors
y <- voom(dge, design.mat, plot=T)

# weighted least squares for each gene
fit = lmFit(y, design.mat)
head(coef(fit)) # coefficients of gene

# specify samples to compare
contrast.mat <- makeContrasts(
  Diff = get(tissue1) - get(tissue2), 
  levels = design.mat
)

# to obtain log fold-changes as contrasts
fit = contrasts.fit(fit, contrast.mat)

# shrink standard errors 
fit = eBayes(fit)

# write all logFC scores
write.table(coef(fit),"logFC.csv",sep="\t", col.names = FALSE)

# most differentially expressed genes
top.genes <- topTable(fit, sort.by = "P", n = Inf)
head(top.genes)
length(which(top.genes$adj.P.Val < 0.05))

