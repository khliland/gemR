library(gemR)
data(MS)
ms <- as.character(MS$MS)
cluster <- as.numeric(MS$cluster)
proteins <- unclass(MS$proteins)
##################################

# The example assumes that the following objects are available in
# R's Global Environment:
# proteins - numeric matrix of proteins
# ms       - character vector of MS ("no"/"yes")
# cluster  - numeric vector of cluster origin

# Organise data as data.frame.
MS.data <- data.frame(proteins = I(proteins),
                     ms = factor(ms),
                     cluster = factor(cluster))

# Step 1 - GLM:
MS.gem <- GEM(proteins ~ ms + cluster, data = MS.data)

# Matrices of effects and ER values can be extracted for custom analyses
E.ms       <- MS.gem$effects$ms
E.cluster  <- MS.gem$effects$cluster
ER.ms      <- MS.gem$ER.values$ms
ER.cluster <- MS.gem$ER.values$cluster

# Step 2 - PLS
ncomp <- 20
# PLS analysis of 'ms' effect
MS.pls.ms <- pls(MS.gem, 'ms', ncomp, validation = "CV",
                 segments = 10, segment.type = "interleaved",
                 jackknife = TRUE, df.used = 1)

# Element extraction
# Cross-validated classifications of 'ms' (# correct)
colSums(MS.pls.ms$classes == as.numeric(MS.data$ms))
# Parsimonious estimate from the above
MS.ncomp       <- 4
# Coefficients for 'ms' effect
MS.coef.ms     <- coef(MS.pls.ms)
# Scores and loadings
MS.scores.ms   <- scores(MS.pls.ms)
MS.loadings.ms <- loadings(MS.pls.ms)
# Jackknife P-value estimates and significant subset
MS.jack.ms     <- MS.pls.ms$jack[,1,4]
MS.signif.ms   <- which(MS.jack.ms < 0.05)

# Example plots
oldpar <- par(mfrow = c(2,2), mar = c(4,4,2,1))
scoreplot(MS.pls.ms, main = "Scores", panel.first = abline(h=0, v=0, col="gray"),
          col = MS.data$ms, pch = as.numeric(MS.data$cluster))
loadingplot(MS.pls.ms, labels = "names", scatter = TRUE, main = "Loadings",
            panel.first = abline(h=0, v=0, col="gray"))
corrplot(MS.pls.ms, main = "Correlation loadings", label = "names")
plot(colSums(MS.pls.ms$classes == as.numeric(MS.data$ms)),
     ylab = "# correct", xlab = "# components", panel.first = grid(),
     main = "Classification accuracy")
par(oldpar)

# Step 2 - Elastic Net
# Elastic Net analysis of the 'cluster' effect
MS.en.cluster <- elastic(MS.gem, 'cluster', validation = "LOO",
                         alpha = 0.5, family = "binomial")

# Element extraction
# Coefficients for 'cluster' effect
MS.coef.cluster <- coef(MS.en.cluster)
# Non-zero coefficients at optimal complexity
nonZeroID <- which(MS.coef.cluster[,1] != 0)
nonZeroLab <- names(nonZeroID)

# Example plots
oldpar <- par(mfrow = c(3,1), mar = c(4,4,2,1))
plot(MS.en.cluster$glmnet, xvar = "lambda")
plot(MS.en.cluster)
par(oldpar)

