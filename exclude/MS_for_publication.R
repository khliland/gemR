library(gemR)
data(MS)
ms <- as.character(MS$MS)
group <- as.numeric(MS$group)
proteins <- unclass(MS$proteins)
sex <- as.character(MS$sex)
age <- as.numeric(MS$age)
##################################

# The example assumes that the following objects are available in
# R's Global Environment:
# proteins - numeric matrix of proteins
# ms       - character vector of MS ("no"/"yes")
# group    - numeric vector of groups discovered by cluster analysis
# sex      - character vector of sex ("M"/"F")
# age      - numeric vector of patient ages

# Step 0 - Prepare data
# Organise data as data.frame.
MS.data <- data.frame(proteins = I(proteins),
                     ms = factor(ms),
                     group = factor(group),
                     sex = factor(sex),
                     age = age)

# Step 1 - GLM:
MS.gem <- GEM(proteins ~ ms + group + sex + age, data = MS.data)

# Matrices of effects and ER values can be extracted for custom analyses.
# Here shown for the 'ms' effect
E.ms       <- MS.gem$effects$ms
ER.ms      <- MS.gem$ER.values$ms

# Step 2 - PLS
ncomp <- 10
# PLS analysis of 'ms' effect
ms.pls <- pls(MS.gem, 'ms', ncomp,
                 jackknife = TRUE)

# Element extraction
# Cross-validated classifications of 'ms' (# correct)
colSums(ms.pls$classes == as.numeric(MS.data$ms))
# Parsimonious estimate of #components
ms.ncomp    <- 2
# Coefficients, scores and loadings for 'ms' effect
ms.coef     <- coef(ms.pls)
ms.scores   <- scores(ms.pls)
ms.loadings <- loadings(ms.pls)
# Jackknife P-value estimates, significant and subset
ms.jack     <- ms.pls$jack[,1,ms.ncomp]
ms.signif   <- ms.jack < 0.05
ms.subset   <- which(ms.signif)

# Example plots
oldpar <- par(mfrow = c(2,2), mar = c(4,4,2,1))
scoreplot(ms.pls, main = "Scores (ms)", panel.first = abline(h=0, v=0, col="gray"),
          col = MS.data$ms, pch = as.numeric(MS.data$group), cex = 0.7)
loadingplot(ms.pls, scatter = TRUE, main = "Loadings (ms)",
            col = c("gray", "black")[ms.signif+1], pch = c(1,15)[ms.signif+1],
            panel.first = abline(h=0, v=0, col="gray"), cex = 0.7)
corrplot(ms.pls, main = "Correlation loadings",
         col = c("gray", "black")[ms.signif+1], pch = c(1,15)[ms.signif+1],
         cex = 0.3)
plot(ms.pls, ylim = c(0,100), cex = 0.7)
par(oldpar)

# Step 2 - Elastic Net
# Elastic Net analysis of the 'group' effect
MS.en.group <- elastic(MS.gem, 'group', validation = "LOO",
                         alpha = 0.5, family = "binomial")

# Element extraction
# Coefficients for 'group' effect
MS.coef.group <- coef(MS.en.group)
# Non-zero coefficients at optimal complexity
nonZeroID <- which(MS.coef.group[-1,1] != 0)
nonZeroLab <- names(nonZeroID)

# Example plots
oldpar <- par(mfrow = c(3,1), mar = c(4,4,2,1))
plot(MS.en.group$glmnet, xvar = "lambda")
plot(MS.en.group)
par(oldpar)


#############################
# PLS repeated on group
group.pls <- pls(MS.gem, 'group', ncomp,
              jackknife = TRUE)

# Element extraction
# Cross-validated classifications of 'group' (# correct)
colSums(group.pls$classes == as.numeric(MS.data$group))
# Parsimonious estimate of #components
group.ncomp    <- 2
# Coefficients, scores and loadings for 'group' effect
group.coef     <- coef(group.pls)
group.scores   <- scores(group.pls)
group.loadings <- loadings(group.pls)
# Jackknife P-value estimates, significant and subset
group.jack     <- group.pls$jack[,1,group.ncomp]
group.signif   <- group.jack < 0.05
group.subset   <- which(group.signif)

# Example plots
oldpar <- par(mfrow = c(2,2), mar = c(4,4,2,1))
scoreplot(group.pls, main = "Scores (group)", panel.first = abline(h=0, v=0, col="gray"),
          col = MS.data$group, pch = as.numeric(MS.data$group), cex = 0.7)
loadingplot(group.pls, scatter = TRUE, main = "Loadings (group)",
            col = c("gray", "black")[group.signif+1], pch = c(1,15)[group.signif+1],
            panel.first = abline(h=0, v=0, col="gray"), cex = 0.7)
corrplot(group.pls, main = "Correlation loadings",
         col = c("gray", "black")[group.signif+1], pch = c(1,15)[group.signif+1],
         cex = 0.3)
plot(group.pls, ylim = c(0,100), cex = 0.7)
par(oldpar)
