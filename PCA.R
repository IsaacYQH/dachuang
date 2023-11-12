# Perform PCA using the correlation matrix
pca_cor <- prcomp(apply(as.matrix(data_raw[,-2:-1]),FUN = function(x){diff(log(x))},MARGIN = 2))
summary(pca_cor)

# Get the eigenvalues and eigenvectors for each PCA
eig_cor <- pca_cor$sdev^2
vect_cor <- pca_cor$rotation

# Draw the scree plot
plot(eig_cor, type = 'b', main = "Scree Plot (Correlation Matrix)")

# Compute the PC scores for each PCA
scores_cor <- pca_cor$x

# Plot all pairwise PC scores in a matrix plot
pairs(scores_cor, main = "Matrix Plot (Correlation Matrix)", col = "red")

X <- cbind(intercept=rep(1,nrow(data_raw)-1), scores_cor[,1:2])
