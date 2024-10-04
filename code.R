#install.packages("BiocManager")
#BiocManager::install("GEOquery")
library("GEOquery")

gse <- getGEO("GSE49126")
gse <- gse[[1]]
show(gse)

head(exprs(gse))
ex <- exprs(gse)
dim(ex)
colnames(ex)

boxplot(ex) #asymmetric

ex1 <- log2(ex)

boxplot(ex1)

channel_medians <- apply(ex1, 2, median)
ex2 <- sweep(ex1, 2, channel_medians, "-")

#norm_log_x <- scale(norm_log_x)

boxplot(ex2)

#ex2 <- na.omit(as.matrix(ex2))

# PCA ----
pca <- prcomp(t(ex2)) #genes on the columns and samples on the rows
summary(pca)
screeplot(pca)

# draw PCA plot
grpcol <- c(rep("lightblue",20), rep("magenta",30))
plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", 
     main="PCA for components 1&2", type="p", pch=10, col=grpcol)
text(pca$x[,1], pca$x[,2], rownames(pca$x), cex=0.75)

#install.packages("rgl")
library("rgl")
scores = as.data.frame(pca$x)
plot3d(scores[,1:3], size=5, col = grpcol)
#text3d(scores[,1:3], texts=c(rownames(scores)), cex=0.7, pos=3)


# K-means ----
#BiocManager::install("useful")
library("useful")
library('ggplot2')

k <- 2
kmeans_result <- kmeans(t(ex2), k)
table(kmeans_result$cluster)
kmeans_data <- data.frame(sample = colnames(ex2), group = factor(kmeans_result$cluster))
kmeans_data$shape <- c(rep(15, 20), rep(16, 30))
# Perform PCA for dimensionality reduction for plotting
pca <- prcomp(t(ex2))

kmeans_data$PC1 <- pca$x[,1]
kmeans_data$PC2 <- pca$x[,2]

ggplot(kmeans_data, aes(x = PC1, y = PC2, color = group, shape = factor(shape))) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(15, 16), labels = c("controls", "cases")) +
  labs(title = "K-means Clustering", x = "Principal Component 1", y = "Principal Component 2") +
  guides(color = guide_legend(override.aes = list(shape = 15)), 
         shape = guide_legend(override.aes = list(color = "black")))

dist_matrix <- dist(t(ex2)) #try with different methods
hc_result <- hclust(dist_matrix, method = "single") #try with different methods
k <- 2
groups <- cutree(hc_result, k=k)
table(groups)
plot(hc_result, hang <- -1, labels=gse$title)
rect.hclust(hc_result, k = 2, which = NULL, x = NULL, h = NULL, border = 2, cluster = NULL) 


# Random Forest ----
#BiocManager::install("ALL")
library("ALL"); data(ALL)
# keep only 30 arrays JUST for computational convenience
e.mat <- (2^ex2)
small.eset <- log2(na.omit(e.mat))
dim(small.eset)
group <- c(rep('C',20),rep('P',30)) # classification, in order

# build RF
#BiocManager::install("randomForest") 
library(randomForest)
set.seed(1234)
rf <- randomForest(x=t(small.eset), y=as.factor(group), ntree=1000)
predict(rf, t(small.eset[, 1:5]))
plot(sort(rf$importance, decreasing=TRUE)) # can also use: varImpPlot(rf)
#extract the most 'important' genes
probe.names <- rownames(rf$importance)
top200 <- probe.names[order(rf$importance, decreasing=TRUE)[1:200]]
write.csv(top200, file = "probes-top200.txt", quote=FALSE, row.names = FALSE, col.names=FALSE)

# Look at variable importance
imp.temp <- abs(rf$importance[,])
t <- order(imp.temp,decreasing=TRUE)
plot(c(1:nrow(small.eset)),imp.temp[t],log='x',cex.main=1.5,
     xlab='gene rank',ylab='variable importance',cex.lab=1.5,
     pch=16,main='ALL subset results')
# Get subset of expression values for 25 most 'important' genes
gn.imp <- names(imp.temp)[t]
gn.25 <- gn.imp[1:25] # vector of top 25 genes, in order
t <- is.element(rownames(small.eset),gn.25)
sig.eset <- small.eset[t,] # matrix of expression values, not necessarily in order

## Make a heatmap
library(RColorBrewer)
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256)
colnames(sig.eset) <- group
csc <- rep(hmcol[50],30)
csc[group=='C'] <- hmcol[200]
heatmap(sig.eset, scale="row", col=hmcol, ColSideColors=csc)
probe.names <- c('43812','11461','15197','31423')
gene_info <- featureData(gse)[probe.names, c("GENE_SYMBOL")]

# LDA + ROC CURVE ----
#BiocManager::install("pROC")
#BiocManager::install("genefilter")
library("genefilter")

ex3<-ex2[,1:40]
f <- factor(c(rep("control",20), rep("affected",20)))
tt40 <- rowttests(ex3,f)
keepers <- which(tt40$p.value<0.1)
ex4 <- ex3[keepers,]
tex4 <- t(ex4)
dat <- cbind(as.data.frame(tex4),f)
colnames(dat)[ncol(dat)] <- "AFFECTED"
n.controls <- 20
n.affected <- 20
train <- sample(1:(n.controls), (n.controls-5))
test <- setdiff(1:(n.controls),train)
test<- c(test, test+20)
train <- c(train, train+20)
library("MASS")
mod <- lda(AFFECTED ~ ., data=dat, prior = c(0.5,0.5),
           subset = train)
par(mar = c(1, 1, 1, 1))
plot(mod) #projections of the two groups, 4 misclassified
mod.values <- predict(mod, dat[train,])
mod.values$class
par(mfrow=c(1,1))
plot(mod.values$x[,1], ylab=c("LDA Axis"))
text(mod.values$x[,1],
     col=c(as.numeric(dat[train,"AFFECTED"])+10))
preds<-predict(mod, dat[test,])
preds$class
table(as.numeric(preds$class),
      as.numeric(dat[test, "AFFECTED"]) )
library("pROC")
roc_lda <- plot.roc(as.numeric(preds$class),
                    as.numeric(dat[test, "AFFECTED"]) ) #not useful with LDA

# CARET ----
#install.packages("ggplot2")
#BiocManager::install("caret")
library("caret")
# Run algorithms using 10-fold cross validation
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
fit.lda <- train(AFFECTED~., data=dat, method="lda",
                 metric=metric, trControl=control)
fit.rf <- train(AFFECTED~., data=dat, method="rf",
                metric=metric, trControl=control)
results <- resamples(list(LDA=fit.lda, RF=fit.rf))
summary(results)
ggplot(results) + labs(y = "Accuracy")
# Run algorithms using 10-fold cross validation, 5 times
control <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
fit.lda.2 <- train(AFFECTED~., data=dat, method="lda",
                   metric=metric, trControl=control)
fit.rf.2 <- train(AFFECTED~., data=dat, method="rf",
                  metric=metric, trControl=control)
results <- resamples(list(LDA=fit.lda.2, RF=fit.rf.2))
ggplot(results) + labs(y = "Accuracy")


# Lasso ----
ex3<-ex2[,1:40]
dat <- t(ex3)
y <- c(rep(0,20),rep(1,20))
f <- factor(y)
#install.packages("glmnet")
library("glmnet")
fit=glmnet(dat,y,standardize=FALSE,family="binomial")
plot(fit, xvar = "lambda", label=TRUE)
cfit=cv.glmnet(dat,y,standardize=FALSE,family="binomial")
plot(cfit)
coef(cfit, s=cfit$lambda.min)
# repeat analysis but by using train + test sample subsets
n.controls<-20
n.affected<-20
train <- sample(1:(n.controls), (n.controls-5))
test <- setdiff(1:(n.controls),train)
test<- c(test, test+20)
train <- c(train, train+20)
fit=glmnet(dat[train,],y[train],standardize=FALSE,family="binomial")
plot(fit)
cfit=cv.glmnet(dat[train,],y[train],standardize=FALSE,family="binomial")
plot(cfit)
predict(fit,dat[test,], type="class", s= cfit$lambda.min)

# plot ROCR curve
#install.packages("ROCR")
library("ROCR")
pred2 <- predict(fit,dat[test,], type="response", s=cfit$lambda.min)
plot(performance(prediction(pred2, y[test]), 'tpr', 'fpr'))
# compute Area Under the Curve (AUC)
auc.tmp <- performance(prediction(pred2, y[test]),"auc")
auc <- as.numeric(auc.tmp@y.values)

# CARET + Lasso ----
ex3<-ex2[,1:40]
f <- factor(c(rep("control",20), rep("affected",20)))
library("genefilter")
tt40 <- rowttests(ex3,f)
keepers <- which(tt40$p.value<0.1)
ex4 <- ex3[keepers,]
tex4 <- t(ex4)
#tex3 <- t(ex3) # disable feature selection
dat <- cbind(as.data.frame(tex4),f)
colnames(dat)[ncol(dat)] <- "AFFECTED"

#install.packages("caret")
library("caret")
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"
fit.lasso <- train(AFFECTED~., data=dat, method="glmnet",
                   family = "binomial",
                   tuneGrid = expand.grid(alpha = 1,
                                          lambda = seq(0,1,by=0.05)),
                   trControl = control,
                   metric = metric)

plot(fit.lasso)
# comparison with other classification methods
fit.lda <- train(AFFECTED~., data=dat, method="lda",
                 metric=metric, trControl=control)
fit.rf <- train(AFFECTED~., data=dat, method="rf",
                metric=metric, trControl=control)
results <- resamples(list(RF=fit.rf, LDA=fit.lda, Lasso=fit.lasso))
summary(results)
ggplot(results) + labs(y = "Accuracy") 


# RSCUDO ----
dat <- ex2[,1:40]
y <- c(rep(0,20),rep(1,20))
f <- factor(y, labels = c("Control",
                          "Affected"))
library("caret")
set.seed(123)
inTrain <- createDataPartition(f, list = FALSE)
trainData <- dat[, inTrain]
testData <- dat[, -inTrain]
#BiocManager::install("rScudo")
library("rScudo")
trainRes <- scudoTrain(trainData, groups = f[inTrain],
                       nTop = 25, nBottom = 25, alpha = 0.05)
trainRes
upSignatures(trainRes)[1:5,1:5]
consensusUpSignatures(trainRes)[1:5, ]
probe.names.c <- c('17901','3808','4838','5098','23751')
gene_info_c <- featureData(gse)[probe.names.c, c("GENE_SYMBOL")]
probe.names.a <- c('3431','40012','25294','34586','28032')
gene_info_a <- featureData(gse)[probe.names.a, c("GENE_SYMBOL")]
trainNet <- scudoNetwork(trainRes, N = 0.2)
scudoPlot(trainNet, vertex.label = NA)
testRes <- scudoTest(trainRes, testData, f[-inTrain],
                     nTop = 25, nBottom = 25)
testNet <- scudoNetwork(testRes, N = 0.2)
scudoPlot(testNet, vertex.label = NA)
# identify clusters on map
library("igraph")
testClust <- igraph::cluster_spinglass(testNet, spins = 2)
plot(testClust, testNet, vertex.label = NA)
# classification
classRes <- scudoClassify(trainData, testData, N = 0.25,
                          nTop = 12, nBottom = 12,
                          trainGroups = f[inTrain], alpha = 0.5)
caret::confusionMatrix(classRes$predicted, f[-inTrain])

# RSCUDO + CARET ----
dat <- ex2[,1:40]
y <- c(rep(0,20),rep(1,20))
f <- factor(y, labels = c("Control",
                          "Affected"))
library("caret")
set.seed(123)
inTrain <- createDataPartition(f, list = FALSE)
trainData <- dat[, inTrain]
testData <- dat[, -inTrain]
model <- scudoModel(nTop = (2:6)*5, nBottom = (2:6)*5,
                    N = 0.25)
control <- caret::trainControl(method = "cv", number = 5,
                               summaryFunction =
                                 caret::multiClassSummary)
cvRes <- caret::train(x = t(trainData), y = f[inTrain],
                      method = model,
                      trControl = control)
# plot map
testRes <- scudoTest(trainRes, testData, f[-inTrain],
                     cvRes$bestTune$nTop,
                     cvRes$bestTune$nBottom5)
testNet <- scudoNetwork(testRes, N = 0.2)
scudoPlot(testNet, vertex.label = NA)

classRes <- scudoClassify(dat[, inTrain], dat[, -inTrain],
                          0.25,
                          cvRes$bestTune$nTop,
                          cvRes$bestTune$nBottom, f[inTrain], alpha = 0.05)
caret::confusionMatrix(classRes$predicted, f[-inTrain])


# Functional annotation ----
consensusUpSignatures(trainRes)[1:25, ]
up.c <- consensusUpSignatures(trainRes)[1:25, ]['Control']
probe.names.c <- unlist(up.c)
gene_info_c <- featureData(gse)[probe.names.c, c("GENE_SYMBOL")]
up.a <- consensusUpSignatures(trainRes)[1:25, ]['Affected']
probe.names.a <- unlist(up.a)
gene_info_a <- featureData(gse)[probe.names.a, c("GENE_SYMBOL")]


# network ----
library(limma)

group <- c(rep("Controls", 20), rep("Affected",30))
design <- model.matrix(~0+group)
colnames(design) <- c("Controls","Affected")
rownames(design)<- colnames(ex2)

fit <- lmFit(ex2, design)
cont.matrix <- makeContrasts(contrasts = "Affected-Controls", levels=design) # i want the up and down reg in tumors respect to contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
stats<-topTable(fit2, number=50)
symbols <- c()
for (i in rownames(stats)){
  row_index <- which(rownames(ex2) == i)
  symbols <- c(symbols,gse[["GSE49126_series_matrix.txt.gz"]]@featureData@data[["GENE_SYMBOL"]][row_index])
}
dif_genes <- stats
dif_genes$Gene.symbol <- symbols
dif_genes$adj.P.Val <- stats$adj.P.Val
dif_genes <- dif_genes[order(dif_genes$adj.P.Val), ]
dif_genes <- na.omit(dif_genes)
dif_genes <- dif_genes[,c('Gene.symbol','logFC','adj.P.Val')]
write.csv(dif_genes, file = "pathfindR1.csv", row.names = FALSE)
