lung_cdata <- read.csv("Normalized_count_Data_final.csv", row.names = 1)
#View(lung_cdata)
lung_meta <- data.frame(subtype=factor(c(rep("IrregularSmoker", 20), rep("RegularSmoker",20))))
rownames(lung_meta) <- colnames(lung_cdata)
View(lung_meta)

install.packages("caret")
install.packages("DALEX")
install.packages("pROC")

library(caret)
library(pROC)
library(DALEX)

#Unsupervised machine learning knn
table(lung_meta$subtype) #how many samples do we have under each subset
#boxplot(lung_cdata[,], las = 2)
#par(oma = c(1,0,0,0))
#boxplot(log10(lung_cdata[,]+1), ylim = c(0,10), las = 2)

#data preprocessing for selecting top variable genes
all.trans <- data.frame(t(lung_cdata)) #transposed matrix is converted to df
#View(all.trans)
dim(all.trans)
SDs = apply(all.trans, 2, sd) #calculate SD over column
#View(SDs)
topPreds = order(SDs, decreasing = T)[1:15000]
all.trans = all.trans[,topPreds]
View(all.trans)
all.trans[1:5,1:5]
#final dataset for ML
all.trans <- merge(all.trans, lung_meta, by ="row.names") #creates a rownames column
#View(all.trans)
rownames(all.trans) <- all.trans$Row.names #making the sample ids rownames again
all.trans <- all.trans[,-1] #dropping the row.names column

#pre-processing step 2
#remove near zero variation in a column
dim(all.trans)
all.zero<- preProcess(all.trans, method = 'nzv', uniqueCut = 5) #out of 40 samples, if at least 5 samples doesn't have different count, it will be removed
all.trans <- predict(all.zero, all.trans)
dim(all.trans) #no gene got filtered

#sd 1, mean=0, centralize
all.center <- preProcess(all.trans, method = "center")
all.trans <- predict(all.center, all.trans) #no gene got filtered

#remove highly correlated
all.corr <- preProcess(all.trans, method = "corr", cutoff =0.75)
all.trans <- predict(all.corr, all.trans) #1438 genes remained

#splitting the dataset (70-30)
intrain <- createDataPartition(y=all.trans$subtype,p=0.7)[[1]]
intrain
#separate the test and training 
train.lung_data <- all.trans[intrain,]
test.lung_data <- all.trans[-intrain,]
dim(train.lung_data)
dim(test.lung_data)
train.lung_data$subtype
test.lung_data$subtype
#training and testing
#using repeated cross validation method
ctrl.lung_data <- trainControl(method = 'repeatedcv', number = 3, repeats = 3)
set.seed(123)  # Set a random seed for reproducibility
knn.lung_data <- train(subtype~.,
                       data = train.lung_data,
                       method = "knn",
                       trControl = ctrl.lung_data,
                       tuneGrid = data.frame(k=1:18))
#the best k
knn.lung_data$bestTune
trainPred <- predict(knn.lung_data, newdata = train.lung_data)
testPred <- predict(knn.lung_data, newdata = test.lung_data)
dim(train.lung_data)
dim(test.lung_data)

##interpretation of model accuracy
#confusion matrix
confusionMatrix(trainPred, train.lung_data$subtype)
confusionMatrix(testPred, test.lung_data$subtype) #58% accuracy in test data prediction

#determine variable importance
explainer.lung_data <- explain(knn.lung_data,
                               label = "knn",
                               data = train.lung_data,
                               y = as.numeric(train.lung_data[,ncol(train.lung_data)]))
importance.lung_data <- feature_importance(explainer.lung_data, n_sample = 50, type = "difference")
head(importance.lung_data$variable)
tail(importance.lung_data$variable)
plot(importance.lung_data)

#using random forest classification
rf.ctrl <- trainControl(method = 'cv') #or bootstrap can be used
rf.lung_data <- train(subtype~.,
                      data = train.lung_data,
                      method = "ranger",
                      trControl = rf.ctrl,
                      importance = 'permutation',
                      tuneGrid = data.frame(mtry=100,
                                            min.node.size=1,
                                            splitrule = "gini"
                      ))
rf.lung_data$finalModel$prediction.error #error rate
plot(varImp(rf.lung_data), top = 20)
