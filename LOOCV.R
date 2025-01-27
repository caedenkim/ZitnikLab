install.packages("caret")

library(caret)

#specify the cross-validation method
ctrl <- trainControl(method = "LOOCV")

cd8positivecytotcell$patient <- factor(cd8positivecytotcell$patient)
#fit a regression model and use LOOCV to evaluate performance
model <- train(response ~ B2M	+ CD44 + VIM + PTPRC + HLA.C + TSC22D3 + CREM + HSP90AA1 + CD2 + HLA.A + CFLAR + NFKBIA + NR3C1 + CXCR4 + ITK + MCL1 + HMGB1 + HLA.E + PRDM1 + PTGER4, data = cd8positivecytotcell, method = "glm", trControl = ctrl) #splitting data and training model at the same time
model #how well your model does
model2 <- train(response ~ , data = maturenktcell, method = "glm", trControl = ctrl)
model2
modelall <- train(response ~ , data = maturenktcell, method = "glm", trControl = ctrl)
modelall

fitControl = trainControl(method = "repeatedcv", number = 10, repeats = 20, summaryFunction = twoClassSummary, 
                          classProbs = TRUE, savePredictions = T) #use this control instead of line 6
model = train(response ~ B2M	+ CD44 + VIM + PTPRC + HLA.C + TSC22D3 + CREM + HSP90AA1 + CD2 + HLA.A + CFLAR + NFKBIA + NR3C1 + CXCR4 + ITK + MCL1 + HMGB1 + HLA.E + PRDM1 + PTGER4, data = cd8positivecytotcell, method = 'glm', metric = 'ROC', trControl = fitControl)

#view summary of LOOCV               
print(model)


