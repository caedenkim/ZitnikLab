library(aod)
library(ggplot2)


maturenktcell <- read.csv('/Users/caeden/Desktop/maturenktcell.csv')
maturenktcell
rownames(maturenktcell) <- maturenktcell$patient
# each row = patient (total 45 rows) 
# each column = gene (total 20 columns)
# each patient has a vector of 20 values represented
View(maturenktcell)
head(maturenktcell) # view the first few rows of data

as.numeric(as.vector(maturenktcell[1,])) # vector of patient PB-16-002

summary(maturenktcell) # 5 number summary for each gene
sapply(maturenktcell, sd) # standard deviation of each gene


# patient_df.txt$patient <- factor(patient_df.txt$patient) # patient is a categorical variable
mylogit <- glm(response ~ HSP90AA1 + B2M + CXCR4 + HLA.E + KLRB1 + NFKBIA + FOS + IRF1 + TAGAP + CCL4 + DUSP1 + MCL1 + JUN + HLA.A + HLA.C + CASP8	+ CD55 + CTNNB1 + TSC22D3 + IL2RB, data = maturenktcell, family = "binomial") # all genes in patient data  
summary(mylogit)
# description of linear predictor and description of error distribution ?
# binomial distribution ? 
summary(mylogit) # see Copy of ML task description google doc




# split data into train, validation, and test set





# load library pCOR
library(pROC)

# predicted data
prediction <- predict(mylogit, maturenktcell, type="response")
prediction

# create roc curve
roc_object <- roc(maturenktcell$patient, prediction)

# calculate area under curve
auc( roc_object )
