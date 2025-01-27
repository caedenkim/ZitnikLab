
library(aod)
library(ggplot2)


cd8positivecytotcell <- read.csv('/Users/caeden/Desktop/cd8-positive, alpha-beta cytotoxic t cell.csv')
cd8positivecytotcell
rownames(cd8positivecytotcell) <- cd8positivecytotcell$patient
# each row = patient (total 45 rows) 
# each column = gene (total 20 columns)
# each patient has a vector of 20 values represented
View(cd8positivecytotcell)
head(cd8positivecytotcell) # view the first few rows of data

as.numeric(as.vector(cd8positivecytotcell[1,])) # vector of patient PB-16-002

summary(cd8positivecytotcell) # 5 number summary for each gene
sapply(cd8positivecytotcell, sd) # standard deviation of each gene


# patient_df.txt$patient <- factor(patient_df.txt$patient) # patient is a categorical variable
mylogit <- glm(response ~ B2M	+ CD44 + VIM + PTPRC + HLA.C + TSC22D3 + CREM + HSP90AA1 + CD2 + HLA.A + CFLAR + NFKBIA + NR3C1 + CXCR4 + ITK + MCL1 + HMGB1 + HLA.E + PRDM1 + PTGER4, data = cd8positivecytotcell, family = "binomial")
summary(mylogit)
mylogit <- glm(response ~ HSP90AA1 + B2M + CXCR4 + HLA.E + KLRB1 + NFKBIA + FOS + IRF1 + TAGAP + CCL4 + DUSP1 + MCL1 + JUN + HLA.A + HLA.C + CASP8 + CD55 + CTNNB1 + TSC22D3 + IL2RB, data = patient_df.txt, family = "binomial")
# description of linear predictor and description of error distribution ?
# binomial distribution ? 
summary(mylogit) # see Copy of ML task description google doc




# split data into train, validation, and test set





# load library pCOR
library(pROC)

# predicted data
prediction <- predict(mylogit, cd8positivealphabetacytotcell, type="response")
prediction

# create roc curve
roc_object <- roc(cd8positivecytotcell$patient, prediction)

# calculate area under curve
auc( roc_object )

