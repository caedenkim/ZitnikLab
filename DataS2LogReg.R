library(aod)
library(ggplot2)


patient_df.txt2 <- read.csv('/Users/caeden/Desktop/datas2patient_df.txt.csv')
patient_df.txt2
rownames(patient_df.txt2) <- patient_df.txt2$patient
# each row = patient (total 45 rows) 
# each column = gene (total 20 columns)
# each patient has a vector of 20 values represented
View(patient_df.txt2)
head(patient_df.txt2) # view the first few rows of data

as.numeric(as.vector(patient_df.txt2[1,])) # vector of patient PB-16-002

summary(patient_df.txt2) # 5 number summary for each gene
sapply(patient_df.txt2, sd) # standard deviation of each gene


# patient_df.txt$patient <- factor(patient_df.txt$patient) # patient is a categorical variable
mylogit <- glm(response ~ CTLA4 + CD274 + CD28 + IL2 + IL10 + FOXP3 + IL6 + CD80 + PTPN11 + LCK + ITGAM + PTPN6 + CD40 + AKT1 + PDCD1 + LAG3 + ITGAX + CSF2 + IFNG , data = patient_df.txt2, family = "binomial") # all genes in patient data  
summary(mylogit)
# description of linear predictor and description of error distribution ?
# binomial distribution ? 
summary(mylogit) # see Copy of ML task description google doc




# split data into train, validation, and test set





# load library pCOR
library(pROC)

# predicted data
prediction <- predict(mylogit, patient_df.txt, type="response")
prediction

# create roc curve
roc_object <- roc(patient_df.txt$patient, prediction)

# calculate area under curve
auc( roc_object )