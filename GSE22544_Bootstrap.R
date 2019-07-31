# this R script shows the process to obtain significant genes, binarize them, and do a 1000x bootstrap
# followed by using a distance matrix/heat map 
# and confusion matrix to obtain the F-score, recall, precision

# Set working directory
setwd("D:/Code/RE/My R scripts")

# ==========================================================================================================
# Reading GSE22544 dataframe into script

new_df_GSE22544 <- read.csv("new_df_GSE22544.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1) # read exprs dataframe
new_df_GSE22544 <- cbind(new_df_GSE22544[, c(4,8,9,13)], new_df_GSE22544[,-c(4,8,9,13)]) # rearranges the normal classes to the first 4 columns
new_df_GSE22544 <- na.omit(new_df_GSE22544) # removes NAs which causes probs with for loop

# (2) Method 1: using for loops to produce binary matrix (but leads to error that says "data are essentially constant")
# ===========================================================================================================

significantgenes <- c() # to initialise the vector
boot_list <- list()
library(progress)
pb <- progress_bar$new(total = 1000)

    for(j in 1:1000){ # 1000 is the number of bootstraps ## cannot be run if done on bigger datasets
      significantgenes <- c()
      boot_normal <- sample(new_df_GSE22544[, 1:4], size=4, replace=TRUE)
      boot_IDC <- sample(new_df_GSE22544[,5:20], size=4, replace=TRUE)
      
     for(i in 1:2000){ # 9994 is the number of genes # if set to 9995, "data are essentially constant" error
      ttest_1  <- try(t.test(boot_normal[i,], boot_IDC[i,]))
      if (ttest_1$p.value <= 0.05) # if p value is <= 0.05(less than or equals),
      {
        significantgenes <- append(significantgenes, 1) # assign the genes (rows) as 1 (significant)
      }
      else # if p value is > 0.05,
      {
        significantgenes <- append(significantgenes, 0) # assign the genes (rows) as 0 (not significant)
      }
     }
      boot_list <- append(boot_list, list(significantgenes))
      pb$tick()
      Sys.sleep(1 / 1000)
    }

# (2) Method 2: Using "genefilter" package's "rowttests" function (runs faster than nested loop & w/o error)
# =================================================================================================

# -------- Prep ----------------
#BiocManager::install("genefilter")
library(genefilter)
class_factor_2 <- as.factor(c(rep("normal", 4), rep("IDC", 4))) # creating a class factor (normal & IDC)

# ------ rowttests function to generate binary matrix -----------
significantgenes_2 <- c()
boot_list_2 <- list()
library(progress)
pb <- progress_bar$new(total = 1000)

for (i in 1:1000){
  m_normal <- as.matrix(sample(new_df_GSE22544[,1:4], size = 4, replace = T)) # samples 4 normal patients into a matrix
  m_IDC <- as.matrix(sample(new_df_GSE22544[,5:20], size = 4, replace = T)) # samples 4 normal patients into a matrix
  m_bind <- cbind(m_normal, m_IDC) # has 8 columns, first 4 columns are sampled normal, last 4 columns are sample IDC
  
  ttest_2 <- rowttests(m_bind, class_factor_2) # rowttests tests all 9994 rows(genes) w/o error
  significantgenes_2 <- as.numeric(ttest_2$p.value < 0.05) # <0.05 forms a boolean output, and it is changed to numeric by as.numeric.
  boot_list_2 <- append(boot_list_2, list(significantgenes_2)) # append each sampling as a list to a list
  
  pb$tick() # for progress bar
  Sys.sleep(1 / 1000) # for progress bar
}
# -------rowttests function to generate binary matrix -----------

# (3) Convert list into matrix
# ==================================================================================================
boot_mat_2 <- matrix(unlist(boot_list_2), ncol = 1000, byrow = FALSE) # converting list into matrix
# dim(boot_mat_2)

# (4) Computing rowsums of each gene to see significance out of 1000
# ===========================================================================
sum_vect <- rowSums(boot_mat_2)
tail(sort(sum_vect), 5) 

hist(sum_vect, freq=FALSE, breaks = 200, 
     main = "Distribution of Significant Samplings in Breast Cancer Genes", 
     xlab = "Significant Samplings", 
     col = "orange",
     xlim=c(0,1000),
     ylim = c(0, 0.020)) # to visually inspect to determine point threshold of significance

# `sum_sig` is assigned the numeric output of the boolean of `i>180`
sum_sig <- as.numeric(sum_vect>180)

# (5) Creating confusion matrix
# =======================================================
first_sample <- as.factor(boot_mat_2[,1]) # actual
sum_sig <- as.factor(sum_sig) # observation

# install.packages("caret")
library(caret) # for confusion matrix function

confusionMatrix(sum_sig, first_sample)
            # Confusion Matrix and Statistics
            
            #             Reference
            # Prediction    0         1
            #           0 9025 (TP)  324 (FP)
            #           1  546 (FN)   99 (TN)
            # 
            # Accuracy : 0.9129          
            # 95% CI : (0.9072, 0.9184)
            # No Information Rate : 0.9577          
            # P-Value [Acc > NIR] : 1               
            # 
            # Kappa : 0.1415          
            # 
            # Mcnemar's Test P-Value : 6.752e-14       
            # 
            # Sensitivity : 0.9430          
            # Specificity : 0.2340          
            # Pos Pred Value : 0.9653          
            # Neg Pred Value : 0.1535          
            # Prevalence : 0.9577          
            # Detection Rate : 0.9030          
            # Detection Prevalence : 0.9355          
            # Balanced Accuracy : 0.5885          
            # 
            # 'Positive' Class : 0               

conf_mat <- (confusionMatrix(sum_sig, first_sample))$table

# (6) Calculating metrics: Precision, Recall, F-Score
# ========================================================================

# Precision: TP/(TP+FP):
precision <- conf_mat[1,1]/sum(conf_mat[1,1:2])
precision # [1] 0.9653439

# Recall: TP/(TP+FN):
recall <- conf_mat[1,1]/sum(conf_mat[1:2,1])
recall # [1] 0.9429527

# F-Score: 2 * precision * recall /(precision + recall):
f_score <- 2 * precision * recall / (precision + recall)
# [1] 0.9540169

# Metrics
# Precision: 0.9653439
# Recall: 0.9429527
# F-score: 0.9540169

# Visualising the Confusion Matrix
cm <-confusionMatrix(sum_sig, first_sample)

# ?confusionMatrix
draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2.5)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 438, 'Positive', cex=2)
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 438, 'Negative', cex=2)
  text(125, 370, 'Predicted', cex=2, srt=90, font=2)
  text(245, 450, 'Actual', cex=2, font=2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(143, 400, 'Positive', cex=2, srt=90)
  text(143, 335, 'Negative', cex=2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=4, font=2, col='black')
  text(195, 335, res[2], cex=4, font=2, col='black')
  text(295, 400, res[3], cex=4, font=2, col='black')
  text(295, 335, res[4], cex=4, font=2, col='black')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  

draw_confusion_matrix(cm)

# (7) Create user-defined function that calculates jaccard coefficient (intersection over union)
# =========================================================================
jaccard_fun <- function (x,y) {   
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}

# Initialising dataframe for imputation
jaccard_df <- as.data.frame(matrix(NA,nrow = 1000, ncol = 1000))
names(jaccard_df) <- paste0('S', 1:1000)
rownames(jaccard_df) <- paste0('S', 1:1000)

library(progress)
pb <- progress_bar$new(total = 1000)

# for loop that will produce the distance heatmap/
for (r in 1:1000) {
  for (c in 1:1000) {
    if (c == r) { # if rows iteration is the same as column iteration,
      jaccard_df[r,c] = 1 # assign as 1
    } else if (c > r) { # if not then when columns is more than rows, add the variables of rows and
      jaccard_df[r,c] <- jaccard_fun(boot_mat_2[,r], boot_mat_2[,c]) # replace with list from above
    }
  }
  pb$tick() # for progress bar
  Sys.sleep(1 / 1000) # for progress bar
}

write.csv(boot_mat_2, file = "D:/Code/RE/My R scripts/boot_mat_2.csv", row.names = T)
write.csv(jaccard_df, file = "D:/Code/RE/My R scripts/jaccard_df.csv", row.names = T)
