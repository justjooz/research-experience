# This R script is for obtaining a dataframe from the GEO database using GEOQuery
#======================================================================================

# (1) Set working directory

      setwd("D:/Code/RE/My R scripts")

# (2) Obtain GEO GSE file
      
      library(GEOquery) # always use the library() function first in order to use a certain package
      
      # GSE22544 <- getGEO("GSE22544", destdir = '.') #run this to get the saved dataset file
      GSE22544 <- getGEO(filename = 'GSE22544_series_matrix.txt.gz')
      GSE22544
      
      View(GSE22544)
      class(GSE22544)
      
      View(exprs(GSE22544))
      
# (3) Putting pheno data and feature data into respective vectors
      pheno_GSE22544 <- phenoData(GSE22544) # pheno data
      feature_GSE22544 <- fData(GSE22544) # feature data
      
# (4) Obtaining expression data from feature data
      exprs_GSE22544 <- exprs(GSE22544) # exprs() function applied to GSE22544
      class(exprs_GSE22544) # shows as matrix
      exprs_GSE22544 <- as.data.frame(exprs_GSE22544) # change matrix to data.frame
      View(exprs_GSE22544) # GSM on the columns, 'at' on the rows
      
# (5) Putting GENE Symbols into a df
      GENE_Symbols <- (feature_GSE22544 [11])
      rownames(GENE_Symbols) <- NULL # removes probe IDs as rownames
      class(GENE_Symbols)
      
# (6) cbind to combine GENE_Symbols df with expression df
      bound_GSE22544 <- cbind(GENE_Symbols, exprs_GSE22544)
      rownames(bound_GSE22544) <- NULL
      
      View(bound_GSE22544)
      ncol(bound_GSE22544)
      
# (7) Creating a class information vector and factor
      
      # using gsub to change the class labels
      pheno_GSE22544@data$characteristics_ch1 <- gsub("tissue:", "", pheno_GSE22544@data$characteristics_ch1)
      pheno_GSE22544@data$characteristics_ch1 <- gsub("breast cancer", "IDC", pheno_GSE22544@data$characteristics_ch1)
      pheno_GSE22544@data$characteristics_ch1 <- gsub("normal breast", "normal", pheno_GSE22544@data$characteristics_ch1)
      
      # making a factor for use later in bootstrap step
      class_factor <- as.factor(pheno_GSE22544@data$characteristics_ch1)
      levels(class_factor)
      levels(class_factor) [1]
      levels(class_factor) [2]
      
      class_factor == levels(class_factor) [1] #  [1]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
      
      
      class_factor == levels(class_factor) [2] #  [1] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
      
      
      # making a vector for insertion into exprs dataframe later in bootstrap step
      class_vector <- c(pheno_GSE22544@data$characteristics_ch1)
      
# (8) Creating a gene symbols vector      
      gene_symbol_GSE22544 <- (feature_GSE22544 [11])
      rownames(gene_symbol_GSE22544) <- NULL
      # View(gene_symbol_GSE22544)
      
      gene_symbol_GSE22544 <- as.vector(gene_symbol_GSE22544[[1]])
      is.vector(gene_symbol_GSE22544)
      
      # writing into a csv file
      # write.csv()
      
write.csv(bound_GSE22544, file = "D:/Code/RE/My R scripts/bound_GSE22544.csv") # writes bound df with gene symbol and expression data
write.csv(gene_symbol_GSE22544, file = "D:/Code/RE/My R scripts/gene_symbol_GSE22544.csv") # writes gene symbol vector
write.csv(class_factor, file = "D:/Code/RE/My R scripts/class_factor.csv") # writes class labels factor
write.csv(class_vector, file = "D:/Code/RE/My R scripts/class_vector.csv")
          
