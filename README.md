# Machine Learning Models for Cancer Prediction Based on Proteomics Data
Author: John Paul Wick

Date: 16 September 2024

# Background

Machine learning offers powerful potential in predicting various outcome labels from high-dimensional data. This makes such models very useful tools for predictions based on omics data in clinical or research settings. 
To train such models, large volumes of well-organized data must be available. This is why in bioinformatics the creation and maintenance of large data repositories is crucial for the development and improvement of pipelines and tools. 

The Cancer Genome Atlas (TCGA) is one such archive of a collection of large datasets which each describe a different type of cancer, organized into projects. 
In this project, I present an R pipeline for developing and testing machine learning models to predict and test presence or absence of breast cancer, using 226 samples downloaded from the TCGA-BRCA project (113 control samples and 113 cancer samples). 
There are four model types included: logistic regression (LR), support vector machine (SVM), neural network (NN), and random forest (RF).

# Pipeline Description

See the Final-Project.md document in this repository for the full pipeline analysis, including plots and results.

The pipeline loads in a SummarizedExperiment object containing the pre-downloaded data from TCGA-BRCA. The data are normalized and split into training data (50 controls and 50 cancer samples) and testing data (the remaining 126 samples).
The training-test split was made early such that the testing data could not in any way influence model development. 

I first narrowed the dataset down to only protein-coding genes, as a way to simplify the model building. Further analysis could be done with other types of genes. 
Through differential expression analysis using DESeq2, variance analysis, t-tests, and correlation analysis, I narrowed the protein-coding genes down to a small set of 29 genes which contained the following characteristics:

1. The 29 genes were in the top percentile of all the nearly 20,000 protein coding genes for having the highest variance

2. The 29 genes all had significant adjusted p-values in the DESeq2 analysis.

3. The 29 genes all had significant t-test p-values between control and cancer samples.

4. The 29 genes were all uncorrelated with themselves (strongest absolute value of correlation coefficient = 0.37), as well as being uncorrelated with the sample metadata attributes (besides the control/cancer annotation).

PCA analysis was also performed in the pipeline, revealing striking clusters when colored by sample type. 
While PCA is most useful during unsupervised exploration of unlabeled data, the results indicated that a large portion of the variance in the data was explained by sample type.
The clusters were also colored according to a few other sample metadata attributes, without strong cluster separation. 

Finally, once the set of 29 genes was selected, predictive models were trained and cross-validated to optimize both numbers of features included as well as other parameters (SVM gamma value, NN number of hidden layers, and RF number of features to sample per tree). 

# Results

The four final models were used to predict sample type (control or cancer) in the 126 testing samples set aside, from the expression of the features selected for model training. 
Assessment of model performance was based on area under the curve (AUC) of Receiver Operating Characateristic (ROC) curves. All four models were excellent, with SVM outperforming the others (LR AUC = 0.92, SVM AUC = 0.98, NN AUC = 0.83, and RF AUC = 0.97). 
For more details and figures on the model results and a discussion, see the final sections of Final-Project.md. 
