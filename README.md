# NBS_classifier_code
A repository containing the code used in 'Classification of Human Chronotype Based on fMRI Network-Based Statistics'.
The Run_classifier_function runs can run both the classifier as well as the optimised version of the code for when the t-statistic threshold is varied. This file will call the matlab file Data_random. This is an array 70x70x38 (N_ROIs x N_ROIs x N_Subj) created using random numbers. The results of this data are meaningless but enabled you to run the data to test it is working. In addition, the file Design_matrix (N_Subj x 2). This can be viewed to check you are inputting the correct information.
The classifier is then called Classifier_function.
The classifier when varying the t-statistic is called NBS_vary_t_stat.

#Tikhonov Partial correlation matrices
The file ... can be used to create Tikhonov parital correlation matrices
This relies on invChol_mex.c (https://www.mathworks.com/matlabcentral/fileexchange/34511-fast-and-accurate-symmetric-positive-definite-matrix-inverse-using-cholesky-decomposition)
As this is a .c file it may need to be compiled before it works this is done by ensuring the path with the file invChol_mex.c is on your file path and then running the lines:
mex -setup
mex -O invChol_mex.c
