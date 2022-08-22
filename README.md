# NBS_classifier_code
A repository containing the code used in 'Classification of Human Chronotype Based on fMRI Network-Based Statistics'.

The Run_classifier_function.m runs the classifier called Classifier_function.m and the file NBS_vary_t_stat.m, which is the optimised version of the classifier for when the t-statistic threshold is varied. This file will call the matlab files Data_random and Design_matrix. Data_random is an array 70x70x38 (N_ROIs x N_ROIs x N_Subj) created using random numbers. Design_matrix (N_Subj x 2) shows how the design matrix must be saved. The results of this data are meaningless but enabled you to run the data to test it is working.

The classifier is called Classifier_function.

The classifier when varying the t-statistic is called NBS_vary_t_stat.

Tikhonov Partial correlation matrices
The file Tikhonov_matrix_create can be used to create Tikhonov parital correlation matrices.
This relies on invChol_mex.c (https://www.mathworks.com/matlabcentral/fileexchange/34511-fast-and-accurate-symmetric-positive-definite-matrix-inverse-using-cholesky-decomposition)
As this is a .c file it may need to be compiled before it works this is done by ensuring the path with the file invChol_mex.c is on your file path and then running the lines:

mex -setup

mex -O invChol_mex.c

This file requires the fMRI timeseries for each region of interest per person.

The other files are supporting files, with charpath.m, distance_wei.m and weight_conversion from the Brain Connectivity Toolbox and NBS_Fast is an editted version of NBS. This editted version removes some of the displayed text. This speeds up the NBS methodology, which is important when NBS is being repeated many times for each test subject. In addition, NBS_run will produce an output directly rather than needing to make nbs a global function in order to access information.

Note that in all cases the results presented in the paper are calculated using the random seed 1.

