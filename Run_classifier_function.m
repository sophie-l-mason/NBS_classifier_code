%% All the inputs for the classifier
addpath(genpath('NBS_Fast')) %adds the NBS toolbox
%load('Data_random.mat'); %load FNs N_ROI x N_ROI x N_Subj
Data=cat(3,rand(70,70,16)+2,rand(70,70,22)+1.9);
load('Design_matrix.mat'); %load design matrix N_Subjx2
contrast=1; %select 1 for [1,-1] or -1 for [-1,1] 
N_perms=50; %Number of permutations used to create the p_val e.e. 1000 or 5000
alpha=0.05; %significance threshold level
path='C:\Users\sophi\Documents\NBS_github' %this is a path for which a folder is created results to be saved

% The classifier
rng(1)%sets the seed used for the random seed generator 'shuffle' uses a different random seed each time
[class_label,Accuracy,P_val, N_edges,Step]=Classifier_function(Data,Design_matrix,contrast,N_perms,alpha,path)

%% Aditional inputs for the varying the t-statistic threshold in the classifier
threshold_min=0.01; %minimum t-statistic threshold to consider
threshold_max=4.5 %maximum t-statistic threshold to consider
Title='Random data' %title for the figure
Y_Label='[ECP $>$ LCP]' %y label for the figure
figure_index=1 %the number for figure window  
Data=cat(3,rand(70,70,16)+2,rand(70,70,22)+1.9)
% The classifier: optimised for varying the t-statistic threshold
rng(1)%sets the seed used for the random seed generator 'shuffle' uses a different random seed each time
[Results,class_label]=NBS_vary_t_stat(Data,Design_matrix,contrast,N_perms,alpha,path,threshold_min,threshold_max,Title,Y_Label,figure_index)

