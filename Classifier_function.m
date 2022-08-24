function [class_label,Accuracy,P_val, N_edges,step]=Classifier_function(Data,Design,contrast,N_perms,alpha,path,random_seed)
% Classifier_function.m
% Sophie Mason - University of Birmingham, 2022
%
% Copyright (c) 2022, Brunno Machado de Campos &  Sophie Mason
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% A function which calculates the leave one out cross validated accuracy
% INPUT:
% Data = NROIxNROIxN_Subj FC matrix
% Design = N_Subjx2 design matrix
% Label_true = N_Subj x 1 i.e.[zeros(16,1);ones(22,1)]; the true classification labels for the subjects
% contrast: 1 = the contrast [1,-1]
%          -1 = the contrast [-1,1]
% N_perms = the number of permuations used to calculate the FWE corrected
%           p-values
% alpha = the threshold for significance below which the dysconnected network is significant
% path = a path given as a char i.e. '~\documents\NBS_store'
% random_seed = the value for the random seed 'shuffle' makes it random while a scalar i.e. 1 will generate the same values each time.
% OUTPUT:
%class_label = N_Subjx1 vector storing the class label assigned by the
%              classifier to each of the subjects
% Accuracy = 1x1 value for the accuracy of the classifier
% P_val = N_Subjx4 store of the p-values for each dysconnected network      
%         created in the order
%         [P_val_E_tE,P_val_L_tL,P_val_E_tL,P_val_L_tE].
%         If the value is NaN then the dysconnected networks was created
%         and it was non-significant. 
%         If the value is -99 then the dysconnected network was not created
% N_edges = N_Subjx4 store of the number of edges for each dysconnected network      
%           created in the order
%           [N_edges_E_tE,N_edges_L_tL,N_edges_E_tL,N_edges_L_tE.
%           If the value is NaN then the dysconnected networks was created
%           and it was non-significant. 
%           If the value is -99 then the dysconnected network was not created

N_Subj=size(Data,3); %size for a NROIxNROIxNSubj FC matrix
N_early=sum(Design(:,1)); %Number of ECPs
disp(['Number of group1: ' num2str(N_early)])
N_late=sum(Design(:,2)); %Number of LCPs
disp(['Number of group2: ' num2str(N_late)])
Label_true=[zeros(N_early,1);ones(N_late,1)]; %true labels assigned from design matrix

%If 99 appears then you know that it has not been evaluated i.e. the
%decision tree did not need to evaluate it
class_label_1=ones(N_Subj,3)*99; %store for class label
N_edges_E_tE=ones(N_Subj,1)*-99; %store for number of edges
N_edges_L_tE=ones(N_Subj,1)*-99; %store for number of edges
N_edges_E_tL=ones(N_Subj,1)*-99; %store for number of edges
N_edges_L_tL=ones(N_Subj,1)*-99; %store for number of edges
P_val_E_tE=ones(N_Subj,1)*-99; %store for p values
P_val_L_tE=ones(N_Subj,1)*-99; %store for p values
P_val_E_tL=ones(N_Subj,1)*-99; %store for p values
P_val_L_tL=ones(N_Subj,1)*-99; %store for p values
threshold_early=ones(N_Subj,1)*-99; %store for threshold_early
threshold_late=ones(N_Subj,1)*-99; %store for threshold_late
class_label=ones(N_Subj,1)*-99; %store for class_label
step=ones(N_Subj,1)*-99; %store for step of the classifier when the label is given
sig_E_tL=ones(N_Subj,1)*-99; %store for sig_E_tL
sig_L_tE=ones(N_Subj,1)*-99; %store for sig_L_tE

LCC_70Nodes=1;% if 1 this will ensure you find the threshold where the
%largest connected component has 70 nodes, if 0 this will ensure you find
%the threshold where the t-statistic matrix goes from 70 to 69 nodes.

mkdir([path '\Store'])
for k=1:size(Data,3)
    rng(random_seed)
    disp(['Loop ' num2str(k) ' out of ' num2str(size(Data,3))]) %Shows how far the loop you are
    clearvars nbs 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--Label the subj an Early--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loo_Design_k=[Design((1:k-1),:); [1 0]; Design((k+1:end),:)]; %label the kth subject an early
    save([path '\Store\loo_Design_early' num2str(k)],'loo_Design_k') %saves the design matrix with k removed
    save([path '\Store\loo_Evening_early' num2str(k)],'Data') %saves the FC matrix with k removed
    UI.design.ui=[path '\Store\loo_Design_early' num2str(k) '.mat']; %uploads the design matrix with k removed into NBS
    UI.matrices.ui=[path '\Store\loo_Evening_early' num2str(k) '.mat']; %uploads the FC matrix with k removed into NBS
    UI.method.ui='Run NBS';
    UI.test.ui='t-test'; %type of statistical test
    UI.size.ui='Intensity'; %considers the sum of the t-statistics weights
    UI.perms.ui=num2str(N_perms); %number of permuations for FWE test
    UI.alpha.ui='1'; %value of alpha - set to 1 so all information for all edges is collected and then alpha is reduced later
    if contrast==1 %type of contrast used
        UI.contrast.ui='[1,-1]'; 
    elseif contrast==-1
        UI.contrast.ui='[-1,1]';
    end
    UI.exchange.ui='';
    UI.node_coor.ui='';
    UI.node_label.ui='';
    UI.thresh.ui='1'; % a value for the t-statistic threshold is needed in order to extract information required later    
    nbs=NBSrun(UI,''); %this runs NBS
    %global nbs; %this makes the information stored in nbs a global variable
    
    %find the percolation threshold early
    t_stat_mat=nbs.NBS.test_stat; %obtain the t statistic matrix   
    ordered_t_test=sort(t_stat_mat(:),'descend'); %All the t-stastics calculated for subject k sorted in order 
    %this considered each t-statistic in turn and stops when it finds once
    %that makes the dsyconnected network lose 1 node
    connected=0;
    count=1;
    while connected==0
        count=count+2;
        test_t_copy=t_stat_mat;
        threshold=ordered_t_test(count,1);
        test_t_copy(test_t_copy<=threshold)=0;
        if LCC_70Nodes==1
            L=weight_conversion( test_t_copy,'lengths');
            D=distance_wei(L);
            charpath_no_k=charpath(D);
            if isinf(charpath_no_k) %connected graph
                connected=0;
            else
                connected=1;
                break;
            end
        end
    end
    threshold_least_redundancy_early=ordered_t_test(count); %the first t-statistic value that removes at least 1 node
    
    %This uses the given alpha value to see if the dysconnected network is
    %significant with a FWE p-value
    UI.alpha.ui=num2str(alpha); %the alpha value given
    UI.thresh.ui=num2str(threshold_least_redundancy_early); %the threshold found
    nbs=NBSrun(UI,'');
    pval=nbs.NBS.pval; %this is the p-value it will be a number if it significant and empty if it is non-significant
    
    %store the important values
    if isempty(pval)
        threshold_early(k,1)=0; %is the network significant Y / N = 1 / 0
        N_edges_E_tE(k,1)=nan; %how many edges in the dysconnected network
        P_val_E_tE(k,1)=nan; %what is the pvalue for the dsyconnected network
    else
        threshold_early(k,1)=1; %is the network significant Y / N = 1 / 0
        P_val_E_tE(k,1)=pval; %what is the pvalue for the dsyconnected network
        %this extracts the number of edges in the network
        sparse_con_mat=nbs.NBS.con_mat{1,1};
        con_mat=full(sparse_con_mat);
        G=graph(con_mat,'upper');
        N_edges_E_tE(k,1)=height(G.Edges); %how many edges in the dysconnected network
    end      
    delete([path '\Store\loo_Design_early' num2str(k) '.mat']) %deletes the stored design matrix when k is removed
    delete([path '\Store\loo_Evening_early' num2str(k) '.mat']) %deletes the stored FC matrix when k is removed 
    clear nbs  %clears the global value
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--Label the subj a late--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loo_Design_k=[Design((1:k-1),:); [0 1]; Design((k+1:end),:)]; %label the kth subject a late in the design matrix
    save([path '\Store\loo_Design_late' num2str(k)],'loo_Design_k') %saves the design matrix with k removed
    save([path '\Store\loo_Evening_late' num2str(k)],'Data') %saves the FC matrix with k removed
    UI.design.ui=[path '\Store\loo_Design_late' num2str(k) '.mat']; %uploads the design matrix with k removed to NBS
    UI.matrices.ui=[path '\Store\loo_Evening_late' num2str(k) '.mat']; %uploads the FC matrix with k removed to NBS
    UI.method.ui='Run NBS';
    UI.test.ui='t-test'; %type of statistical test
    UI.size.ui='Intensity'; %use sum of edge weights rather than sum of number of edges
    UI.perms.ui=num2str(N_perms); %number of permutations in FWE p value test
    UI.alpha.ui='1'; %significance level - 1 is used to ensure we collected all information actual alpha value is used below
    if contrast==1 %type of contrast used
        UI.contrast.ui='[1,-1]';
    elseif contrast==-1
        UI.contrast.ui='[-1,1]';
    end
    UI.exchange.ui='';
    UI.node_coor.ui='';
    UI.node_label.ui='';
    UI.thresh.ui='1'; % a t-statistic value is needed for nbs to run and to collect the data needed 1 is arbitrary, the percolation value is used later
    nbs=NBSrun(UI,''); %runs NBS 
    
    %find the threshold of least redundancy late
    t_stat_mat=nbs.NBS.test_stat; %obtain the t statistic matrix 
    ordered_t_test=sort(t_stat_mat(:),'descend'); %All the t-stastics calculated for subject k sorted in order 
    %this considered each t-statistic in turn and stops when it finds once
    %that makes the dsyconnected network lose 1 node
    connected=0;
    count=1;
    while connected==0
        count=count+2;
        test_t_copy=t_stat_mat;
        threshold=ordered_t_test(count,1);
        test_t_copy(test_t_copy<=threshold)=0;
        if LCC_70Nodes==1
            L=weight_conversion(test_t_copy,'lengths');
            D=distance_wei(L);
            charpath_no_k=charpath(D);
            if isinf(charpath_no_k) %connected graph
                connected=0;
            else
                connected=1;
                break;
            end
        end
    end
    threshold_least_redundancy_late=ordered_t_test(count); %the first t-statistic value that removes at least 1 node
    
    UI.alpha.ui=num2str(alpha); %the alpha value given
    UI.thresh.ui=num2str(threshold_least_redundancy_late); %the threshold found
    nbs=NBSrun(UI,'');    
    pval=nbs.NBS.pval;%this is the p-value it will be a number if it significant and empty if it is non-significant

    %store the important values
    if isempty(pval)
        threshold_late(k,1)=0; %is the network significant Y / N = 1 / 0
        N_edges_L_tL(k,1)=nan; %how many edges in the dysconnected network
        P_val_L_tL(k,1)=nan;  %what is the pvalue for the dsyconnected network
    else
        threshold_late(k,1)=1; %is the network significant Y / N = 1 / 0
        P_val_L_tL(k,1)=pval; %what is the pvalue for the dsyconnected network
        %this extracts the number of edges in the network
        sparse_con_mat=nbs.NBS.con_mat{1,1};
        con_mat=full(sparse_con_mat);
        G=graph(con_mat,'upper');
        N_edges_L_tL(k,1)=height(G.Edges);
    end
    delete([path '\Store\loo_Design_late' num2str(k) '.mat'])
    delete([path '\Store\loo_Evening_late' num2str(k) '.mat'])
    clear nbs
    
    %%%%%%%%%%%%%%%%%%%%%--Do the classifier step 1--%%%%%%%%%%%%%%%%%%%%%%    
    if (threshold_early(k,1)==0)&(threshold_late(k,1)==0)
        class_label_1(k,1)=nan;
        class_label(k,1)=nan;
        step(k,1)=1;
    elseif (threshold_early(k,1)==1)&(threshold_late(k,1)==0)
        class_label_1(k,1)=0;
        class_label(k,1)=0;
        step(k,1)=1;
    elseif (threshold_early(k,1)==0)&(threshold_late(k,1)==1)
        class_label_1(k,1)=1;
        class_label(k,1)=1;
        step(k,1)=1;
    else
        class_label_1(k,1)=99;
        class_label(k,1)=99;
    end
    
    %%%%%%%%%%%--Do the classifier step 2--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (threshold_early(k,1)==1)&(threshold_late(k,1)==1)
        
        %%%%%%%%%%%%%%Label the subject early use late threshold %%%%%%%%%%
        loo_Design_k=[Design((1:k-1),:); [1 0]; Design((k+1:end),:)]; %remove the kth subject
        %     %The save is because I couldnt get the runUI to work for just a matlab
        %     %variable only a saved matlab variable
        save([path '\Store\loo_Evening_early_tL' num2str(k)],'Data')
        save([path '\Store\loo_Design_early_tL' num2str(k)],'loo_Design_k')
        UI.design.ui=[path '\Store\loo_Design_early_tL' num2str(k) '.mat'];
        UI.matrices.ui=[path '\Store\loo_Evening_early_tL' num2str(k) '.mat'];
        UI.method.ui='Run NBS';
        UI.test.ui='t-test';
        UI.size.ui='Intensity';
        UI.perms.ui=num2str(N_perms);
        UI.alpha.ui=num2str(alpha);
        if contrast==1
            UI.contrast.ui='[1,-1]';
        elseif contrast==-1
            UI.contrast.ui='[-1,1]';
        end
        UI.exchange.ui='';
        UI.node_coor.ui='';
        UI.node_label.ui='';
        UI.thresh.ui=num2str(threshold_least_redundancy_late);        
        nbs=NBSrun(UI,'');              
        pval=nbs.NBS.pval;
        
        %store important values
        if isempty(pval)
            sig_E_tL(k,1)=0;
            N_edges_E_tL(k,1)=nan;
            P_val_E_tL(k,1)=nan;
        else
            sig_E_tL(k,1)=1;
            sparse_con_mat=nbs.NBS.con_mat{1,1};
            con_mat=full(sparse_con_mat);
            G=graph(con_mat,'upper');
            N_edges_E_tL(k,1)=height(G.Edges);
            P_val_E_tL(k,1)=pval(1,1);
        end
        delete([path '\Store\loo_Design_early_tL' num2str(k) '.mat'])
        delete([path '\Store\loo_Evening_early_tL' num2str(k) '.mat'])       
        clear nbs
        %%%%%%%%%%%%%%Label the subject late use early threshold%%%%%%%%%%
        loo_Design_k=[Design((1:k-1),:); [0 1]; Design((k+1:end),:)]; %remove the kth subject
        save([path '\Store\loo_Evening_late_tE' num2str(k)],'Data')
        save([path '\Store\loo_Design_late_tE' num2str(k)],'loo_Design_k')
        UI.design.ui=[path '\Store\loo_Design_late_tE' num2str(k) '.mat'];
        UI.matrices.ui=[path '\Store\loo_Evening_late_tE' num2str(k) '.mat'];
        UI.method.ui='Run NBS';
        UI.test.ui='t-test';
        UI.size.ui='Intensity';
        UI.perms.ui=num2str(N_perms);
        UI.alpha.ui=num2str(alpha);
        if contrast==1
            UI.contrast.ui='[1,-1]';
        elseif contrast==-1
            UI.contrast.ui='[-1,1]';
        end
        UI.exchange.ui='';
        UI.node_coor.ui='';
        UI.node_label.ui='';
        UI.thresh.ui=num2str(threshold_least_redundancy_early); %Note I tried 2.2 to start with but for every single LOOCV it went on to select 2.1 as the significiant threshold        
        nbs=NBSrun(UI,'');    
        pval=nbs.NBS.pval;
        
        %store important values
        if isempty(pval)
            sig_L_tE(k,1)=0;
            N_edges_L_tE(k,1)=nan;
            P_val_L_tE(k,1)=nan;
        else
            sig_L_tE(k,1)=1;
            sparse_con_mat=nbs.NBS.con_mat{1,1};
            con_mat=full(sparse_con_mat);
            G=graph(con_mat,'upper');
            N_edges_L_tE(k,1)=height(G.Edges);
            P_val_L_tE(k,1)=pval(1,1);
        end
        delete([path '\Store\loo_Design_late_tE' num2str(k) '.mat'])
        delete([path '\Store\loo_Evening_late_tE' num2str(k) '.mat'])
        clear nbs
        
        if (sig_E_tL(k,1)==0)&(sig_L_tE(k,1)==0)
            class_label_1(k,2)=nan;
            class_label(k,1)=nan;
            step(k,1)=2;
        elseif (sig_E_tL(k,1)==1)&(sig_L_tE(k,1)==0)
            class_label_1(k,2)=0;
            class_label(k,1)=0;
            step(k,1)=2;
        elseif (sig_E_tL(k,1)==0)&(sig_L_tE(k,1)==1)
            class_label_1(k,2)=1;
            class_label(k,1)=1;
            step(k,1)=2;
        else
            class_label_1(k,2)=99;
            class_label(k,1)=99;
        end
        
        %%%%%%%%%%%%%---Do the classifier step 3----%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (sig_E_tL(k,1)==1)&(sig_L_tE(k,1)==1)
            step(k,1)=3;
            if (N_edges_E_tE(k,1)==N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)==N_edges_L_tL(k,1))
                class_label_1(k,3)=nan;
                class_label(k,1)=nan;
            elseif (N_edges_E_tE(k,1)>N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)<N_edges_L_tL(k,1))
                class_label_1(k,3)=nan;
                class_label(k,1)=nan;
            elseif (N_edges_E_tE(k,1)<N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)>N_edges_L_tL(k,1))
                class_label_1(k,3)=nan;
                class_label(k,1)=nan;
            elseif (N_edges_E_tE(k,1)>N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)>N_edges_L_tL(k,1))
                class_label_1(k,3)=0;
                class_label(k,1)=0;
            elseif (N_edges_E_tE(k,1)<N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)<N_edges_L_tL(k,1))
                class_label_1(k,3)=1;
                class_label(k,1)=1;
            elseif (N_edges_E_tE(k,1)>N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)==N_edges_L_tL(k,1))
                class_label_1(k,3)=0.01;
                class_label(k,1)=0.01;
            elseif (N_edges_E_tE(k,1)<N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)==N_edges_L_tL(k,1))
                class_label_1(k,3)=1.01;
                class_label(k,1)=1.01;
            elseif (N_edges_E_tE(k,1)==N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)<N_edges_L_tL(k,1))
                class_label_1(k,3)=1.01;
                class_label(k,1)=1.01;
            elseif (N_edges_E_tE(k,1)==N_edges_L_tE(k,1))&(N_edges_E_tL(k,1)>N_edges_L_tL(k,1))
                class_label_1(k,3)=0.01;
                class_label(k,1)=0.01;
            end
        end
    end
end
delete([path '\Store'])   
P_val=[P_val_E_tE,P_val_L_tL,P_val_E_tL,P_val_L_tE];
N_edges=[N_edges_E_tE,N_edges_L_tL,N_edges_E_tL,N_edges_L_tE];
Accuracy=sum(floor(class_label)==Label_true)/(N_Subj)*100;
class_label=floor(class_label);
disp(['Accuracy of the classifier: ' num2str(Accuracy)])
end
