function [Results,class_label]=NBS_classifier_vary_t_stat(Data,Design_matrix,contrast,N_perms,alpha,path,threshold_min,threshold_max,Title,y_label,figure_index.random_seed)
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
%A function which calculates the leave one out cross validated accuracy
% INPUT:
% Data = NROIxNROIxN_Subj FC matrix
% Design_matrix = N_Subjx2 Design_matrix matrix
% contrast: 1 = the contrast [1,-1] ECP>LCP
%          -1 = the contrast [-1,1] ECP<LCP
% N_perms = the number of permuations used to calculate the FWE corrected
%           p-values e.g.1000 or 5000
% alpha = the threshold for significance below which the dysconnected
%         network is significant e.g 0.05
% path = a path given as a char i.e. '~\documents\NBS_store'
% threshold_min = the minimum t-statistic threshold you want to consider
%                 e.g. 0.01
% threshold_max = the maximum t-statistic threshold you want to consider
%                 e.g. 4.5
% random_seed= the value for the random seed 'shuffle' makes it random while a scalar i.e 1 will generate the same values each time.
% OUTPUT:
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
% class_label = N_Subjx1 vector storing the class label assigned by the
%              classifier to each of the subjects

N_Subj=size(Data,3); %size for a NROIxNROIxNSubj FC matrix
N_early=sum(Design_matrix(:,1)); %Number of ECPs
disp(['Number of Earlies: ' num2str(N_early)])
N_late=sum(Design_matrix(:,2)); %Number of LCPs
disp(['Number of Lates: ' num2str(N_late)])
Label_true=[zeros(N_early,1);ones(N_late,1)]; %true labels assigned from Design_matrix matrix
test_thresholds=(threshold_min:0.01:threshold_max);
N_thresh=size(test_thresholds,2);
class_label=ones(N_Subj,N_thresh)*99;
N_edges_E_tE=ones(N_Subj,N_thresh)*-99;
N_edges_L_tL=ones(N_Subj,N_thresh)*-99;
P_val_E_tE=ones(N_Subj,N_thresh)*-99;
P_val_L_tL=ones(N_Subj,N_thresh)*-99;
%p_val_all=ones(N_Subj,2)*-99;
step=ones(N_Subj,N_thresh)*-99;
Accuracy=ones(N_thresh,1)*-99;
Unclassified=ones(N_thresh,1)*-99;
Inaccuracy=ones(N_thresh,1)*-99;

for k=1:N_Subj
    rng(random_seed)
    disp(['Subject: ' num2str(k) ' out of ' num2str(N_Subj)])
    %Label the subject an early
    loo_Evening_k=Data; %(:,:,[1:22,24:38]);
    loo_Design_k=[Design_matrix((1:k-1),:); [1 0]; Design_matrix((k+1:end),:)]; %remove the kth subject
    %     %The save is because I couldnt get the runUI to work for just a matlab
    %     %variable only a saved matlab variable
    save([path '\Store\loo_Evening_early' num2str(k)],'loo_Evening_k')
    save([path '\Store\loo_Design_early' num2str(k)],'loo_Design_k')
    UI.design.ui=[path '\Store\loo_Design_early' num2str(k) '.mat'];
    UI.matrices.ui=[path '\Store\loo_Evening_early' num2str(k) '.mat'];
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
    
    count=0;
    for thresh=test_thresholds
        count=count+1;
        UI.thresh.ui=num2str(thresh); 
        nbs=NBSrun(UI,'');       
        %global nbs
        %nbs.NBS;
        pval=nbs.NBS.pval;
        con_matrix=nbs.NBS.con_mat;
        component=find(pval==min(pval));
        if size(component,2)>1
        clearvars N_edges_comp
            for comp_index=size(component,2)
                comp=component;
                sparse_con_mat=con_matrix{1,comp_index};
                con_mat=full(sparse_con_mat);
                G=graph(con_mat,'upper');
                N_edges_comp(1,comp_index)=height(G.Edges);
            end
            component=find(N_edges_comp==max(N_edges_comp));
            if component>1
                component=component(1,1);
            end
        end
        if isempty(pval)
            %threshold_early(k,count)=0;
            N_edges_E_tE(k,count)=0;
            P_val_E_tE(k,count)=nan;          
        else
            %threshold_early(k,count)=1;
            sparse_con_mat=nbs.NBS.con_mat{1,1};
            con_mat=full(sparse_con_mat);
            G=graph(con_mat,'upper');
            N_edges_E_tE(k,count)=height(G.Edges);
            P_val_E_tE(k,count)=pval(1,1);     
        end
    end
    delete([path '\Store\loo_Design_early' num2str(k) '.mat'])
    delete([path '\Store\loo_Evening_early' num2str(k) '.mat'])
    
    %Label the subj a late
    loo_Evening_k=Data;
    loo_Design_k=[Design_matrix((1:k-1),:); [0 1]; Design_matrix((k+1:end),:)];
    %     %The save is because I couldnt get the runUI to work for just a matlab
    %     %variable only a saved matlab variable
    save([path '\Store\loo_Evening_late' num2str(k)],'loo_Evening_k')
    save([path '\Store\loo_Design_late' num2str(k)],'loo_Design_k')
    UI.design.ui=[path '\Store\loo_Design_late' num2str(k) '.mat'];
    UI.matrices.ui=[path '\Store\loo_Evening_late' num2str(k) '.mat'];
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
    count=0;
    for thresh=test_thresholds
        count=count+1;
        UI.thresh.ui=num2str(thresh);
        nbs=NBSrun(UI,'');
        %global nbs
        %nbs.NBS
        pval=nbs.NBS.pval;
        con_matrix=nbs.NBS.con_mat;
        
        component=find(pval==min(pval));
        if size(component,2)>1
        clearvars N_edges_comp
            for comp_index=size(component,2)
                sparse_con_mat=con_matrix{1,comp_index};
                con_mat=full(sparse_con_mat);
                G=graph(con_mat,'upper');
                N_edges_comp(1,comp_index)=height(G.Edges);
            end
            component=find(N_edges_comp==max(N_edges_comp));
            if component>1
                component=component(1,1);
            end
        end
        if isempty(pval)
            %threshold_late(k,count)=0;
            N_edges_L_tL(k,count)=0;
            P_val_L_tL(k,count)=nan;
        else
            %threshold_late(k,count)=1;
            sparse_con_mat=nbs.NBS.con_mat{1,1};
            con_mat=full(sparse_con_mat);
            G=graph(con_mat,'upper');
            N_edges_L_tL(k,count)=height(G.Edges);
            P_val_L_tL(k,count)=pval(1,1);
        end
    end
    delete([path '\Store\loo_Design_late' num2str(k) '.mat'])
    delete([path '\Store\loo_Evening_late' num2str(k) '.mat'])
    %p_val_all(k,1)=P_val_E_tE(k,count);
    %p_val_all(k,2)=P_val_L_tL(k,count);
end

count=0;
for thresh=test_thresholds
    count=count+1;
    for k=1:N_Subj
        %Do the classifier step 1
        if (N_edges_E_tE(k,count)==0)&(N_edges_L_tL(k,count)==0)
            %class_label_1(k,count)=nan;
            class_label(k,count)=nan;
            step(k,count)=1;
        elseif (N_edges_E_tE(k,count)~=0)&(N_edges_L_tL(k,count)==0)
            %class_label_1(k,count)=0;
            class_label(k,count)=0;
            step(k,count)=1;
        elseif (N_edges_E_tE(k,count)==0)&(N_edges_L_tL(k,count)~=0)
            %class_label_1(k,count)=1;
            class_label(k,count)=1;
            step(k,count)=1;
        else
            %class_label_1(k,count)=99;
            class_label(k,count)=99;
        end
        %Do the classifier step 2
        if (N_edges_E_tE(k,count)~=0)&(N_edges_L_tL(k,count)~=0)
            step(k,count)=2;
            if (N_edges_E_tE(k,count)==N_edges_L_tL(k,count))
                %class_label_1(k,3)=nan;
                class_label(k,count)=nan;
            elseif (N_edges_E_tE(k,count)>N_edges_L_tL(k,count))
                %class_label_1(k,3)=0;
                class_label(k,count)=0;
            elseif (N_edges_E_tE(k,count)<N_edges_L_tL(k,count))
                %class_label_1(k,3)=1;
                class_label(k,count)=1;                
            end
        end
        class_label(k,count)
    end
Accuracy(count,1)=(sum((class_label(:,count)==Label_true))/N_Subj) *100;
Unclassified(count,1)=(sum(isnan(class_label(:,count))/N_Subj) *100);
Inaccuracy(count,1)=100-Accuracy(count,1)-Unclassified(count,1);
end
Results=table(test_thresholds',Accuracy,Unclassified,Inaccuracy);

%% Plot results Accuracy bar chart
f=figure(figure_index);clf;
f.Units='centimeters';
f.Color='white';
f.Position=[10,10,15,10];
binsize=5;
% cd('\\adf.bham.ac.uk\eps\prhome21\S\SLM585\PHD\First Year\Chronotype_research\NBS\Full_pipeline\Afternoon1_-1\Threshold_list')
% set(groot,'defaultAxesTickLabelInterpreter','latex')
% A_table=load('Table_accuracy_un_missclassified.mat');
stats=Results.Accuracy;
stats_bin=ones(3,90)*-90;
for bin=1:90
    stats_bin(1,bin)=mean(stats(1+binsize*(bin-1):binsize+binsize*(bin-1),1));
end
stats=Results.Inaccuracy;
for bin=1:90
    stats_bin(2,bin)=mean(stats(1+binsize*(bin-1):binsize+binsize*(bin-1),1));
end
stats=Results.Unclassified;
for bin=1:90
    stats_bin(3,bin)=mean(stats(1+binsize*(bin-1):binsize+binsize*(bin-1),1));
end
hold on
bar(stats_bin',1,'stacked','EdgeColor',[.2 .2 .2],'LineWidth',0.0001)
xticks([1 11 21 31 41 51 61 71 81 91])
xticklabels((0:0.5:4.5))
xlim([0.5 90])
ylim([0,100])
%xline(15.23,'--k')
xlabel('Threshold Values','interpreter', 'latex')
h1 = text(-10, 13,Y_Label,'interpreter', 'latex');
set(h1, 'rotation', 90)
title(Title,'interpreter', 'latex')
hold off
