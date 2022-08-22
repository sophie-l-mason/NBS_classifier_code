% Brunno Machado de Campos   &  Sophie Mason
% University of Campinas        University of Birmingham 
% 
% Copyright (c) 2020,
% Brunno Machado de Campos and Sophie Mason
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
%
% The TPCM final variable will be a 3D matrix where the 3th dimension are
% the subjects (Nseeds X Nseeds X Nsubjects)

plotB = 1; % 1 or 0 to plot or not graphs

% Adding the time series of all volunteers:
[file,path] = uigetfile({'*.mat','Matlab files'},...
    'Select all the Seeds_TimeSeries.mat files',...
    'MultiSelect','on');

% Obtaining the Data Matrix "X" of all subjects
for Subj = 1:size(file,2)
    Xtmp = load([path,filesep,file{Subj}]);
    X(:,:,Subj) = Xtmp.MatrixSeeds;
    % Demeaning
    X_mean = mean(X(:,:,Subj),1);
    X_mean_matrix = repmat(X_mean,size(Xtmp.MatrixSeeds,1),1);
    % Normalising by the STD
    standard_dev = std(X(:,:,Subj),1);
    standard_dev_matrix = repmat(standard_dev,size(Xtmp.MatrixSeeds,1),1);
    X(:,:,Subj) = (X(:,:,Subj) - X_mean_matrix)./(standard_dev_matrix);   
end

if plotB
    figure;imagesc(X(:,:,1)); title('Subject 1 Data Matrix:'); xlabel('ROIs (n)');
    ylabel('Time Points (tp)'); ax1 = gca; 
end

% Defining the matrix sizes:
% 70, n is the number of ROIs
n = size(X,2);

% 180, tp is the number of timepoints
tp = size(X,1);
% The  for subject 1 would be a matrix of replicated vectors. 
% The vector elements are the temporal average of each ROI.

file=zeros(size(X,3),size(X,3));

% mean at the first dimension, the time points (tp) dimension
meanX = mean(X(:,:,1),1); 

% The Subject 1 final mean matrix (), based on the replicated vector (VEC):
meanX = repmat(meanX,tp,1);


if plotB
    figure;imagesc(meanX);
    title('Subject 1 Replicated Mean Vector (matrix Xmean)');
    xlabel('ROIs (n)'); ylabel('Time Points (tp)');
end

% Calculating the Empirical Covariance matrix of all volunteers, using the above define  (meanX)
E = zeros(n,n,size(file,2)); % Prelocation;
for Subj = 1:size(file,2)
    meanX = repmat(mean(X(:,:,Subj),1),tp,1);
    E(:,:,Subj) = (1/tp) .* ((X(:,:,Subj) - meanX)' * (X(:,:,Subj) - meanX));
end

if plotB
    figure;imagesc(E(:,:,1)); title('Subject 1 Empirical Covariance Matrix');
    xlabel('ROIs (n)'); ylabel('ROIs (n)');
end

% Pure identity matrix "I" [nXn]
I = eye(n); 

% Estimating alpha (group information). Note that I am using the "invChol_mex"
% function for Cholesky Decomposition:

% mean (subjects) Tikhonov precision matrix considering alpha = 0;
for i = 1:size(E,3)
    W_a0(:,:,i) = invChol_mex(E(:,:,i)) + 0.*I; 
end
Wmean_a0 = mean(W_a0,3);

Gamma_sub = zeros(n,n,size(file,2));

% Considering 4 floating points.
h1 = waitbar(0,'Estimating optimal alpha value');
Alph = 0.0001:0.0001:1;
Attempts = 1;
for Loop = Alph
    waitbar(Attempts/numel(Alph),h1)
    for Subj = 1:size(file,2)
        % Individual Tikhonov precision matrix considering alpha = Alph minus
        % the mean (subjects) Tikhonov precision matrix considering alpha = 0;
        PMtx = E(:,:,Subj) + Loop.*I;
        Gamma_sub(:,:,Subj) = invChol_mex(PMtx) - Wmean_a0;
    end
    % squaring the subtracted matrices and summing (subject wisely)
    %     Gamma(:,:) = (sum(Gamma_sub,3)).^2; 
    Gamma(:,:) = sum((Gamma_sub.^2),3); 
    
    % getting the upper triangular matrix (i<j) --> Exclude the Diagonal
    GammaUpTri = triu(Gamma(:,:),1);
    
    % summing the triangular matrix elements
    SGammaUpTri = sum(sum(GammaUpTri)); 
    
    % square-root of the final value for the specific Loop alpha.
    Temp_Alphas(Attempts) = sqrt(SGammaUpTri); 
    Attempts = Attempts + 1;
end
close(h1)

% Finding the alpha which resulted in the lowest value:
Opt_Alpha = Alph(find(min(Temp_Alphas)==Temp_Alphas));
Opt_Alpha_iter(1,iter)=Opt_Alpha;

if isequal(Opt_Alpha,1)
    fprintf('The optimum alpha is the maximum of the test range [%.5f,%.5f] and may not be correct, the range may need to be increased\n',Alph(1),Alph(end));
    newlimits = inputdlg({sprintf('The optimum alpha is the maximum of the test range [%.5f,%.5f] and may not be correct, the range may need to be increased. Add a new limit. E.g.: 0 100\n',Alph(1),Alph(end))},...
        'Attention');
    newlimits = str2num(cell2mat(newlimits));
    Gamma_sub = zeros(n,n,size(file,2));
    clear Gamma Temp_Alphas GammaUpTri SGammaUpTri
    h1 = waitbar(0,'Estimating optimal alpha value');
    Alph = newlimits(1):0.0001:newlimits(2);
    Attempts = 1;
    for Loop = Alph
        waitbar(Attempts/numel(Alph),h1)
        for Subj = 1:size(file,2)
            PMtx = E(:,:,Subj) + Loop.*I;
            Gamma_sub(:,:,Subj) = invChol_mex(PMtx) - Wmean_a0;
        end
        Gamma(:,:) = sum((Gamma_sub.^2),3); 
        GammaUpTri = triu(Gamma(:,:),1);
        SGammaUpTri = sum(sum(GammaUpTri)); 
        Temp_Alphas(Attempts) = sqrt(SGammaUpTri); 
        Attempts = Attempts + 1;
    end
    close(h1)
    Opt_Alpha = Alph(find(min(Temp_Alphas)==Temp_Alphas));
end

[l,c]=size(Alph);
test=zeros(c-1,1);
for a=1:c-1
    if Temp_Alphas(a)>Temp_Alphas(a+1)
        test(a,1) = 1;
    end
end

if test==ones(c-1,1)
    fprintf('monotonically decreasing! – the optimal choice of alpha may not be correct, check the data is normalised!\n')
else
    fprintf('NOT monotonically decreasing :)\n')
end

if plotB
    figure;plot(Temp_Alphas); title('Function outputs'); xlabel('Attempt');
    ylabel('Resultant value');
end

if plotB
    txt = ['\uparrowargmin(alpha) = ',num2str(Opt_Alpha)]; 
    text(find(min(Temp_Alphas)==Temp_Alphas),Temp_Alphas(find(min(Temp_Alphas)==Temp_Alphas))-10,txt);  
    xlim([-0.1*numel(Alph) 1.1*numel(Alph)])
    ylim([0.90*min(Temp_Alphas),1.1*max(Temp_Alphas)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotB
    fs=10;
fn='times';
wd=15;
ht=20;
figure;plot(Temp_Alphas); 
    hx=xlabel({ '$\alpha$'}, 'interpreter','latex');
    hy=ylabel({ '$f(\alpha)$'}, 'interpreter','latex');
     xticks([1,2000,4000,6000,8000,10000]);
     xticklabels({'0', '0.2', '0.4', '0.6' ,'0.8','1'});
     xlim([-0.01*numel(Alph) 1.01*numel(Alph)])
    set(hx,'fontsize',fs); set(hx,'fontname',fn);
    set(hy,'fontsize',fs); set(hy,'fontname',fn);
    set(gca,'fontsize',fs); set(gca,'fontname',fn);
    box on; set(gca,'tickdir','out');

end
%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimating the Tikhonov partial correlation matrix for each included volunteer:
h2 = waitbar(0,'Estimating Tikhonov partial correlation for all subjects');
for Subj = 1:size(file,2)
    waitbar(Subj/size(file,2),h2)
    % Tikhonov covariance matrix [nXn]
    Ea(:,:,Subj) = E(:,:,Subj) + Opt_Alpha.*I;  

    % Tikhonov precision matrix [nXn]
    W(:,:,Subj) = invChol_mex(Ea(:,:,Subj)); 
    for pi = 1:n
        for pj = 1:n
            if ~isequal(pi,pj)
                % Tikhonov partial correlation matrix
                TPCM(pi,pj,Subj) = (-1.*W(pi,pj,Subj))/sqrt(W(pi,pi,Subj)*W(pj,pj,Subj));
            else
                TPCM(pi,pj,Subj) = NaN;
            end
        end
    end
end
close(h2)

if plotB
    figure;imagesc(TPCM(:,:,1)); title('Subject 1 Tikhonov partial correlation matrix'); 
    xlabel('ROIs (n)'); ylabel('ROIs (n)');
    ax4 = gca;
    chart2 = ax4.Children(1);datatip(chart2,27,23);
end

AllSubj3D = TPCM;