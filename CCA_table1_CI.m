clear all; close all; clc

% Set and load inputs:
Ddir = 'input_CCA_table1_withWarp.mat';
Dir = '/home/fs0/janineb/scratch/';
Ddir = Ddir(1:end-4);
S = fullfile(Dir,'HCP','netmat_simulations');
Nperm = 0;

% Load Steve's reference result
U_ICA = load(fullfile('/vols/Scratch/janineb/','HCP','CCA','files','ICA200_MSMall_PartialNetmat_JBgrotU.txt'));

% Set up paths
path(path,fullfile(Dir,'matlab'))
path(path,fullfile(Dir,'HCP','CCA','files'))
path(path,fullfile(Dir,'matlab','cifti-matlab'))
path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
addpath ~steve/NETWORKS/FSLNets;
path(path,fullfile(Dir,'HCP','netmat_simulations','palm-alpha106'));

brainsAll = load(fullfile(S,'Results',Ddir));
A1 = fieldnames(brainsAll);
CIs = zeros(size(A1,1),6); % R0(1) Rsur(1) CIlow CIhigh Fudgefactor min_least_squares
for x = 14%:size(A1,1)
    fprintf('running CCA on %s.%s\n',Ddir,A1{x});
    
    % Load data
    brains = brainsAll.(A1{x});
    if size(brains,1)==441
        load(fullfile('../CCA','files','Permutation_100000_MattData.mat'),'uu2');
    elseif size(brains,1)==790
        load(fullfile('../CCA','files','Permutation_task_100000.mat'),'uu2');
    else
        load(fullfile('../CCA','files','Permutation_100000.mat'),'uu2');
    end
    X_true = brains;
    Y_true = uu2;
    
    % Run main CCA
    [~,~,R0,~,~,~] = canoncorr(X_true,Y_true); CIs(x,1) = R0(1);
    
    %%% Estimate fudge factor
    Ftry = 0:0.01:0.2;
    Rtry = zeros(100,length(Ftry));
    for n = 1:length(Ftry)
        for i=1:100
            % Force subject covariances to be the same
            XY = mvnrnd(zeros(size([X_true Y_true]))',[X_true Y_true]*[X_true Y_true]')';
            XY = mvnrnd(XY,cov([X_true Y_true]));
            X1 = XY(:,1:size(X_true,2)); Y1 = XY(:,size(X_true,2)+1:end);
            
            % Different subject  covariances
            X = mvnrnd(zeros(size(X_true))',X_true*X_true')';
            Y = mvnrnd(zeros(size(Y_true))',Y_true*Y_true')';
            XY = mvnrnd([X Y],cov([X_true Y_true]));
            X2 = XY(:,1:size(X_true,2)); Y2 = XY(:,size(X_true,2)+1:end);
            
            % Average the two with fudge factor
            F = Ftry(n);
            X = mean(cat(3,F*X1,X2),3);
            Y = mean(cat(3,F*Y1,Y2),3);
            
            % Run CCA
            [~,~,R,~,~,~] = canoncorr(X,Y);
            Rtry(i,n) = R(1);
        end
    end
    %figure; plot(Ftry,(mean(Rtry)-R0(1)).^2); set(gca,'xtick',Ftry)
    [f,F] = min((mean(Rtry)-R0(1)).^2); F = Ftry(F); 
    CIs(x,5) = F; CIs(x,6) = f;
    
    %%% Obtain confidence intervals:
    
    Rsur = zeros(1000,1);
    for i=1:1000
        fprintf('surrogate %d out of %d\n',i,1000);
        % Force subject covariances to be the same
        XY = mvnrnd(zeros(size([X_true Y_true]))',[X_true Y_true]*[X_true Y_true]')';
        XY = mvnrnd(XY,cov([X_true Y_true]));
        X1 = XY(:,1:size(X_true,2)); Y1 = XY(:,size(X_true,2)+1:end);
        
        % Different subject  covariances
        X = mvnrnd(zeros(size(X_true))',X_true*X_true')';
        Y = mvnrnd(zeros(size(Y_true))',Y_true*Y_true')';
        XY = mvnrnd([X Y],cov([X_true Y_true]));
        X2 = XY(:,1:size(X_true,2)); Y2 = XY(:,size(X_true,2)+1:end);
        
        % Average the two with fudge factor
        X = mean(cat(3,F*X1,X2),3);
        Y = mean(cat(3,F*Y1,Y2),3);
        
        % Run CCA
        [~,~,R,~,~,~] = canoncorr(X,Y);
        Rsur(i) = R(1);
    end
    CIs(x,2) = mean(Rsur);
    CIs(x,3:4) = prctile(Rsur,[2.5 97.5]);
     
end

% save
save(fullfile('Results','results_CCA_table1_CIs.mat'),'CIs','A1');




