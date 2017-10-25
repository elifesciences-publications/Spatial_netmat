clear all; close all; clc

% Load full datasets:
brainsAll = load('Results/input_CCA_table1.mat');
X_true_all = brainsAll.PFM50_spatial; clear brainsAll;
uu = load('../CCA/files/Permutation_100000.mat','uu2');
Y_true_all = uu.uu2; clear uu;

% Create reduced datasets:
load('Results/Data_Matt/ts_real_sim_subs.mat','subs_all');
subs_remove = setdiff(1:820,subs_all);
X_true = [X_true_all(1:168,:); ones(1,size(X_true_all,2)); X_true_all(169:end,:)];
X_true(unique([subs_remove 169]),:) = [];
Y_true = [Y_true_all(1:168,:); ones(1,size(Y_true_all,2)); Y_true_all(169:end,:)];
Y_true(unique([subs_remove 169]),:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run true CCA on reduced numbers dataset:

[~,~,R0,~,~,~] = canoncorr(X_true,Y_true); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimise fudge factor for surrogate data

Ftry = 0:0.01:0.15;
Rtry = zeros(100,length(Ftry));
for n = 1:length(Ftry)
    fprintf('Outer loop %d out of %d\n',n,length(Ftry))
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
figure; plot(Ftry,(mean(Rtry)-R0(1)).^2); set(gca,'xtick',Ftry)
[~,F] = min((mean(Rtry)-R0(1)).^2); F = Ftry(F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain confidence intervals:

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
mean(Rsur)
prctile(Rsur,[2.5 97.5])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check against bootstrap without replacement:

Rboot = zeros(1000,1);
for i = 1:1000
    fprintf('bootstrap %d out of %d\n',i,1000);
    s = randsample(1:819,441,false);
    X = X_true_all(s,:); Y = Y_true_all(s,:);
    [~,~,R,~,~,~] = canoncorr(X,Y);
    Rboot(i) = R(1);
end
mean(Rboot)
prctile(Rboot,[2.5 97.5])
R0(1)

figure; 
hist(Rsur,50); h = findobj(gca,'Type','patch'); set(h,'FaceColor','r','EdgeColor','r','facealpha',0.5);
hold on
hist(Rboot,50); h = findobj(gca,'Type','patch'); set(h,'FaceColor','b','EdgeColor','b','facealpha',0.5);




