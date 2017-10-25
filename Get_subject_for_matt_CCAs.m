clear all; close all; clc

% Set up paths
Dir = '/home/fs0/janineb/scratch/';
path(path,fullfile(Dir,'matlab'))
path(path,fullfile(Dir,'HCP','CCA','files'))
path(path,fullfile(Dir,'matlab','cifti-matlab'))
path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
addpath ~steve/NETWORKS/FSLNets;
Nperm = 100000;
%Nperm = 0;
S = fullfile(Dir,'HCP','CCA');

% Load PFM spatial
load('Results/input_CCA_table1.mat','PFM50_spatial');

% Load Glasser spatial
load('Results/input_CCA_table1.mat','Glasser_maps');

% Load Glasser netmat
load('Results/input_CCA_Glasser.mat','FccaORIG_imputed','PccaORIG_imputed');
Glasser_Fnet = FccaORIG_imputed; Glasser_Pnet = PccaORIG_imputed;
clear FccaORIG_imputed PccaORIG_imputed

% Remove subjects not in Glasser
load('Results/Data_Matt/ts_real_sim_subs.mat','subs_all');
PFM50_spatial = [PFM50_spatial(1:168,:); zeros(1,size(PFM50_spatial,2)); PFM50_spatial(169:end,:)];
subs_remove = setdiff(1:820,subs_all); 
PFM50_spatial(unique([169 subs_remove]),:) = [];

% Load canonical CCA result against behaviour
JBgrotU_ICA = load(fullfile(S,'files','ICA200_MSMall_PartialNetmat_JBgrotU.txt'));
JBgrotU_ICA = [JBgrotU_ICA(1:168,:); ones(1,size(JBgrotU_ICA,2)); JBgrotU_ICA(169:end,:)];
load('Results/Data_Matt/ts_real_sim_subs.mat','subs_all');
subs_remove = setdiff(1:820,subs_all);
JBgrotU_ICA(unique([subs_remove 169]),:) = [];

%%%%%%%%%%%%%%%%%%% CCA PFM maps against Glasser maps %%%%%%%%%%%%%%%%%%%
% Run CCA
[mapsA,mapsB,mapsR,mapsU,mapsV,mapsS] = canoncorr(PFM50_spatial,Glasser_maps);
R = corr([JBgrotU_ICA mapsU]); R = abs(R(2:end,1));
[mapsRica, mapsIica] = max(R);
fprintf('\nPFM maps against Glasser maps \nR=%1.2f \n Ruu=%1.2f \nRuu_i=%1.0f', mapsR(1),mapsRica,mapsIica);

% Run permutations
Nkeep = 100;
if Nperm ~= 0
    load(fullfile(S,'files',sprintf('Permutation_%d_MattData.mat',Nperm)));
    mapsR(101) = mean(mapsR);
    grotRp=zeros(Nperm,Nkeep+1); clear grotRpval;
    for j=1:Nperm
        [grotAr,grotBr,grotRp(j,1:end-1),grotUr,grotVr,grotstatsr]=canoncorr(PFM50_spatial,Glasser_maps(PAPset(:,j),:)); grotRp(j,end)=mean(grotRp(j,1:end-1));
    end
    for I=1:Nkeep;
        mapsP(I)=(1+sum(grotRp(2:end,1)>=mapsR(I)))/Nperm;
    end
    %%%%% Print output
    fprintf('\np=%1.6f \nnumber significant=%d', mapsP(1), sum(mapsP<0.05));
end

%%%%%%%%%%%%%%%%%%% CCA PFM maps against Glasser full netmats %%%%%%%%%%%%%%%%%%%
% Run CCA
[fnetA,fnetB,fnetR,fnetU,fnetV,fnetS] = canoncorr(PFM50_spatial,Glasser_Fnet);
R = corr([JBgrotU_ICA mapsU]); R = abs(R(2:end,1));
[fnetRica, fnetIica] = max(R);
fprintf('\nPFM maps against Glasser maps \nR=%1.2f \n Ruu=%1.2f \nRuu_i=%1.0f', fnetR(1),fnetRica,fnetIica);

% Run permutations
Nkeep = 100;
if Nperm ~= 0
    load(fullfile(S,'files',sprintf('Permutation_%d_MattData.mat',Nperm)));
    fnetR(101) = mean(fnetR);
    grotRp=zeros(Nperm,Nkeep+1); clear grotRpval;
    for j=1:Nperm
        [grotAr,grotBr,grotRp(j,1:end-1),grotUr,grotVr,grotstatsr]=canoncorr(PFM50_spatial,Glasser_Fnet(PAPset(:,j),:)); grotRp(j,end)=mean(grotRp(j,1:end-1));
    end
    for I=1:Nkeep;
        fnetP(I)=(1+sum(grotRp(2:end,1)>=fnetR(I)))/Nperm;
    end
    %%%%% Print output
    fprintf('\np=%1.6f \nnumber significant=%d', fnetP(1), sum(fnetP<0.05));
end

%%%%%%%%%%%%%%%%%%% CCA PFM maps against Glasser partial netmats %%%%%%%%%%%%%%%%%%%
% Run CCA
[pnetA,pnetB,pnetR,pnetU,pnetV,pnetS] = canoncorr(PFM50_spatial,Glasser_Pnet);
R = corr([JBgrotU_ICA mapsU]); R = abs(R(2:end,1));
[pnetRica, pnetIica] = max(R);
fprintf('\nPFM maps against Glasser maps \nR=%1.2f \n Ruu=%1.2f \nRuu_i=%1.0f', pnetR(1),pnetRica,pnetIica);

% Run permutations
Nkeep = 100;
if Nperm ~= 0
    load(fullfile(S,'files',sprintf('Permutation_%d_MattData.mat',Nperm)));
    pnetR(101) = mean(pnetR);
    grotRp=zeros(Nperm,Nkeep+1); clear grotRpval;
    for j=1:Nperm
        [grotAr,grotBr,grotRp(j,1:end-1),grotUr,grotVr,grotstatsr]=canoncorr(PFM50_spatial,Glasser_Pnet(PAPset(:,j),:)); grotRp(j,end)=mean(grotRp(j,1:end-1));
    end
    for I=1:Nkeep;
        pnetP(I)=(1+sum(grotRp(2:end,1)>=pnetR(I)))/Nperm;
    end
    %%%%% Print output
    fprintf('\np=%1.6f \nnumber significant=%d', pnetP(1), sum(pnetP<0.05));
end

save('Get_subjects_for_matt_CCA_results.mat',...
    'pnetA','pnetB','pnetR','pnetU','pnetV','pnetS','pnetRica','pnetIica','pnetP',...
    'fnetA','fnetB','fnetR','fnetU','fnetV','fnetS','fnetRica','fnetIica','fnetP',...
    'mapsA','mapsB','mapsR','mapsU','mapsV','mapsS','mapsRica','mapsIica','mapsP',...
    '-v7.3');

