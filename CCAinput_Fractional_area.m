clear all; close all; clc

Nkeep = 100;

% Set paths
S = '/home/fs0/janineb/scratch/HCP/CCA/CCA_new_PFMs'; 
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/')
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
addpath ~steve/NETWORKS/FSLNets;
load(fullfile(S,'..','files','Permutation_100000_MattData.mat'),'conf');
load('Results/Fractional_Area.mat');

% Get correct subjects
subs = dir('/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/Subjects');
subs = subs(3:end); subs(169) = [];
S = zeros(size(subs,1),2);
for n = 1:size(subs,1); S(n,1) = str2double(subs(n).name); end
subs_missing = [];
for n = 1:length(subs_Matt)
    N = find(S(:,1)==subs_Matt(n));
    if N; S(N,2) = 1; else subs_missing = [subs_missing; n]; end
end
Fractional_area_L(subs_missing,:) = [];
Fractional_area_R(subs_missing,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = [Fractional_area_L Fractional_area_R];

% Remove confounds
F = demean(F);    
F = demean(F-conf*(pinv(conf)*F));  
% Save input data for CCA
[Fractional_area,~,~]=nets_svds(F,Nkeep); 

save('Results/input_CCA_table1_fractional_area.mat','Fractional_area','-v7.3')


