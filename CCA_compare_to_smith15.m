clear all; close all; clc

% Set up paths
Dir = '/home/fs0/janineb/scratch/';
path(path,fullfile(Dir,'matlab'))
path(path,fullfile(Dir,'HCP','CCA','files'))
path(path,fullfile(Dir,'matlab','cifti-matlab'))
path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
addpath ~steve/NETWORKS/FSLNets;

% Run CCA to get subject weights
Nperm = 0;
load('Results/input_CCA_200.mat');
brains = PccaORIG;
[~,CCA.JBgrotU,CCA.JBgrotV,CCA.JBgrotR,CCA.Rmax,CCA.Imax] = S_CCA(brains,Nperm);

% Load Steve's subject weights
load('~steve/newCCA6b.mat','grotU');
id500 = load('/vols/Data/HCP/Phase2/scripts500/vars/vars.txt'); id500 = id500(:,1);
id900 = load('/vols/Data/HCP/Phase2/scripts900/vars/vars.txt'); id900 = id900(:,1); 
I = find(id500)==id900(169); id500(I) = []; grotU(I,:) = []; id900(169) = [];
I = zeros(length(id500),1);
toremove = [];
for n = 1:length(id500); 
    A = find(id900==id500(n));
    if A
        I(n) = A; 
    else
        fprintf('removed subject %d\n',id500(n))
        toremove = [toremove; n];
    end
    clear A
end
id500(toremove) = [];
grotU(toremove,:) = [];
I(toremove) = [];
        
% Compare
[R,p] = corr(grotU(:,1),CCA.JBgrotU(I,1))

