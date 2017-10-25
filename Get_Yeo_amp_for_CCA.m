clear all; close all; clc

% Set data directory
addpath ~steve/NETWORKS/FSLNets;
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/')
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
S = fullfile('/vols/Scratch/janineb/','HCP','CCA');
Nkeep = 100;
load(fullfile(S,'files','Permutation_100000.mat'),'conf');
conf = demean(conf);

% Load Yeo parcellation
Yeo = ft_read_cifti('/vols/Scratch/janineb/HCP/netmat_simulations/Yeo_ROIs_final.dtseries.nii');
Yeo = Yeo.dtseries; Yeo(isnan(Yeo(:,1))==1,:) = [];

% Get amplitudes
subs = dir('/vols/Scratch/janineb/PROFUMO/PFM50_MSMall900/Model7_S.pfm/Subjects/*');
subs = subs(3:end);
runs = {'1_LR','1_RL','2_LR','2_RL'};
Yeo_amp = zeros(length(subs),size(Yeo,2));
for s = 1:length(subs);
    msk_tmp2 = [];
    fprintf('running subject %d \n', s)
    for x = 1:4;
        % Get timeseries from original data
        D = ft_read_cifti(sprintf('/vols/Data/HCP/Phase2/subjects900/%s/MNINonLinear/Results/rfMRI_REST%s/rfMRI_REST%s_Atlas_MSMAll_hp2000_clean.dtseries.nii',subs(s).name,runs{x},runs{x}));
        D = D.dtseries; D(isnan(D(:,1))==1,:) = [];
        msk_tmp2 = [msk_tmp2; nets_demean(D' * Yeo)];
    end
    Yeo_amp(s,:) = std(msk_tmp2);
end

% Remove confounds
Yeo_amp(169,:) = [];
Yeo_amp = demean(Yeo_amp);    
Yeo_amp = demean(Yeo_amp-conf*(pinv(conf)*Yeo_amp)); 

% Run svd
[Yeo_amp,~,~]=nets_svds(Yeo_amp,Nkeep); 

% Save rsults
load('Results/input_CCA_table1.mat')
save('Results/input_CCA_table1.mat','ICA200_amp','ICA200_spatial','ICA25_amp','ICA25_spatial','PFM50_amp','PFM50_spatial','Yeo_amp')

