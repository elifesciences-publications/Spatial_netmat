clear all; close all; clc

% Create separate contiguous parcels from Yeo's 17 network parcellation
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
S = '/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/files/';
surfaceL = fullfile(S,'HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii');
surfaceR = fullfile(S,'HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii');
YEO = fullfile(S,'HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/RSN-networks.32k_fs_LR.dlabel.nii');
system(sprintf('wb_command -cifti-all-labels-to-rois %s 2 %s',YEO,'Yeo_ROIs.dscalar.nii'));
system(sprintf('wb_command -cifti-find-clusters %s 0.5 20 0.5 0.5 COLUMN Yeo_ROIs_clustered.dscalar.nii -left-surface %s -right-surface %s', 'Yeo_ROIs.dscalar.nii', surfaceL, surfaceR));

Y = ft_read_cifti('Yeo_ROIs_clustered.dscalar.nii');
Yeo_all = [];
a = fieldnames(Y);
Tn = strncmp('x',a,1);
for n = find(Tn)'
    tt = a{n};
    M = Y.(tt);
    if isa(M,'double')
        u = unique(M);
        if length(u)>1
            fprintf('including %s \n', tt);
            u = setdiff(u,0);
            for m = 1:length(u)
                O = zeros(64984,1);
                O(M==u(m)) = 1;
                Yeo_all = [Yeo_all O];
            end
        end
    end
end

T = ft_read_cifti(fullfile(S,'HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/','Q1-Q6_R440.MyelinMap_BC.32k_fs_LR.dscalar.nii'));
Yeo_all(isnan(T.myelinmap_bc)==1,:) = nan;
PFMmaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
Yeo_all = [Yeo_all; zeros(size(PFMmaps.dtseries,1)-size(Yeo_all,1),size(Yeo_all,2))];
PFMmaps.hdr.dim(6) = size(Yeo_all,2); PFMmaps.time = 1:size(Yeo_all,2);
PFMmaps.dtseries = Yeo_all;
ft_write_cifti('Yeo_ROIs_final',PFMmaps,'parameter','dtseries');

