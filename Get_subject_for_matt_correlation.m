clear all; close all; clc

% Set up subjects
load('Results/Fractional_Area.mat');
load('../CCA/files/Permutation_100000_MattData.mat');
%load('/vols/Scratch/janineb/HCP/CCA/CCA_new_PFMs/Results/JBgrotU_S820_M50_Aug16_spatial_features.mat');
subs = dir('/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/Subjects');
subs = subs(3:end); subs(169) = [];
S = zeros(size(subs,1),2);
for n = 1:size(subs,1); S(n,1) = str2double(subs(n).name); end
subs_missing = [];
for n = 1:length(subs_Matt)
    N = find(S(:,1)==subs_Matt(n));
    if N; S(N,2) = 1; else subs_missing = [subs_missing; n]; end
end
subs_Matt(subs_missing) = [];

% Load CCA results
load('Results/CCA_MMP_PFM.mat','Uall')
U = mean(Uall(:,[2 6]),2);
corr([U Uall(:,[2 6])])

% Get all subject MMPs
MGall = zeros(64984,360,length(subs_Matt));
for s = 1:length(subs_Matt)
    fprintf('running subject %d \n', s)
    if exist(fullfile('Results','Data_Matt',sprintf('%d_rfMRI_REST_Atlas_MSMAll_2_d41_WRN_DeDrift_hp2000_clean_tclean_nobias_vn_BC_CorticalAreas_dil_Individual.ptseries.nii',subs_Matt(s))),'file');
        MG = ft_read_cifti(fullfile('Results','Data_Matt',sprintf('%d_rfMRI_REST_Atlas_MSMAll_2_d41_WRN_DeDrift_hp2000_clean_tclean_nobias_vn_BC_CorticalAreas_dil_Individual.ptseries.nii',subs_Matt(s))));
        MG = MG.brainordinate.parcellation;
        for n = 1:360
            MGall(MG==n,n,s) = 1;
        end
    end
end

% Correlate against CCA
MGcorr = zeros(size(MGall,1),360); MGp = zeros(size(MGall,1),360);
for n = 1:360
    [MGcorr(:,n), MGp(:,n)] = corr(squeeze(MGall(:,n,:))',U);
end
% MGcorrFWE = zeros(1000,1);
% for x = 1:1000
%     fprintf('MMP permutation %d\n',x);
%     for n = 1:360
%         m = corr(squeeze(MGall(:,n,:))',U(PAPset(:,x)));
%         MGcorrFWE(x) = nanmax(m(:));
%     end
% end
[~,I] = nanmax(abs(MGcorr),[],2);
MGcorrmax = zeros(size(MGall,1),1); MGpmax = zeros(size(MGall,1),1);
for n = 1:size(MGall,1)
    MGcorrmax(n) = MGcorr(n,I(n))*-1;
    MGpmax(n) = MGp(n,I(n));
end
[~,~,p_fdr] = fdr(MGpmax); sig = MGcorrmax(p_fdr<0.05); 
fprintf('Minimum absolute correlation for FDR significance (MMPs): %1.4f\n',min(abs(sig)));
example = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
example.dtseries = nan(size(example.dtseries,1),2); example.time = 1:2; example.hdr.dim(6) = 2;
example.dtseries(1:64984,1) = MGcorrmax;
example.dtseries(1:64984,2) = abs(MGcorrmax);
ft_write_cifti('Results/MMP_correlation_N441',example,'parameter','dtseries')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get all subjects PFMs
P = '/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel';
maps = PFM_loadSubjectSpatialMaps(P,1:50);
S = [S(1:168,:); [169 0]; S(169:end,:)];
maps(S(:,2)==0,:,:) = [];
PFMcorr = zeros(91282,50); PFMp = zeros(91282,50);
for i = 1:size(maps,3)
    [PFMcorr(:,i), PFMp(:,i)] = corr(U,maps(:,:,i));
end
% PFMcorrFWE = zeros(1000,1);
% for x = 1:1000
%     fprintf('PFM permutation %d\n',x);
%     for i = 1:size(maps,3)
%         m = corr(U(PAPset(:,x)),maps(:,:,i));
%         PFMcorrFWE(x) = nanmax(m(:));
%     end
% end
[~,I] = max(abs(PFMcorr),[],2);
PFMcorrmax = zeros(91282,1); PFMpmax = zeros(91282,1);
for n = 1:91282
    PFMcorrmax(n) = PFMcorr(n,I(n))*-1;
    PFMpmax(n) = PFMp(n,I(n));
end
[~,~,p_fdr] = fdr(PFMpmax); sig = PFMcorrmax(p_fdr<0.0001); 
fprintf('Minimum absolute correlation for FDR significance (PFMs): %1.4f\n',min(abs(sig)));
example = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
example.dtseries = repmat(example.dtseries(:,1),1,2); example.time = 1:2; example.hdr.dim(6) = 2;
example.dtseries(isnan(example.dtseries(:,1))==0,1) = PFMcorrmax;
example.dtseries(isnan(example.dtseries(:,1))==0,2) = abs(PFMcorrmax);
ft_write_cifti('Results/PFM_correlation_N441',example,'parameter','dtseries')

%fprintf('\nFWE thresholds: \nMMPs = %1.4f \nPFMs = %1.4f\n',prctile(MGcorrFWE,95),prctile(PFMcorrFWE,95));

