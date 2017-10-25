clear all; close all; clc

Group_parcels = 0; % Set to 1 to use group Glasser parcellation rather than subject

% Set paths:
path(path,fullfile('/vols/Scratch/janineb','matlab','cifti-matlab'))
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
addpath ~steve/NETWORKS/FSLNets;
path(path,'/home/fs0/janineb/scratch/matlab')
PFMdir = '/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/';

% Load Yeo parcellation
Yeo = ft_read_cifti('/vols/Scratch/janineb/HCP/netmat_simulations/Yeo_ROIs_final.dtseries.nii');
Yeo = Yeo.dtseries; Yeo(isnan(Yeo(:,1))==1,:) = [];

% Create simulated data and perform weighted regression of ICA maps onto this
drTC_ORIG = [];
subs = dir('/vols/Scratch/janineb/PROFUMO/PFM50_MSMall900/Model7_S.pfm/Subjects/*');
subs = subs(3:end);
runs = {'1_LR','1_RL','2_LR','2_RL'};
ts = PFM_loadTimeCourses_pfmNEW(PFMdir,0.72,1,1,0,[]);
PFMmaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
PFMmaps = PFMmaps.dtseries; dt_isnan = isnan(PFMmaps(:,1)); clear PFMmaps

for s = 1:ts.Nsubjects
    fprintf('running subject %d \n', s)
    % Get subject timeseries
    msk_tmp2 = [];
    for x = 1:4;
        % Get timeseries from original data
        D = ft_read_cifti(sprintf('/vols/Data/HCP/Phase2/subjects900/%s/MNINonLinear/Results/rfMRI_REST%s/rfMRI_REST%s_Atlas_MSMAll_hp2000_clean.dtseries.nii',subs(s).name,runs{x},runs{x}));
        D = D.dtseries; D(isnan(D(:,1))==1,:) = [];
        msk_tmp2 = [msk_tmp2; nets_demean(D' * Yeo)];
    end
    drTC_ORIG = [drTC_ORIG; msk_tmp2];
end
clear D i dr_tmp dd s TC newTC Snet Sstd

% Calculate original and new subject netmats
tsORIG = ts; tsORIG.Nnodes = 108; tsORIG.NnodesOrig = 108; tsORIG.DD = 1:108;
tsORIG.ts = drTC_ORIG; clear drTC_ORIG

% Calculate partial netmats 
P.netORIG = nets_netmats(tsORIG,-1,'ridgep',0.01);
[~,A] = nets_groupmean(P.netORIG,0);
A = reshape(A,1,tsORIG.Nnodes*tsORIG.Nnodes);

% Loop over range of rho's to find optimum
rho_range = 0.01:0.02:0.5;
rho_corrs = zeros(tsORIG.Nsubjects,length(rho_range));
for rhoN = 1:length(rho_range)
    rho = rho_range(rhoN); 
    P.netORIG = nets_netmats(tsORIG,-1,'ridgep',rho);
    for s = 1:tsORIG.Nsubjects
        rho_corrs(s,rhoN) = corr(P.netORIG(s,:)',A');
    end        
end

% R-to-z rho corrs & plot results
rho_corrs = 0.5*log((1+rho_corrs)./(1-rho_corrs));
figure; plot(rho_range,mean(rho_corrs));
xlabel('rho')
ylabel('mean z-transformed r with average netmat')
set(gca,'xtick',rho_range)
[n,i] = max(mean(rho_corrs));
title(sprintf('optimising rho for Yeo partial correlation; max %1.2f, rho %1.2f',n,rho_range(i)))
set(gcf,'Position',[500 500 1000 500],'PaperPositionMode','auto')
print(gcf,'-dtiff','-r300','Results/Yeo_rho_optimisation.tif')


