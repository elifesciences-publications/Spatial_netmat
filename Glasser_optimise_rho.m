clear all; close all; clc

Group_parcels = 0; % Set to 1 to use group Glasser parcellation rather than subject

% Set paths:
path(path,fullfile('/vols/Scratch/janineb','matlab','cifti-matlab'))
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
addpath ~steve/NETWORKS/FSLNets;
path(path,'/home/fs0/janineb/scratch/matlab')
PFMdir = '/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/';

% Create simulated data and perform weighted regression of ICA maps onto this
subs_all = [];
drTC_ORIG = [];
missing_nodes = zeros(819,15);
subs = dir('/vols/Scratch/janineb/PROFUMO/PFM50_MSMall900/Model7_S.pfm/Subjects/*');
subs = subs(3:end);
runs = {'1_LR','1_RL','2_LR','2_RL'};
ts = PFM_loadTimeCourses_pfmNEW(PFMdir,0.72,1,1,0,[]);
PFMmaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
PFMmaps = PFMmaps.dtseries; dt_isnan = isnan(PFMmaps(:,1)); clear PFMmaps

for s = 1:ts.Nsubjects
    fprintf('running subject %d \n', s)
    
    % Get subject parcellation
    if exist(fullfile('Results','Data_Matt',sprintf('%s_rfMRI_REST_Atlas_MSMAll_2_d41_WRN_DeDrift_hp2000_clean_tclean_nobias_vn_BC_CorticalAreas_dil_Individual.ptseries.nii',subs(s).name)),'file');
        subs_all = [subs_all; s];
        if Group_parcels == 1;
            MG = ft_read_cifti('/vols/Data/HCP/workbench/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group.32k_fs_LR.dlabel.nii');
            MG = MG.indexmax;
        else
            MG = ft_read_cifti(fullfile('Results','Data_Matt',sprintf('%s_rfMRI_REST_Atlas_MSMAll_2_d41_WRN_DeDrift_hp2000_clean_tclean_nobias_vn_BC_CorticalAreas_dil_Individual.ptseries.nii',subs(s).name)));
            MG = MG.brainordinate.parcellation;
        end
        MG(dt_isnan(1:size(MG,1))==1) = [];
        I = setdiff(1:360,unique(MG));
        if I; missing_nodes(s,1:length(I)) = I; end
        
        % Make into ROIs
        MGnew = zeros(91282,360);
        for n = 1:360
            MGnew(MG==n,n) = 1;
        end
        MG = MGnew;
        
        % Get subject timeseries
        msk_tmp2 = [];
        for x = 1:4;        
            % Get timeseries from original data
            D = ft_read_cifti(sprintf('/vols/Data/HCP/Phase2/subjects900/%s/MNINonLinear/Results/rfMRI_REST%s/rfMRI_REST%s_Atlas_MSMAll_hp2000_clean.dtseries.nii',subs(s).name,runs{x},runs{x}));
            D = D.dtseries; D(isnan(D(:,1))==1,:) = [];
            msk_tmp2 = [msk_tmp2; nets_demean(D' * MG)];
        end
        drTC_ORIG = [drTC_ORIG; msk_tmp2];
    end
end
clear D i dr_tmp dd s TC newTC Snet Sstd

% Calculate original and new subject netmats
tsORIG = ts; tsORIG.Nnodes = 360; tsORIG.NnodesOrig = 360; tsORIG.DD = 1:360;
tsORIG.Nsubjects = length(subs_all); 
tsORIG.Ntimepoints = size(drTC_ORIG,1);
tsORIG.ts = drTC_ORIG; clear drTC_ORIG

% Calculate partial netmats and replace missing edges with group mean at
% covariance stage:
P.netORIG = zeros(length(subs_all),tsORIG.Nnodes*tsORIG.Nnodes);
rho = 0.01; output_precision=0; N = tsORIG.Nnodes; do_rtoz = -1;
cov_group_mean = zeros(tsORIG.Nsubjects,tsORIG.Nnodes*tsORIG.Nnodes);
for s = 1:tsORIG.Nsubjects
    grot = tsORIG.ts((s-1)*tsORIG.NtimepointsPerSubject+1:s*tsORIG.NtimepointsPerSubject,:);
    grot = cov(grot);  
    grot(grot==0) = nan;
    cov_group_mean(s,:) = reshape(grot,1,N*N);
end
cov_group_mean = nanmean(cov_group_mean);
for s = 1:tsORIG.Nsubjects
    grot = tsORIG.ts((s-1)*tsORIG.NtimepointsPerSubject+1:s*tsORIG.NtimepointsPerSubject,:);
    grot = cov(grot);  
    scale_factor = pinv(cov_group_mean(:))*grot(:);
    grot(grot==0) = scale_factor.*cov_group_mean(grot==0);
    grot = grot/sqrt(mean(diag(grot).^2));
    grot = nearestSPD(grot);
    grot = inv(grot+rho*eye(N));
    if output_precision==0
        grot=-grot; grot=(grot ./ repmat(sqrt(abs(diag(grot))),1,N)) ./ repmat(sqrt(abs(diag(grot)))',N,1);  grot(eye(N)>0)=0;
    end
    P.netORIG(s,:) = reshape(grot,1,N*N);
end
P.netORIG = 0.5*log((1+P.netORIG)./(1-P.netORIG))* (-do_rtoz);
[~,A] = nets_groupmean(P.netORIG,0);
A = reshape(A,1,N*N);

% Loop over range of rho's to find optimum
rho_range = 0.01:0.02:0.5;
rho_corrs = zeros(length(subs_all),length(rho_range));
for rhoN = 1:length(rho_range)
    P.netORIG = zeros(length(subs_all),tsORIG.Nnodes*tsORIG.Nnodes);
    rho = rho_range(rhoN); output_precision=0; N = tsORIG.Nnodes; do_rtoz = -1;
    fprintf('doing rho %d (%1.2f) out of %d\n',rhoN,rho,length(rho_range));
    cov_group_mean = zeros(tsORIG.Nsubjects,tsORIG.Nnodes*tsORIG.Nnodes);
    for s = 1:tsORIG.Nsubjects
        grot = tsORIG.ts((s-1)*tsORIG.NtimepointsPerSubject+1:s*tsORIG.NtimepointsPerSubject,:);
        grot = cov(grot);
        grot(grot==0) = nan;
        cov_group_mean(s,:) = reshape(grot,1,N*N);
    end
    cov_group_mean = nanmean(cov_group_mean);
    for s = 1:tsORIG.Nsubjects
        grot = tsORIG.ts((s-1)*tsORIG.NtimepointsPerSubject+1:s*tsORIG.NtimepointsPerSubject,:);
        grot = cov(grot);
        scale_factor = pinv(cov_group_mean(:))*grot(:);
        grot(grot==0) = scale_factor.*cov_group_mean(grot==0);
        grot = grot/sqrt(mean(diag(grot).^2));
        grot = nearestSPD(grot);
        grot = inv(grot+rho*eye(N));
        if output_precision==0
            grot=-grot; grot=(grot ./ repmat(sqrt(abs(diag(grot))),1,N)) ./ repmat(sqrt(abs(diag(grot)))',N,1);  grot(eye(N)>0)=0;
        end
        P.netORIG(s,:) = reshape(grot,1,N*N);
    end
    P.netORIG = 0.5*log((1+P.netORIG)./(1-P.netORIG))* (-do_rtoz);
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
title(sprintf('optimising rho for Glasser partial correlation; max %1.2f, rho %1.2f',n,rho_range(i)))
set(gcf,'Position',[500 500 1000 500],'PaperPositionMode','auto')
print(gcf,'-dtiff','-r300','Results/Glasser_rho_optimisation.tif')


