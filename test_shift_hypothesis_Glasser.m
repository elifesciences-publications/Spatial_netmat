clear all; close all; clc

% Set inputs:
group_amplitude = 1; % Set to 1 to switch off using subject-based amplitude stuff
group_maps = 0; % Set to 1 to switch off using subject PFM maps in simulation
map_thresh = -95; % Set to 1 to threshold subject maps to remove netmat-like map differences
map_bin = 1; % Set to 1 to threshold and binarize subject maps
subject_nets = 0; % Set to 1 to use subject netmats rather than group netmats
Group_parcels = 0; % Set to 1 to use group Glasser parcellation rather than subject
rho = 0.23; % estimated using Glasser_optimise_rho.m 

% Set paths:
path(path,fullfile('/vols/Scratch/janineb','matlab','cifti-matlab'))
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
addpath ~steve/NETWORKS/FSLNets;
path(path,'/home/fs0/janineb/scratch/matlab')

% Load PFM group maps
PFMmaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
PFMmaps = PFMmaps.dtseries; dt_isnan = isnan(PFMmaps(:,1)); PFMmaps(isnan(PFMmaps(:,1))==1,:) = [];
if group_amplitude == 1
    load('group_amplitude_inputs.mat')
end

% Load ground truth netmats created using test_shift_hypothesis_create_ground_truth.m
load Ground_truth_netmat_PFMnew.mat

% Load PFM subject maps
PFMdir = '/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/';
Subject_maps = PFM_loadSubjectSpatialMaps(PFMdir,1:50);

% Load PFM subject timeseries
ts = PFM_loadTimeCourses_pfmNEW(PFMdir,0.72,1,1,0,[]);

% Create simulated data and perform weighted regression of ICA maps onto this
subs_all = [];
drTC_NEW = [];
drTC_ORIG = [];
missing_nodes = zeros(500,10);
subs = dir('/vols/Scratch/janineb/PROFUMO/PFM50_MSMall900/Model7_S.pfm/Subjects/*');
subs = subs(3:end);
runs = {'1_LR','1_RL','2_LR','2_RL'};

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
        TCall = ts.ts((s-1)*ts.NtimepointsPerSubject+1:s*ts.NtimepointsPerSubject,:);
        msk_tmp = [];
        msk_tmp2 = [];
        
        for x = 1:4;
            % Get run timeseries
            TC = TCall((x-1)*1200+1:x*1200,:);
            
            % Get subject correlations
            sCov = cov(TC);
            sStd = sqrt(diag(sCov));
            sCorr = diag(1./sStd) * sCov * diag(1./sStd);
            
            % And make a new set of time series without any temporal information
            newTC = TC * diag( 1 ./ sStd );         % Variance normalise
            newTC = newTC * sCorr ^-0.5;            % Remove correlations
            if subject_nets == 1
                newTC = newTC * sCorr^0.5;
            else
                newTC = newTC * GTnet^0.5;          % Add in ground truth correlations
            end
            newTC = newTC * diag( 1./std(newTC) );  % Variance normalise again
            if group_amplitude == 1
                sStd = sStd_group;
            end
            newTC = newTC * diag( sStd );           % Add back in the original variances
            
            % Load posterior noise variance and amplitude weights needed for creating full dataset
            if group_amplitude == 0
                a = importfile(fullfile(PFMdir,'Subjects', subs(s).name,'Runs',runs{x},'NoisePrecision.post','GammaPosterior.txt'));
                noiseStd = sqrt(a(2)/a(1));
                H = h5read(fullfile(PFMdir,'Subjects',subs(s).name,'Runs',runs{x},'ComponentWeightings.post','Means.hdf5'),'/dataset');
            elseif group_amplitude == 1
                noiseStd = noiseStd_group;
                H = H_group;
            end
            
            % Select correct maps and apply threshold if desired
            if group_maps == 0
                M = squeeze(Subject_maps(s,:,:));
            elseif group_maps == 1
                M = PFMmaps;
            end
            if map_thresh > 0
                Mnew = zeros(size(M));
                Mnew(M<-map_thresh) = M(M<-map_thresh);
                Mnew(M>map_thresh) = M(M>map_thresh);
                M = Mnew;
            end
            if map_thresh < 0
                Mnew = zeros(size(M));
                for p = 1:50
                    Uthr(s,p) = prctile(M(:,p),-map_thresh);
                    Lthr(s,p) = prctile(M(:,p),100+map_thresh);
                    Mnew(M(:,p)<Lthr(s,p),p) = M(M(:,p)<Lthr(s,p),p);
                    Mnew(M(:,p)>Uthr(s,p),p) = M(M(:,p)>Uthr(s,p),p);
                end
                M = Mnew;
            end
            if map_bin ~= 0
                Mnew = zeros(size(M));
                Mnew(M<0) = -1;
                Mnew(M>0) = 1;
                M = Mnew;
            end
            
            % Take outer product and add noise to create full dataset
            D = M * diag(H) * newTC' + noiseStd * randn(91282, 1200);
            
            % Undo variance normalisation that PFM pipeline applies (to avoid issues with relative strength of subcortical modes)
            if group_amplitude == 0
                norm = h5read(fullfile(PFMdir,'..','Preprocessing',subs(s).name,sprintf('%s_Normalisation.hdf5',runs{x})),'/dataset')';
            elseif group_amplitude == 1
                norm = norm_group';
            end
            D = bsxfun(@times, D, 1./norm');
            
            % Run masking against Yeo parcellation
            msk_tmp = [msk_tmp; nets_demean(D' * MG)];
            
            % Get timeseries from original data
            D = ft_read_cifti(sprintf('/vols/Data/HCP/Phase2/subjects900/%s/MNINonLinear/Results/rfMRI_REST%s/rfMRI_REST%s_Atlas_MSMAll_hp2000_clean.dtseries.nii',subs(s).name,runs{x},runs{x}));
            D = D.dtseries; D(isnan(D(:,1))==1,:) = [];
            msk_tmp2 = [msk_tmp2; nets_demean(D' * MG)];
        end
        drTC_NEW = [drTC_NEW; msk_tmp];
        drTC_ORIG = [drTC_ORIG; msk_tmp2];
    end
end
clear D i dr_tmp dd s TC newTC Snet Sstd
        
% Initialise output
OUTPUT_all_nodes = zeros(2,5);
OUTPUT_cortex = zeros(2,5);
OUTPUT_subcortex = zeros(2,5);
% Row 1 = full and Row 2 = partial

% Calculate original and new subject netmats
tsNEW = ts; tsNEW.Nnodes = 360; tsNEW.NnodesOrig = 360; tsNEW.DD = 1:360;
tsNEW.Nsubjects = length(subs_all); 
tsNEW.Ntimepoints = size(drTC_NEW,1);
tsNEW.ts = drTC_NEW; clear drTC_NEW
tsORIG = tsNEW; 
tsORIG.ts = drTC_ORIG; clear drTC_ORIG
F.netNEW = nets_netmats(tsNEW,-1,'corr');
F.netORIG = nets_netmats(tsORIG,-1,'corr');

% Calculate partial netmats and replace missing edges with group mean at
% covariance stage:
P.netORIG = zeros(size(F.netORIG)); P.netNEW = zeros(size(F.netNEW));
output_precision=0; N = tsORIG.Nnodes; do_rtoz = -1;
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

% Calculate partial netmats and replace missing edges with group mean at
% covariance stage:
P.netNEW = zeros(size(F.netNEW)); 
output_precision=0; N = tsNEW.Nnodes; do_rtoz = -1;
cov_group_mean = zeros(tsNEW.Nsubjects,tsNEW.Nnodes*tsNEW.Nnodes);
for s = 1:tsNEW.Nsubjects
    grot = tsNEW.ts((s-1)*tsNEW.NtimepointsPerSubject+1:s*tsNEW.NtimepointsPerSubject,:);
    grot = cov(grot);  
    grot(grot==0) = nan;
    cov_group_mean(s,:) = reshape(grot,1,N*N);
end
cov_group_mean = nanmean(cov_group_mean);
for s = 1:tsNEW.Nsubjects
    grot = tsNEW.ts((s-1)*tsNEW.NtimepointsPerSubject+1:s*tsNEW.NtimepointsPerSubject,:);
    grot = cov(grot);  
    scale_factor = pinv(cov_group_mean(:))*grot(:);
    grot(grot==0) = scale_factor.*cov_group_mean(grot==0);
    grot = grot/sqrt(mean(diag(grot).^2));
    grot = nearestSPD(grot);
    grot = inv(grot+rho*eye(N));
    if output_precision==0
        grot=-grot; grot=(grot ./ repmat(sqrt(abs(diag(grot))),1,N)) ./ repmat(sqrt(abs(diag(grot)))',N,1);  grot(eye(N)>0)=0;
    end
    P.netNEW(s,:) = reshape(grot,1,N*N);
end
P.netNEW = 0.5*log((1+P.netNEW)./(1-P.netNEW))* (-do_rtoz);

% Replace missing edges with group mean
NanEdges = isnan(F.netORIG);
A = nanmean(P.netNEW); A = repmat(A,tsNEW.Nsubjects,1); P.netNEW(isnan(F.netNEW)==1) = A(isnan(F.netNEW)==1);
A = nanmean(F.netNEW); A = repmat(A,tsNEW.Nsubjects,1); F.netNEW(isnan(F.netNEW)==1) = A(isnan(F.netNEW)==1);
A = nanmean(P.netORIG); A = repmat(A,tsORIG.Nsubjects,1); P.netORIG(isnan(F.netORIG)==1) = A(isnan(F.netORIG)==1);
A = nanmean(F.netORIG); A = repmat(A,tsORIG.Nsubjects,1); F.netORIG(isnan(F.netORIG)==1) = A(isnan(F.netORIG)==1);
clear A

% Normalise netmats
F.netNEW_norm = F.netNEW ./ repmat(std(F.netNEW,[],2),1,size(F.netNEW,2));
F.netORIG_norm = F.netORIG ./ repmat(std(F.netORIG,[],2),1,size(F.netORIG,2));
P.netNEW_norm = P.netNEW ./ repmat(std(P.netNEW,[],2),1,size(P.netNEW,2));
P.netORIG_norm = P.netORIG ./ repmat(std(P.netORIG,[],2),1,size(P.netORIG,2));

% Calculate group average netmats
[F.ZnetNEW,F.MnetNEW] = nets_groupmean(F.netNEW,0);
[F.ZnetORIG,F.MnetORIG] = nets_groupmean(F.netORIG,0);
[P.ZnetNEW,P.MnetNEW] = nets_groupmean(P.netNEW,0);
[P.ZnetORIG,P.MnetORIG] = nets_groupmean(P.netORIG,0);
[F.ZnetNEW_norm,F.MnetNEW_norm] = nets_groupmean(F.netNEW_norm,0);
[F.ZnetORIG_norm,F.MnetORIG_norm] = nets_groupmean(F.netORIG_norm,0);
[P.ZnetNEW_norm,P.MnetNEW_norm] = nets_groupmean(P.netNEW_norm,0);
[P.ZnetORIG_norm,P.MnetORIG_norm] = nets_groupmean(P.netORIG_norm,0);

% Correlate subject*subject correlation matrices
F.corr_NEW = corr(F.netNEW'); F.corr_ORIG = corr(F.netORIG');
A = F.corr_NEW(:); A(eye(size(F.corr_NEW))==1) = []; B = F.corr_ORIG(:); B(eye(size(F.corr_ORIG))==1) = [];
OUTPUT_all_nodes(1,1) = corr(A,B);
P.corr_NEW = corr(P.netNEW'); P.corr_ORIG = corr(P.netORIG');
A = P.corr_NEW(:); A(eye(size(F.corr_NEW))==1) = []; B = P.corr_ORIG(:); B(eye(size(F.corr_NEW))==1) = [];
OUTPUT_all_nodes(2,1) = corr(A,B);

% Subtract & regress netmats
GRP = repmat(F.MnetNEW(:)',tsNEW.Nsubjects,1);
F.netNEWsubtr = F.netNEW - GRP;
F.netNEWregr = F.netNEW - demean(GRP) * pinv(demean(GRP))*demean(F.netNEW);
GRP = repmat(F.MnetORIG(:)',tsNEW.Nsubjects,1);
F.netORIGsubtr = F.netORIG - GRP;
F.netORIGregr = F.netORIG - demean(GRP) * pinv(demean(GRP))*demean(F.netORIG);
GRP = repmat(P.MnetNEW(:)',tsNEW.Nsubjects,1);
P.netNEWsubtr = P.netNEW - GRP;
P.netNEWregr = P.netNEW - demean(GRP) * pinv(demean(GRP))*demean(P.netNEW);
GRP = repmat(P.MnetORIG(:)',tsNEW.Nsubjects,1);
P.netORIGsubtr = P.netORIG - GRP;
P.netORIGregr = P.netORIG - demean(GRP) * pinv(demean(GRP))*demean(P.netORIG);
clear GRP
GRP = repmat(F.MnetNEW_norm(:)',tsNEW.Nsubjects,1);
F.netNEWsubtr_norm = F.netNEW_norm - GRP;
F.netNEWregr_norm = F.netNEW_norm - demean(GRP) * pinv(demean(GRP))*demean(F.netNEW_norm);
GRP = repmat(F.MnetORIG_norm(:)',tsNEW.Nsubjects,1);
F.netORIGsubtr_norm = F.netORIG_norm - GRP;
F.netORIGregr_norm = F.netORIG_norm - demean(GRP) * pinv(demean(GRP))*demean(F.netORIG_norm);
GRP = repmat(P.MnetNEW_norm(:)',tsNEW.Nsubjects,1);
P.netNEWsubtr_norm = P.netNEW_norm - GRP;
P.netNEWregr_norm = P.netNEW_norm - demean(GRP) * pinv(demean(GRP))*demean(P.netNEW_norm);
GRP = repmat(P.MnetORIG_norm(:)',tsNEW.Nsubjects,1);
P.netORIGsubtr_norm = P.netORIG_norm - GRP;
P.netORIGregr_norm = P.netORIG_norm - demean(GRP) * pinv(demean(GRP))*demean(P.netORIG_norm);
clear GRP

% Compare difference between individual and group netmats
[OUTPUT_all_nodes(1,2)] = compare_nets(F.netNEWsubtr,F.netORIGsubtr);
[OUTPUT_all_nodes(1,3)] = compare_nets(F.netNEWregr,F.netORIGregr);
[OUTPUT_all_nodes(1,4)] = compare_nets(F.netNEWsubtr_norm,F.netORIGsubtr_norm);
[OUTPUT_all_nodes(1,5)] = compare_nets(F.netNEWregr_norm,F.netORIGregr_norm);
[OUTPUT_all_nodes(2,2)] = compare_nets(P.netNEWsubtr,P.netORIGsubtr);
[OUTPUT_all_nodes(2,3)] = compare_nets(P.netNEWregr,P.netORIGregr);
[OUTPUT_all_nodes(2,4)] = compare_nets(P.netNEWsubtr_norm,P.netORIGsubtr_norm);
[OUTPUT_all_nodes(2,5)] = compare_nets(P.netNEWregr_norm,P.netORIGregr_norm);

% Concatenate output
OUTPUT = [OUTPUT_all_nodes OUTPUT_cortex OUTPUT_subcortex];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Look at CCA results for original and new netmats

S = '/home/fs0/janineb/scratch/HCP/CCA/';
Nkeep = 100;
Nperm = 100000;
% [uu2,names,varskeep,varsd,PAPset,conf] = vars2uu2(Nperm,setdiff(1:820,subs_all));
% save(fullfile(S,'files','Permutation_100000_MattData'),'uu2','names','varskeep','varsd','PAPset','conf');
load(fullfile(S,'files','Permutation_100000_MattData.mat'),'conf'); conf = demean(conf);

FccaNEW = demean(F.netNEW); FccaNEW(subs_all==169,:) = []; NanEdges(subs_all==169,:) = [];
FccaNEW = demean(FccaNEW-conf*(pinv(conf)*FccaNEW)); A = FccaNEW;
[FccaNEW,~,~]=nets_svds(FccaNEW,Nkeep);
A(NanEdges==1) = nan;
COV=zeros(size(A,1));
for i=1:size(A,1)
  for j=1:size(A,1)
    grot=A([i j],:); grot=cov(grot(:,sum(isnan(grot))==0)'); COV(i,j)=grot(1,2);
  end
end
COV2=nearestSPD(COV);
[FccaNEW_imputed,~]=eigs(COV2,Nkeep);

FccaORIG = demean(F.netORIG); FccaORIG(subs_all==169,:) = [];
FccaORIG = demean(FccaORIG-conf*(pinv(conf)*FccaORIG)); A = FccaORIG;
[FccaORIG,~,~]=nets_svds(FccaORIG,Nkeep);
A(NanEdges==1) = nan;
COV=zeros(size(A,1));
for i=1:size(A,1)
  for j=1:size(A,1)
    grot=A([i j],:); grot=cov(grot(:,sum(isnan(grot))==0)'); COV(i,j)=grot(1,2);
  end
end
COV2=nearestSPD(COV);
[FccaORIG_imputed,~]=eigs(COV2,Nkeep);

PccaNEW = demean(P.netNEW); PccaNEW(subs_all==169,:) = [];
PccaNEW = demean(PccaNEW-conf*(pinv(conf)*PccaNEW)); A = PccaNEW;
[PccaNEW,~,~]=nets_svds(PccaNEW,Nkeep);
A(NanEdges==1) = nan;
COV=zeros(size(A,1));
for i=1:size(A,1)
  for j=1:size(A,1)
    grot=A([i j],:); grot=cov(grot(:,sum(isnan(grot))==0)'); COV(i,j)=grot(1,2);
  end
end
COV2=nearestSPD(COV);
[PccaNEW_imputed,~]=eigs(COV2,Nkeep);

PccaORIG = demean(P.netORIG); PccaORIG(subs_all==169,:) = [];
PccaORIG = demean(PccaORIG-conf*(pinv(conf)*PccaORIG)); A = PccaORIG;
[PccaORIG,~,~]=nets_svds(PccaORIG,Nkeep);
A(NanEdges==1) = nan;
COV=zeros(size(A,1));
for i=1:size(A,1)
  for j=1:size(A,1)
    grot=A([i j],:); grot=cov(grot(:,sum(isnan(grot))==0)'); COV(i,j)=grot(1,2);
  end
end
COV2=nearestSPD(COV);
[PccaORIG_imputed,~]=eigs(COV2,Nkeep);

A = '';
if group_maps == 1; A = [A '_using_group_maps']; end
if group_amplitude ==1; A = [A '_using_group_amps']; end
if map_thresh ~= 0; A = sprintf('%s_map_thresh_%d',A,map_thresh); end
if map_bin ~= 0; A = sprintf('%s_map_thresh_%d_bin',A,map_bin); end
if subject_nets == 1; A = sprintf('%s_using_subject_netmats',A); end
if Group_parcels == 1; A = sprintf('_group_parcellation_%s',A); end
save(sprintf('Results/results_Glasser%s.mat',A),'OUTPUT','F','P','subs_all','missing_nodes','-v7.3')
save(sprintf('Results/input_CCA_Glasser%s.mat',A),'FccaNEW_imputed','FccaORIG_imputed','PccaNEW_imputed','PccaORIG_imputed')   
save(sprintf('Results/input_CCA_Glasser_notimputed%s.mat',A),'FccaNEW','FccaORIG','PccaNEW','PccaORIG')   




