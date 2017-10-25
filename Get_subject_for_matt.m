clear all; close all; clc

% Set up paths
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/')
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
path(path,'/vols/Scratch/janineb/HCP/CCA/files/')
Nsubs = 25;

% Load PFM spatial CCA results & subjects
load('/vols/Scratch/janineb/HCP/CCA/CCA_new_PFMs/Results/JBgrotU_S820_M50_Aug16_spatial_features.mat');
subs = dir('/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/Subjects');
subs = subs(3:end); subs(169) = [];
S = zeros(size(subs,1),2);
for n = 1:size(subs,1); S(n,1) = str2double(subs(n).name); end

% Order subjects according to U
load(fullfile('/home/fs0/janineb/scratch/HCP/CCA/','files','Permutation_100000.mat'),'varsd','names');
[U,I] = sort(JBgrotU(:,1));
subs = subs(I); S = S(I,:);
V = JBgrotV(I,1);
varsd = varsd(I,:);

% See what subject are in Matt's parcellation
for n = 1:size(S,1)
    if exist(fullfile('Results','Data_Matt',sprintf('%6.0f_rfMRI_REST_Atlas_MSMAll_2_d41_WRN_DeDrift_hp2000_clean_tclean_nobias_vn_BC_CorticalAreas_dil_Individual.ptseries.nii',S(n,1))),'file');
        S(n,2) = 1;
    end
end

% See what direction the correlation with behaviour is
[beh_corrs,beh_names] = CCA_get_behavioural_mode(JBgrotV,1);

% Check behaviour
THC = strcmp('THC',names); THC = find(THC==1);
PicVocab = strcmp('PicVocab_Unadj',names); PicVocab = find(PicVocab==1);
THCout = [mean(varsd(1:Nsubs,THC)) mean(varsd(Nsubs+1:end,THC))]
PicVocabout = [mean(varsd(1:Nsubs,PicVocab)) mean(varsd(Nsubs+1:end,PicVocab))]

% Save results
dlmwrite('Results/Subjects_matt_pos_to_neg.csv',S,'delimiter',',','precision',7)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Sanity Check

% Average maps across same 25 subjects matt used
Sall = S;
varsd(S(:,2)==0,:) = [];
S(S(:,2)==0,:) = [];
S = [S(1:Nsubs,1); S(end-Nsubs+1:end,1)];
varsd = [varsd(1:Nsubs,:); varsd(end-Nsubs+1:end,:)];
Ddir = 'S820_M50_Aug16.pfm';
d = 50;
PFMdir = '/home/fs0/samh/scratch/HCP/';
P = fullfile(PFMdir,Ddir,'FinalModel'); Ddir = Ddir(1:end-4);
for i = 1:size(S,1)
    means = h5read(fullfile(P,sprintf('Subjects/%6.0f',S(i)),'SpatialMaps.post','Gaussian','Means.hdf5'),'/dataset');
    probs = h5read(fullfile(P,sprintf('Subjects/%6.0f',S(i)),'SpatialMaps.post','Gaussian','MembershipProbabilities.hdf5'),'/dataset'); 
    maps = means .* probs;
    if i == 1
        subject_maps = zeros(size(S,1),size(maps,1),size(maps,2));
    end
    subject_maps(i,:,:) = maps;
end
Top = squeeze(mean(subject_maps(1:Nsubs,:,:)));
Bottom = squeeze(mean(subject_maps(Nsubs+1:end,:,:)));
GroupMaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
GroupMaps = GroupMaps.dtseries; GroupMaps(isnan(GroupMaps(:,1))==1,:) = [];

% Save out as cifti's
DMN = ft_read_cifti(fullfile('/home/fs0/janineb/scratch','HCP','DMN','DMN_mask','DMNconn_firstOrder_thres.dtseries.nii'));
DMN.dtseries = repmat(DMN.dtseries(:,1),1,150); DMN.time = 1:150; DMN.hdr.dim(6) = 150;
c = 1;
for n = 1:50
    DMN.dtseries(isnan(DMN.dtseries(:,1))==0,c) = GroupMaps(:,n); c=c+1;
    DMN.dtseries(isnan(DMN.dtseries(:,1))==0,c) = Top(:,n); c=c+1;
    DMN.dtseries(isnan(DMN.dtseries(:,1))==0,c) = Bottom(:,n); c=c+1;
end
ft_write_cifti('Results/Subjects_matt_PFMs',DMN,'parameter','dtseries')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Average correlations for my stuff
CorrJ = ft_read_cifti('/vols/Scratch/janineb/HCP/CCA/CCA_new_PFMs/Results/Maps/CCAmap_S820_M50_Aug16.dtseries.nii');
CorrJ = CorrJ.dtseries; CorrJ(isnan(CorrJ(:,1))==1,:) = [];
CorrJ = max(abs(CorrJ'));

% Average correlations for Matt's stuff:
D = ft_read_cifti('Results/Diff.dscalar.nii');
A = fieldnames(D);
I = strfind(A,'roi'); I = cellfun(@isempty,I);
DD = zeros(64984,length(find(I==0))); c = 1;
for n = 1:length(A)
    if I(n)==0
        a = D.(A{n}); a(isnan(a)) = 0;
        DD(:,c) = a; c = c+1;
    end
end
Diff = sum(abs(DD'));

% FTB across my stuff
[~,TopFTB] = max(abs(Top'));
[~,BottomFTB] = max(abs(Bottom'));

% FTB across Matt's stuff
D = ft_read_cifti('Results/Bottom.dscalar.nii');
A = fieldnames(D);
I = strfind(A,'roi'); I = cellfun(@isempty,I);
DD = zeros(64984,length(find(I==0))); c = 1;
for n = 1:length(A)
    if I(n)==0
        a = D.(A{n}); a(isnan(a)) = 0;
        DD(:,c) = a; c = c+1;
    end
end
[~,BottomFTB_M] = max(abs(DD'));
D = ft_read_cifti('Results/Top.dscalar.nii');
A = fieldnames(D);
I = strfind(A,'roi'); I = cellfun(@isempty,I);
DD = zeros(64984,length(find(I==0))); c = 1;
for n = 1:length(A)
    if I(n)==0
        a = D.(A{n}); a(isnan(a)) = 0;
        DD(:,c) = a; c = c+1;
    end
end
[~,TopFTB_M] = max(abs(DD'));

% Save out as cifti's
DMN = ft_read_cifti(fullfile('/home/fs0/janineb/scratch','HCP','DMN','DMN_mask','DMNconn_firstOrder_thres.dtseries.nii'));
DMN.dtseries = repmat(DMN.dtseries(:,1),1,6); DMN.time = 1:6; DMN.hdr.dim(6) = 6;
DMN.dtseries(isnan(DMN.dtseries(:,1))==0,1) = CorrJ;
DMN.dtseries(1:size(DD,1),2) = Diff'; DMN.dtseries(size(DD,1)+1:end,2) = zeros(size(DMN.dtseries,1)-size(DD,1),1);
DMN.dtseries(isnan(DMN.dtseries(:,1))==0,3:4) = [TopFTB' BottomFTB'];
DMN.dtseries(1:size(DD,1),5:6) = [TopFTB_M' BottomFTB_M']; DMN.dtseries(size(DD,1)+1:end,5:6) = zeros(size(DMN.dtseries,1)-size(DD,1),2);
ft_write_cifti('Results/Subjects_matt_summary',DMN,'parameter','dtseries')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get netmats
rho = 0.23;
drTC = [];
missing_nodes = zeros(500,10);
runs = {'1_LR','1_RL','2_LR','2_RL'};
PFMmaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
PFMmaps = PFMmaps.dtseries; dt_isnan = isnan(PFMmaps(:,1)); PFMmaps(isnan(PFMmaps(:,1))==1,:) = [];
clear PFMmaps
for s = 1:size(S,1)
    MG = ft_read_cifti(fullfile('Results','Data_Matt',sprintf('%d_rfMRI_REST_Atlas_MSMAll_2_d41_WRN_DeDrift_hp2000_clean_tclean_nobias_vn_BC_CorticalAreas_dil_Individual.ptseries.nii',S(s))));
    MG = MG.brainordinate.parcellation; MG(dt_isnan(1:size(MG,1))==1) = [];
    I = setdiff(1:360,unique(MG));
    if I; missing_nodes(s,1:length(I)) = I; end
    MGnew = zeros(91282,360);
    for n = 1:360; MGnew(MG==n,n) = 1; end; MG = MGnew;
    msk_tmp2 = [];
    for x = 1:4;
        D = ft_read_cifti(sprintf('/vols/Data/HCP/Phase2/subjects900/%d/MNINonLinear/Results/rfMRI_REST%s/rfMRI_REST%s_Atlas_MSMAll_hp2000_clean.dtseries.nii',S(s),runs{x},runs{x}));
        D = D.dtseries; D(isnan(D(:,1))==1,:) = [];
        msk_tmp2 = [msk_tmp2; nets_demean(D' * MG)];
    end
    drTC = [drTC; msk_tmp2];
end

ts.ts = drTC;
ts.Nnodes = 360; ts.NnodesOrig = 360; ts.DD = 1:360;
ts.Nsubjects = length(S); 
ts.Ntimepoints = size(drTC,1);
ts.NtimepointsPerSubject = 4800;

full = nets_netmats(ts,-1,'corr');

partial = zeros(size(full)); 
output_precision=0; N = ts.Nnodes; do_rtoz = -1;
cov_group_mean = zeros(ts.Nsubjects,ts.Nnodes*ts.Nnodes);
for s = 1:ts.Nsubjects
    grot = ts.ts((s-1)*ts.NtimepointsPerSubject+1:s*ts.NtimepointsPerSubject,:);
    grot = cov(grot);  
    grot(grot==0) = nan;
    cov_group_mean(s,:) = reshape(grot,1,N*N);
end
cov_group_mean = nanmean(cov_group_mean);
for s = 1:ts.Nsubjects
    grot = ts.ts((s-1)*ts.NtimepointsPerSubject+1:s*ts.NtimepointsPerSubject,:);
    grot = cov(grot);  
    scale_factor = pinv(cov_group_mean(:))*grot(:);
    grot(grot==0) = scale_factor.*cov_group_mean(grot==0);
    grot = grot/sqrt(mean(diag(grot).^2));
    grot = nearestSPD(grot);
    grot = inv(grot+rho*eye(N));
    if output_precision==0
        grot=-grot; grot=(grot ./ repmat(sqrt(abs(diag(grot))),1,N)) ./ repmat(sqrt(abs(diag(grot)))',N,1);  grot(eye(N)>0)=0;
    end
    partial(s,:) = reshape(grot,1,N*N);
end
partial = 0.5*log((1+partial)./(1-partial))* (-do_rtoz);

% Replace missing edges with group mean
A = nanmean(partial); A = repmat(A,ts.Nsubjects,1); partial(isnan(full)==1) = A(isnan(full)==1);
A = nanmean(full); A = repmat(A,ts.Nsubjects,1); full(isnan(full)==1) = A(isnan(full)==1);
clear A

top_partial = partial(1:Nsubs,:);
bottom_partial = partial(Nsubs+1:end,:);
top_full = full(1:Nsubs,:);
bottom_full = full(Nsubs+1:end,:);
save('Results/Subjects_matt_netmats.mat','top_partial','bottom_partial','top_full','bottom_full');

[Zfull_top,Mfull_top] = nets_groupmean(top_full,0);
[Zpartial_top,Mpartial_top] = nets_groupmean(top_partial,0);
[Zfull_bottom,Mfull_bottom] = nets_groupmean(bottom_full,0);
[Zpartial_bottom,Mpartial_bottom] = nets_groupmean(bottom_partial,0);

figure;
subplot(2,2,1); imagesc(Zfull_top,[-9 9]); title('full top 25');
subplot(2,2,2); imagesc(Zfull_bottom,[-9 9]); title('full bottom 25');
subplot(2,2,3); imagesc(Zpartial_top,[-9 9]); title('partial top 25');
subplot(2,2,4); imagesc(Zpartial_bottom,[-9 9]); title('partial bottom 25');

figure;
subplot(2,2,1); imagesc(Mfull_top,[-1.2 1.2]); title('full top 25');
subplot(2,2,2); imagesc(Mfull_bottom,[-1.2 1.2]); title('full bottom 25');
subplot(2,2,3); imagesc(Mpartial_top,[-0.05 0.05]); title('partial top 25');
subplot(2,2,4); imagesc(Mpartial_bottom,[-0.05 0.05]); title('partial bottom 25');
