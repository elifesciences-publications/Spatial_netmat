clear all; close all; clc

load('Results/Data_Matt/ts_real_sim.mat','drTC_ORIG');
load('Results/Data_Matt/ts_real_sim_subs.mat');
S = '/home/fs0/janineb/scratch/HCP/CCA/';
Nkeep = 100;
Nperm = 100000;
% [uu2,names,varskeep,varsd,PAPset,conf] = vars2uu2(Nperm,setdiff(1:820,subs_all));
% save(fullfile(S,'files','Permutation_100000_MattData'),'uu2','names','varskeep','varsd','PAPset','conf');
load(fullfile(S,'files','Permutation_100000_MattData.mat'),'conf'); 
conf = demean(conf);
Pconf = pinv(conf);

% CCA inputs netmats
ts.Nnodes = 360; ts.NnodesOrig = 360; ts.DD = 1:360;
ts.Nsubjects = length(subs_all); 
ts.Ntimepoints = size(drTC_ORIG,1);
ts.NtimepointsPerSubject = 4800;
ts.ts = drTC_ORIG; clear drTC_ORIG
F = nets_netmats(ts,-1,'corr');
A = nanmean(F); A = repmat(A,ts.Nsubjects,1);
F(isnan(F)==1) = A(isnan(F)==1);
F(subs_all==169,:) = [];
F = demean(F);    
F = demean(F-conf*(pinv(conf)*F));
[Glasser_Fnetmat,~,~]=nets_svds(F,Nkeep); 

% CCA inputs amplitudes
[ts_stats,all_stats] = nets_stats(ts);
temp_amp = all_stats.std;
temp_amp(subs_all==169,:) = [];
temp_amp = demean(temp_amp);    
temp_amp = demean(temp_amp-conf*(pinv(conf)*temp_amp));
[Glasser_amp,~,~]=nets_svds(temp_amp,Nkeep);
clear ts ts_stats all_stats F A temp_amp 

% CCA inputs maps
PFMmaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
PFMmaps = PFMmaps.dtseries; dt_isnan = isnan(PFMmaps(:,1)); clear PFMmaps
subs = dir('/vols/Scratch/janineb/PROFUMO/PFM50_MSMall900/Model7_S.pfm/Subjects/*');
subs = subs(3:end);
maps = zeros(442,91282,360); Isub = 1;
for s = 1:820
    fprintf('running subject %d \n', s)
    
    % Get subject parcellation
    if exist(fullfile('Results','Data_Matt',sprintf('%s_rfMRI_REST_Atlas_MSMAll_2_d41_WRN_DeDrift_hp2000_clean_tclean_nobias_vn_BC_CorticalAreas_dil_Individual.ptseries.nii',subs(s).name)),'file');
        MG = ft_read_cifti(fullfile('Results','Data_Matt',sprintf('%s_rfMRI_REST_Atlas_MSMAll_2_d41_WRN_DeDrift_hp2000_clean_tclean_nobias_vn_BC_CorticalAreas_dil_Individual.ptseries.nii',subs(s).name)));
        MG = MG.brainordinate.parcellation;
        MG(dt_isnan(1:size(MG,1))==1) = [];
        I = setdiff(1:360,unique(MG));
        MGnew = zeros(91282,360);
        for n = 1:360
            MGnew(MG==n,n) = 1;
        end
        maps(Isub,:,:) = MGnew; 
        Isub = Isub+1;
    end
end
maps(subs_all==169,:,:) = []; 
C = zeros(size(maps,1),size(maps,1));
for i = 1:size(maps,3)
    fprintf('processing map %d (out of %d)\n',i,size(maps,3));
    m = maps(:,:,i)';
    m = demean(m,2);
    m = (m'-conf*(Pconf*m'))';    
    C = C + m' * m;
end
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
Glasser_maps = V(:,1:Nkeep);

% Save results
clear C conf D dt_isnan i I inds Isub m maps MG MGnew missing_nodes ...
    n Nkeep Nperm Pconf s S subs subs_all V
load('Results/input_CCA_table1.mat')

