clear all; close all; clc

% Set up paths
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/')
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
path(path,'/vols/Scratch/janineb/HCP/CCA/files/')
path(path,'/vols/Scratch/janineb/matlab/')
Nsubs = 25;
load('Results/Fractional_Area.mat');
load(fullfile('/home/fs0/janineb/scratch/HCP/CCA/','files','Permutation_100000.mat'),'varsd','names');

% Fix ROI names
Names = strrep(Names,'_ROI','');
Names_bilateral = Names; for n = 1:length(Names); N = Names{n}; Names_bilateral{n} = N(3:end); end
Names360 = [Names; Names]; for n = 1:length(Names); N = Names{n}; Names360{n} = ['L_' N(3:end)]; end

% Load PFM spatial CCA results & subjects
load('/vols/Scratch/janineb/HCP/CCA/CCA_new_PFMs/Results/JBgrotU_S820_M50_Aug16_spatial_features.mat');
subs = dir('/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/Subjects');
subs = subs(3:end); subs(169) = [];
S = zeros(size(subs,1),2);
for n = 1:size(subs,1); S(n,1) = str2double(subs(n).name); end

% Get gender
FMRIB='/vols/Data/HCP/Phase2/scripts900'; 
Steve_names = importdata(fullfile(FMRIB,'vars','column_headers.txt'));
Steve_names(1) = {'ID'}; 
vars = load(fullfile(FMRIB,'vars','vars.txt'));
sex = nan(length(subs_Matt),1);
for n = 1:length(subs_Matt);
    N = find(vars(:,1)==subs_Matt(n));
    if N
        sex(n) = vars(N,3); 
    end
end

% Get correct subjets
subs_missing = [];
for n = 1:length(subs_Matt)
    N = find(S(:,1)==subs_Matt(n));
    if N; S(N,2) = 1; else subs_missing = [subs_missing; n]; end
end
%JBgrotU(S(:,2)==0,:) = [];
%JBgrotV(S(:,2)==0,:) = [];
%UV = mean([JBgrotU(:,1) JBgrotV(:,1)],2);
varsd(S(:,2)==0,:) = [];
subs_Matt(subs_missing) = []; sex(subs_missing) = [];
Fractional_area_L(subs_missing,:) = []; Fractional_area_R(subs_missing,:) = [];
Total_area_L(subs_missing,:) = []; Total_area_R(subs_missing,:) = [];

% Load CCA results
load('Results/CCA_MMP_PFM.mat','Uall')
U = mean(Uall(:,[2 6]),2);
corr([U Uall(:,[2 6])])

% Create input for CCA
S = '/home/fs0/janineb/scratch/HCP/CCA/';
Nkeep = 100;
Nperm = 100000;
load(fullfile(S,'files','Permutation_100000_MattData.mat'),'conf'); 
conf = demean(conf);
Pconf = pinv(conf);
T = [Fractional_area_R Fractional_area_L];
T = demean(T);    
T = demean(T-conf*(pinv(conf)*T));
[Glasser_fractional_area,~,~]=nets_svds(T,Nkeep); 
save('Results/input_CCA_Glasser_fractional_area.mat','Glasser_fractional_area');

% Correlate fractional area against CCA weights (see Emails David vE 6 August 2017)
[RL,pL] = corr(U,Fractional_area_L);
[RR,pR] = corr(U,Fractional_area_R);
[~,~,p_fdr] = fdr([pL pR]); pL_fdr = p_fdr(1:180); pR_fdr = p_fdr(181:end);

% regress out sex and then do again
U = demean(U); sex = demean(sex);
Fractional_area_L = demean(Fractional_area_L); Fractional_area_R = demean(Fractional_area_R);
U = U-sex*(pinv(sex)*U);
Fractional_area_L = Fractional_area_L-sex*(pinv(sex)*Fractional_area_L);
Fractional_area_R = Fractional_area_R-sex*(pinv(sex)*Fractional_area_R);
[RL_nosex,pL_nosex] = corr(U,Fractional_area_L);
[RR_nosex,pR_nosex] = corr(U,Fractional_area_R);
[~,~,p_nosex_fdr] = fdr([pL_nosex pR_nosex]); pL_nosex_fdr = p_nosex_fdr(1:180); pR_nosex_fdr = p_nosex_fdr(181:end);

% Report some results
A = setdiff(find(p_fdr<0.05),find(p_nosex_fdr<0.05));
fprintf('Regions only significant if not controlling for gender: \n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',Names360{A});

A = setdiff(find(p_nosex_fdr<0.05),find(p_fdr<0.05));
fprintf('Regions only significant after controlling for gender: \n%s\n%s\n',Names360{A});

fprintf('\nLooking further into results after controlling for gender: \n\n');

A = intersect(find(pR_nosex_fdr<0.05)', find(pL_nosex_fdr<0.05)');
fprintf('Regions with bilateral FDR-corrected significance: \n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',Names_bilateral{A});
FDRout = Names_bilateral(A);
FDRout(1:length(A),2) = num2cell(RL_nosex(A));
FDRout(1:length(A),3) = num2cell(RR_nosex(A));

A = setdiff(find(pL_nosex_fdr<0.05)', find(pR_nosex_fdr<0.05)'); 
fprintf('Regions with unilateral FDR-corrected significance: \n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',Names360{A});
FDRout = [FDRout; [Names360(A) num2cell(RL_nosex(A)') num2cell(RR_nosex(A)')]];
A = setdiff(find(pR_nosex_fdr<0.05)', find(pL_nosex_fdr<0.05)'); A = A+180;
fprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n',Names360{A});
FDRout = [FDRout; [Names360(A) num2cell(RL_nosex(A-180)') num2cell(RR_nosex(A-180)')]];

clear vars varsd Steve_names
save('Results/Fractional_Area_results.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

example = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
dt_isnan = isnan(example.dtseries(:,1));

MG = ft_read_cifti('/vols/Data/HCP/workbench/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group.32k_fs_LR.dlabel.nii');
MG = MG.indexmax;
MG(dt_isnan(1:size(MG,1))==1) = [];
MGnew = zeros(91282,3);
for n = 1:180
    MGnew(MG==n,1) = RR_nosex(n);
    if pR_nosex_fdr(n)<0.05
        MGnew(MG==n,2) = RR_nosex(n);
        MGnew(MG==n,3) = 1;
    end
    MGnew(MG==n+180,1) = RL_nosex(n);
    if pL_nosex_fdr(n)<0.05
        MGnew(MG==n+180,2) = RL_nosex(n);
        MGnew(MG==n+180,3) = 1;
    end     
end
example.dtseries = repmat(example.dtseries(:,1),1,3); example.time = 1:3; example.hdr.dim(6) = 3;
example.dtseries(isnan(example.dtseries(:,1))==0,:) = MGnew;
ft_write_cifti('Results/Fractional_Area_results',example,'parameter','dtseries')


