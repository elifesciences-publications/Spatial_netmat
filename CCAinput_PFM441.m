clear all; close all; clc

% Set data directory
Ddir = 'S820_M50_Aug16.pfm';
d = 50;
Nkeep = 100;

% Set paths
S = '/home/fs0/janineb/scratch/HCP/CCA/CCA_new_PFMs'; 
PFMdir = '/home/fs0/samh/scratch/HCP/';
P = fullfile(PFMdir,Ddir,'FinalModel');
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/')
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
addpath ~steve/NETWORKS/FSLNets;
load(fullfile(S,'..','files','Permutation_100000_MattData.mat'),'conf');

% Get subjects
subs = dir('/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/Subjects');
load('Results/subs_Matt.mat');
subs = subs(3:end);
S = zeros(size(subs,1),2);
for n = 1:size(subs,1); S(n,1) = str2double(subs(n).name); end; clear subs
for n = 1:length(subs_Matt)
    N = find(S(:,1)==subs_Matt(n));
    S(N,2) = 1;
end
subs_remove = find(S(:,2)==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Spatial covariance and correlation:
noise = [];
signal = setdiff(1:d,noise);

maps = PFM_loadSubjectSpatialMaps(P,signal);
maps(subs_remove,:,:) = [];
C = zeros(size(maps,1),size(maps,1));

% Prepare confounds
conf = demean(conf);
Pconf = pinv(conf);

% For each mode deconfound and create features input
for i = 1:size(maps,3)
    fprintf('processing map %d (out of %d)\n',i,size(maps,3));
    m = maps(:,:,i)';
    m = demean(m,2);
    m = (m'-conf*(Pconf*m'))';
    
    C = C + m' * m;
end

% Save input data for CCA 
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
PFM441_spatial = V(:,1:Nkeep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Doing netmats \n')
ts = PFM_loadTimeCourses_pfmNEW(P,0.72,1,1,0,[]);
Pnetmat_tsPFMhf = nets_netmats(ts,1,'ridgep',0.01); Pnetmat_tsPFMhf(subs_remove,:) = [];
Fnetmat_tsPFMhf = nets_netmats(ts,1,'corr'); Fnetmat_tsPFMhf(subs_remove,:) = [];
% Remove confounds
Pnetmat_tsPFMhf = demean(Pnetmat_tsPFMhf);    
Pnetmat_tsPFMhf = demean(Pnetmat_tsPFMhf-conf*(pinv(conf)*Pnetmat_tsPFMhf));  
Fnetmat_tsPFMhf = demean(Fnetmat_tsPFMhf);    
Fnetmat_tsPFMhf = demean(Fnetmat_tsPFMhf-conf*(pinv(conf)*Fnetmat_tsPFMhf));  
% Save input data for CCA
[PFM441_Pnetmat,~,~]=nets_svds(Pnetmat_tsPFMhf,Nkeep); 
[PFM441_Fnetmat,~,~]=nets_svds(Fnetmat_tsPFMhf,Nkeep); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Doing amplitude weights \n')
runs = {'LR','RL'};
subs = dir(fullfile(P,'Subjects'));
A = zeros(820,d*4);
for i = 1:size(subs,1)-2
   I = 1;
    for day = 1:2
        for r = 1:2
            A(i,(I-1)*d+1:I*d) = h5read(fullfile(P,'Subjects',subs(i+2).name,'Runs',sprintf('%d_%s',day,runs{r}),'ComponentWeightings.post','Means.hdf5'),'/dataset');
            I = I+1;
        end
    end
end
A(subs_remove,:) = [];
% Remove confounds
A = demean(A);    
A = demean(A-conf*(pinv(conf)*A));  
% Save input data for CCA
[PFM441_amplitude,~,~]=nets_svds(A,Nkeep);  

save('Results/input_CCA_table1_PFM_441.mat','PFM441_amplitude','PFM441_Fnetmat','PFM441_Pnetmat','PFM441_spatial','-v7.3')


