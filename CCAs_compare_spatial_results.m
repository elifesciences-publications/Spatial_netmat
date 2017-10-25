clear all; close all; clc

% Set and load inputs:
Dir = '/home/fs0/janineb/scratch/';
Nperm = 0;

% Set up paths
path(path,fullfile(Dir,'matlab'))
path(path,fullfile(Dir,'HCP','CCA','files'))
path(path,fullfile(Dir,'matlab','cifti-matlab'))
path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
addpath ~steve/NETWORKS/FSLNets;
path(path,fullfile(Dir,'HCP','netmat_simulations','palm-alpha106'));

% Cut down to correct subjects
load('Results/subs_Matt.mat');
PFMsubs = dir('/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/Subjects');
PFMsubs = PFMsubs(3:end); Task_subs = PFMsubs;
subs_remove = [122 160 162 169 200 248 250 260 286 292 295 312 315 320 337 353 367 425 457 463 607 669 679 684 690 691 701 760 764 797];
Task_subs(subs_remove) = []; PFMsubs(169) = [];
subs_Matt(373) = [];
Spfm = zeros(size(PFMsubs,1),2); for n = 1:size(PFMsubs,1); Spfm(n,1) = str2double(PFMsubs(n).name); end
Stask = zeros(size(Task_subs,1),2); for n = 1:size(Task_subs,1); Stask(n,1) = str2double(Task_subs(n).name); end
for n = 1:length(subs_Matt)
    N = find(Spfm(:,1)==subs_Matt(n)); if N; Spfm(N,2) = 1; end
    N = find(Stask(:,1)==subs_Matt(n)); if N; Stask(N,2) = 1; else Smmp=subs_Matt(n); end
end

% Load CCA input data
brainsAll = load('Results/input_CCA_table1.mat');
brainsAll = rmfield(brainsAll,'Glasser_amp');
brainsAll = rmfield(brainsAll,'ICA200_amp');
brainsAll = rmfield(brainsAll,'ICA25_amp');
brainsAll = rmfield(brainsAll,'PFM50_Fnetmat');
brainsAll = rmfield(brainsAll,'PFM50_Pnetmat');
brainsAll = rmfield(brainsAll,'PFM50_amp');
brainsAll = rmfield(brainsAll,'Yeo_amp');
load('Results/input_CCA_table1_fractional_area.mat');
brainsAll.area = Fractional_area;
A1 = fieldnames(brainsAll);
U = zeros(440,2*length(A1)); I = 1;
for x = 1:size(A1,1)
    fprintf('running CCA on %s\n',A1{x});
    brains = brainsAll.(A1{x});
    [~,u,~,~,~,~] = S_CCA(brains,Nperm);
    if size(u,1)==819; u(Spfm(:,2)==0,:) = []; end
    if size(u,1)==790; u(Stask(:,2)==0,:) = []; end
    if size(u,1)==441; u(373,:) = []; end
    U(:,I:I+1) = u(:,1:2); I = I+2;
end

figure
imagesc(corr(U),[-1 1])
hline(2.5:2:18,'k')
vline(2.5:2:18,'k')
set(gca,'xtick',1.5:2:18,'xticklabel',A1,'ytick',1.5:2:18,'yticklabel',A1)


