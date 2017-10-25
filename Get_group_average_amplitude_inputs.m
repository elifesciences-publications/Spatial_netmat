clear all; close all; clc

% Set paths:
path(path,fullfile('/vols/Scratch/janineb','matlab','cifti-matlab'))
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
addpath ~steve/NETWORKS/FSLNets;
path(path,'/home/fs0/janineb/scratch/matlab')
PFMdir = '/home/fs0/samh/scratch/PROFUMO/S820_M50_Aug16.pfm/FinalModel/';

ts = PFM_loadTimeCourses_pfmNEW(PFMdir,0.72,1,1,0,[]);

% Create simulated data and perform weighted regression of ICA maps onto this
subs = dir('/vols/Scratch/janineb/PROFUMO/PFM50_MSMall900/Model7_S.pfm/Subjects/*');
subs = subs(3:end);
runs = {'1_LR','1_RL','2_LR','2_RL'};
sStd = zeros(50,820*4);
noiseStd = zeros(1,820*4);
H = zeros(50,820*4);
norm = zeros(91282,820*4);
I = 1;
for s = 1:ts.Nsubjects
    fprintf('running subject %d \n', s)
    % Get subject timeseries
    TCall = ts.ts((s-1)*ts.NtimepointsPerSubject+1:s*ts.NtimepointsPerSubject,:);   
    for x = 1:4;
        % Get run timeseries
        TC = TCall((x-1)*1200+1:x*1200,:);
        
        % Get subject correlations
        sCov = cov(TC);
        sStd(:,I) = sqrt(diag(sCov));
        
        % Load posterior nosie variance and amplitude weights needed for creating full dataset
        a = importfile(fullfile(PFMdir,'Subjects', subs(s).name,'Runs',runs{x},'NoisePrecision.post','GammaPosterior.txt'));
        noiseStd(I) = sqrt(a(2)/a(1)); 
        H(:,I) = h5read(fullfile(PFMdir,'Subjects',subs(s).name,'Runs',runs{x},'ComponentWeightings.post','Means.hdf5'),'/dataset');
        
        % Undo variance normalisation that PFM pipeline applies (to avoid issues with relative strength of subcortical modes)
        norm(:,I) = h5read(fullfile(PFMdir,'..','Preprocessing',subs(s).name,sprintf('%s_Normalisation.hdf5',runs{x})),'/dataset');
        I = I+1;
    end
end

sStd_group = mean(sStd,2);
noiseStd_group = mean(noiseStd);
H_group = mean(H,2);
norm_group = mean(norm,2);

save('group_amplitude_inputs.mat','sStd_group','noiseStd_group','H_group','norm_group');


