% Set paths
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/vols/Scratch/janineb/HCP/DMN/DMN_functions')
P = '/vols/Scratch/HCP/tfMRI/Q900/';
MSM = '_MSMAll';
addpath('~steve/matlab/icasso122','~steve/matlab/FastICA_25'); 

% Load task maps
load('/vols/Scratch/janineb/HCP/DMN/CompareMethods/CCA/CCA_input_data/Task_subject_spatial_correlation.mat','Names'); Names = Names;
subs = dir('/vols/Scratch/janineb/PROFUMO/PFM50_MSMall900/Model7_S.pfm/Subjects/');
subs = subs(3:end);
Task_maps = nan(91282,86,size(subs,1));
subs_missing = [];
for i = 1:size(subs,1)
        if exist(fullfile(P,subs(i).name,'MNINonLinear','Results',...
                sprintf('%s_All_task_maps_%s.dtseries.nii',MSM(2:end),subs(i).name)),'file');
            maps = ft_read_cifti(fullfile(P,subs(i).name,'MNINonLinear','Results',...
                sprintf('%s_All_task_maps_%s.dtseries.nii',MSM(2:end),subs(i).name)));
            maps = maps.dtseries; maps(isnan(maps(:,1))==1,:) = [];
            Task_maps(:,:,i) = maps;
        else
            subs_missing = [subs_missing; i];
        end
end
subs_missing = sort([subs_missing; 169]); 
Task_maps(:,:,subs_missing) = []; clear i maps

% Average task data
Task_maps_mean = mean(Task_maps,3);

% Run ICA
Ncomponents = 15;
nonlin='tanh'; % pow3 / tanh / gauss / skew    % tanh seems best overall at least for multi-subjects
approach='symm'; % symm / defl
[icaS,icaA,icaW] = fastica(Task_maps_mean','approach',approach,'g',nonlin,'epsilon',1e-11,'maxNumIterations',3000,'lastEig',Ncomponents);
Task_ICA = icaW*Task_maps_mean';
M = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
M.time = 1:Ncomponents; M.hdr.dim(6) = Ncomponents; M.dtseries = repmat(M.dtseries(:,1),1,Ncomponents);
M.dtseries(isnan(M.dtseries(:,1))==0,:) = icaS'; %Task_ICA';
ft_write_cifti('Task_ICA15',M,'parameter','dtseries');

% Run PCA
[u,s,v] = svd(Task_maps_mean);
s = s(s>0);
figure; plot(s)

% Save results
save(sprintf('Results/task_ICA_%d.mat',Ncomponents),'icaA','icaS','icaW','u','s','v')
