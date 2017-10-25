clear all; close all; clc

path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')
addpath ~steve/NETWORKS/FSLNets;

% Load PFM group maps
PFMmaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
PFMmaps = PFMmaps.dtseries; PFMmaps(isnan(PFMmaps(:,1))==1,:) = [];

% Load PFM subject maps
PFMdir = '/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/';
Subject_maps = PFM_loadSubjectSpatialMaps(PFMdir,1:50);

% Threshold
map_thresh = 1;
Mfixed = zeros(size(Subject_maps));
Mfixed(Subject_maps<-map_thresh) = Subject_maps(Subject_maps<-map_thresh);
Mfixed(Subject_maps>map_thresh) = Subject_maps(Subject_maps>map_thresh);
Mfixed(Mfixed<0) = -1;
Mfixed(Mfixed>0) = 1;
    
map_thresh = -95;
Mprct = zeros(size(Subject_maps));
Uthr = prctile(shiftdim(Subject_maps,1),-map_thresh);
Lthr = prctile(shiftdim(Subject_maps,1),100+map_thresh);
Mprct(Subject_maps<shiftdim(repmat(Lthr,size(Subject_maps,2),1,1),2)) = Subject_maps(Subject_maps<shiftdim(repmat(Lthr,size(Subject_maps,2),1,1),2));
Mprct(Subject_maps>shiftdim(repmat(Uthr,size(Subject_maps,2),1,1),2)) = Subject_maps(Subject_maps>shiftdim(repmat(Uthr,size(Subject_maps,2),1,1),2));
Mprct(Mprct<0) = -1;
Mprct(Mprct>0) = 1;

% Compare overlap and rank subjects
DMN = sum(squeeze(Mfixed(:,:,14)).*repmat(PFMmaps(:,14)',820,1),2);
[~,i] = sort(DMN,'ascend');
cifti_example = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
cifti_example.dtseries = repmat(cifti_example.dtseries(:,1),1,20); cifti_example.hdr.dim(6) = 20; cifti_example.time = 1:20;
cifti_example.dtseries(isnan(cifti_example.dtseries(:,1))==0,1) = PFMmaps(:,14);
cifti_example.dtseries(isnan(cifti_example.dtseries(:,1))==0,2:end) = squeeze(Mfixed(i(13:43:820),:,14))';
ft_write_cifti('test',cifti_example,'parameter','dtseries');
